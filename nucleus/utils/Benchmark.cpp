/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2025 Lucas Dworschak
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "Benchmark.h"

#include <nucleus/camera/PositionStorage.h>
#include <nucleus/version.h>

#include <QDesktopServices>
#include <QFile>
#include <QStandardPaths>
#include <QString>
#include <QTimer>
#include <QUrl>

#include <filesystem>
#include <set>

namespace nucleus::utils {

Benchmark::Benchmark(QString name, unsigned id, std::vector<std::string> positions)
    : m_name(name)
    , m_id(id)
    , m_load_count(0)
    , m_record_data(false)
    , m_activated(false)
    , m_current_position(0)
    , m_positions(positions)
{
    qDebug() << "Benchmark: created" << name << id;

    m_continue_timer = std::make_unique<QTimer>(this);
    m_continue_timer->setSingleShot(true);
    m_warm_up_timer = std::make_unique<QTimer>(this);
    m_warm_up_timer->setSingleShot(true);

    connect(m_continue_timer.get(), &QTimer::timeout, this, &Benchmark::next_place);
    connect(m_warm_up_timer.get(), &QTimer::timeout, this, &Benchmark::update_record_state);
}

void Benchmark::register_scheduler(QString name)
{
    m_is_finished[name] = false;
    m_has_quads_requested[name] = true;
}

void Benchmark::activate(unsigned id)
{
    // if there is a timer running -> stop it
    m_continue_timer->stop();
    m_warm_up_timer->stop();

    // reset any data
    m_current_position = 0;
    m_record_data = false;

    if (id != m_id) { // -> another benchmark was selected
        m_activated = false;
        return;
    }

    m_activated = true;

    qDebug() << "Benchmark: starting" << m_name << m_id;

    m_continue_timer->start(100);
}

void Benchmark::increase_load_count(const QString& scheduler_name)
{
    m_is_finished[scheduler_name] = false;
    m_load_count += 1;
    m_record_data = false;

    if (m_continue_timer->isActive()) {
        // new tiles arrived while recording -> stop the timer and remove all data fromt he current location
        m_continue_timer->stop();

        m_timings[m_positions[m_current_position - 1]].clear();
    }

    if (m_warm_up_timer->isActive()) {
        // new data is being loaded -> we have to wait a bit longer
        m_warm_up_timer->stop();
    }
}

void Benchmark::decrease_load_count(const QString& scheduler_name)
{
    m_is_finished[scheduler_name] = true;
    m_load_count -= 1;

    start_if_finished_loading();
}

void Benchmark::start_if_finished_loading()
{
    if (!m_activated || m_load_count != 0 || m_warm_up_timer->isActive())
        return; // not finished yet

    for (const auto& entry : m_is_finished) {
        if (!m_is_finished[entry.first] || m_has_quads_requested[entry.first]) {
            // not finished yet
            return;
        }
    }

    qDebug() << "benchmark: " << m_positions[m_current_position - 1];

    // all are finished -> we can finally start with the benchmark
    // start the warm-up timer to ensure that everything is on the gpu
    m_warm_up_timer->start(warm_up_time);
}

void Benchmark::check_requested_quads(const QString& scheduler_name, const QVariantMap& new_stats)
{
    m_has_quads_requested[scheduler_name] = new_stats["n_quads_requested"] != 0;
}

void Benchmark::next_place()
{
    m_record_data = false;

    // reset dicts
    for (const auto& entry : m_is_finished) {
        m_is_finished[entry.first] = false;
        m_has_quads_requested[entry.first] = true;
    }

    m_previous_time = QDateTime::currentDateTime();

    if (m_current_position >= m_positions.size()) {
        // no more locations to benchmark -> we are finished
        m_warm_up_timer->stop();

        m_activated = false;

        create_report();

        return;
    }

    // change to next location
    emit camera_definition_set_by_user(nucleus::camera::PositionStorage::instance()->get(m_positions[m_current_position++]));
}

void Benchmark::create_report()
{

    auto comp = [](const QString& a, const QString& b) {
        // make sure that gpu_total and cpu_total are the first two values
        if (a != b) {
            if (a == "gpu_total")
                return true;
            if (b == "gpu_total")
                return false;
            if (a == "cpu_total")
                return true;
            if (b == "cpu_total")
                return false;
        }

        return a < b;
    };

    std::set<QString, decltype(comp)> metrics;

    // <location, <timer_name, values>>
    std::unordered_map<std::string, std::map<QString, std::vector<float>>> reports;

    for (const auto& [location, timer_reports] : m_timings) {
        reports[location] = std::map<QString, std::vector<float>>();

        for (const auto& report : timer_reports) {
            for (const auto& entry : report) {
                reports[location][entry.name].push_back(entry.value);
                metrics.insert(entry.name);
            }
        }
    }

    // open file for writing
    // the report will be written in a similar directory as the tile cache -> linux: /home/<user>/.cache/AlpineMaps.org/AlpineApp/benchmarks/
    const auto base_path = std::filesystem::path(QStandardPaths::writableLocation(QStandardPaths::CacheLocation).toStdString());
    std::filesystem::create_directories(base_path);
    const auto path = base_path / "benchmarks";
    std::filesystem::create_directories(path);

    const auto datetime = QDateTime::currentDateTime();
    QFile file(path / (m_name.toStdString() + "_" + nucleus::version() + "_" + datetime.toString("hh-mm-ss").toStdString() + ".csv"));
    const auto success = file.open(QIODeviceBase::WriteOnly);
    Q_UNUSED(success);
    assert(success);
    QTextStream out_stream(&file);

    out_stream << datetime.toString("yyyy-MM-dd-hh-mm-ss") << "\t" << QString::fromStdString(nucleus::version()) << "\n";
    out_stream << "Location\t";
    for (const auto& metric : metrics) {
        out_stream << metric << "\t avg \t min \t max \t";
    }

    out_stream << "\n";

    // map<metric, map<location, avg>>
    std::unordered_map<QString, std::unordered_map<std::string, float>> avg_values_per_location;

    for (const auto& [location, timers] : reports) {
        out_stream << QString::fromStdString(location) << "\t";

        for (const auto& metric : metrics) {
            float avg = 0.0;
            float min = 999999.9;
            float max = 0.0;

            const auto& values = timers.at(metric);

            for (const auto& value : values) {
                if (value < min)
                    min = value;
                if (value > max)
                    max = value;
                avg += value;
            }
            avg /= values.size();

            avg_values_per_location[metric][location] = avg;

            out_stream << "\t" << avg << "\t" << min << "\t" << max << "\t";
        }
        out_stream << "\n";
    }
    out_stream.flush(); // Ensure all data is written to the file.
    file.close(); // Close the file.

    {
        // append to cummulative file
        QFile file(path / (m_name.toStdString() + "_combined" + ".csv"));
        const auto is_new_file = !file.exists();
        const auto success = file.open(QIODeviceBase::WriteOnly | QIODevice::Append);
        Q_UNUSED(success);
        assert(success);
        QTextStream out_stream(&file);

        if (is_new_file) {
            // create headers for the combined file
            out_stream << "\t";
            for (const auto& metric : metrics) {
                bool first = true;
                for (size_t i = 0; i < avg_values_per_location[metric].size(); i++) {
                    if (first) {
                        first = false;
                        out_stream << "avg " << metric << "\t";
                    } else {
                        out_stream << "\t";
                    }
                }
            }
            out_stream << "\nversion\t";
            for (const auto& metric : metrics) {
                for (const auto& [location, avg] : avg_values_per_location[metric]) {
                    out_stream << QString::fromStdString(location) << "\t";
                }
            }
        }

        out_stream << "\n" << QString::fromStdString(nucleus::version()) << "\t";
        for (const auto& metric : metrics) {
            for (const auto& [location, avg] : avg_values_per_location[metric]) {
                out_stream << avg << "\t";
            }
        }

        QDesktopServices::openUrl(QUrl::fromLocalFile(file.fileName())); // Open with system handler.
    }

    qDebug() << "Benchmark: finished" << m_name << m_id;
}

/**
 * We are calling this method with a timer
 * if the timer was finished and we are still at wait_count 0 we can finally benchmark
 */
void Benchmark::update_record_state()
{
    if (m_load_count == 0) {
        // we are now recording -> record for the current location -> wait and then go to the next place
        // qDebug() << "benchmark: start recording";
        m_record_data = true;
        m_continue_timer->start(time_per_location);
    }
}

// void Benchmark::receive_measurements(QList<nucleus::timing::TimerReport> values)
void Benchmark::receive_measurements(QList<nucleus::timing::TimerReport> values)
{
    if (!m_activated)
        return;

    if (m_record_data && m_load_count == 0) {
        m_timings[m_positions[m_current_position - 1]].push_back(values);
    }
}

} // namespace nucleus::utils
