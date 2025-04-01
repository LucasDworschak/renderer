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

#pragma once

#include <QDateTime>
#include <QObject>
#include <QTimer>

#include <nucleus/camera/Definition.h>
#include <nucleus/timing/TimerManager.h>

namespace nucleus::utils {

/**
 * Currently the available benchmarks to choose from are hardcoded in StatsWindow.qml
 * -> you have to declare it there and create it in RenderingContext
 */
class Benchmark : public QObject {
    Q_OBJECT
public:
    Benchmark(QString name, unsigned id, std::vector<std::string> positions);
public slots:
    void activate(unsigned id);
    void increase_load_count();
    void decrease_load_count();

    void receive_measurements(QList<nucleus::timing::TimerReport> values);

private slots:
    void next_place();
    void update_record_state();

signals:
    void camera_definition_set_by_user(const nucleus::camera::Definition&);

private:
    const QString m_name;
    const unsigned m_id;

    int m_load_count;
    bool m_record_data;
    bool m_activated;
    QDateTime m_previous_time;

    size_t m_current_position;
    std::vector<std::string> m_positions;

    std::unique_ptr<QTimer> m_continue_timer;
    std::unique_ptr<QTimer> m_buffer_timer;

    std::unordered_map<std::string, std::vector<QList<nucleus::timing::TimerReport>>> m_timings;
    std::unordered_map<std::string, float> m_total_time_per_location;

    void create_report();

    static constexpr int time_per_location = 10000;
    static constexpr int buffer_time = 5000; // we need a sufficiently high time to make sure that no other tiles arrive anymore
};

} // namespace nucleus::utils
