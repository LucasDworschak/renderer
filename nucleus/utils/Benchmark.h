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

#include <QObject>
#include <nucleus/timing/TimerManager.h>

namespace nucleus::utils {

class Benchmark : public QObject {
    Q_OBJECT
public:
    Benchmark(QString name);
public slots:
    void start_benchmark();
    void increase_load_count();
    void decrease_load_count();

    void receive_measurements(QList<nucleus::timing::TimerReport> values);

signals:
    void benchmark_finished();

private:
    QString m_name;
    int m_wait_count;
    bool m_benchmark_ready;

    void update_benchmark_state();
};

} // namespace nucleus::utils
