/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2022 Adam Celarek
 * Copyright (C) 2023 Jakob Lindner
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

#include "LinearCameraAnimation.h"

#include <QEasingCurve>

#include "AbstractDepthTester.h"

namespace nucleus::camera {

LinearCameraAnimation::LinearCameraAnimation(Definition start, Definition end)
    : m_start(start.model_matrix())
    , m_end(end.model_matrix())
{
    m_current_duration = 0;
    m_stopwatch.restart();
}

std::optional<Definition> LinearCameraAnimation::update(Definition camera, AbstractDepthTester*)
{
    if (m_current_duration >= m_total_duration) {
        return {};
    }

    auto dt = m_stopwatch.lap().count();
    if (m_current_duration + dt > m_total_duration) { // last step
        dt = m_total_duration - m_current_duration;
    }

    const float t_eased = ease_in_out((m_current_duration + float(dt)) / float(m_total_duration));

    m_current_duration += dt;

    const auto mix_factor = t_eased; //m_current_duration / float(m_total_duration);
    const auto new_matrix = m_start * double(1 - mix_factor) + m_end * double(mix_factor);
    camera.set_model_matrix(new_matrix);

    return camera;
}

float LinearCameraAnimation::ease_in_out(float t)
{
    QEasingCurve c(QEasingCurve::Type::OutExpo);
    return c.valueForProgress(t);
}
}
