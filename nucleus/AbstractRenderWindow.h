/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2023 Adam Celarek
 * Copyright (C) 2024 Gerald Kimmersdorfer
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
#include <unordered_set>

#include <glm/glm.hpp>

#include "nucleus/tile/types.h"
#include "utils/ColourTexture.h"

class QOpenGLFramebufferObject;

namespace nucleus {
namespace tile::utils {
    class AabbDecorator;
    using AabbDecoratorPtr = std::shared_ptr<AabbDecorator>;
}
namespace camera {
    class Definition;
    class AbstractDepthTester;
}

class AbstractRenderWindow : public QObject {
    Q_OBJECT
public:
    virtual void initialise_gpu() = 0;
    virtual void resize_framebuffer(int width, int height) = 0;
    virtual void paint(QOpenGLFramebufferObject* framebuffer = nullptr) = 0;
    virtual void destroy() = 0;
    [[nodiscard]] virtual camera::AbstractDepthTester* depth_tester() = 0;
    [[nodiscard]] virtual utils::ColourTexture::Format ortho_tile_compression_algorithm() const = 0;

public slots:
    virtual void update_camera(const camera::Definition& new_definition) = 0;
    virtual void update_debug_scheduler_stats(const QString& stats) = 0;
    virtual void pick_value(const glm::dvec2& screen_space_coordinates) = 0;

signals:
    void update_requested();
    void value_picked(uint32_t value);
};

}
