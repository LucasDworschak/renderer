/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Lucas Dworschak
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

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

#include <glm/glm.hpp>

#include "StyleExpression.h"

namespace mapbox::vector_tile {
class feature;
}

namespace nucleus::vector_layer {

class StyleFilter {
public:
    StyleFilter() { }

    void add_filter(uint32_t style_index, uint32_t layer_index, std::shared_ptr<StyleExpressionBase> filter, glm::uvec2 zoom_range);

    std::pair<uint32_t, uint32_t> indices(unsigned zoom, const mapbox::vector_tile::feature& feature) const;

private:
    // zoom level -> vector<style_index,StyleExpression>
    std::unordered_map<unsigned, std::vector<std::tuple<uint32_t, uint32_t, std::shared_ptr<StyleExpressionBase>>>> m_filter;
};

} // namespace nucleus::vector_layer
