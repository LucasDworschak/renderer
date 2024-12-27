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

#include "StyleFilter.h"

#include <mapbox/vector_tile.hpp>

namespace nucleus::vector_layer {

void StyleFilter::add_filter(uint32_t style_index, std::shared_ptr<StyleExpressionBase> filter, glm::uvec2 zoom_range)
{
    for (unsigned i = zoom_range.x; i < zoom_range.y; i++) {
        if (!m_filter.contains(i)) {
            m_filter[i] = std::vector<std::pair<uint32_t, std::shared_ptr<StyleExpressionBase>>>();
        }
        m_filter[i].push_back(std::make_pair(style_index, filter));
    }
}

uint32_t StyleFilter::style_index(unsigned zoom, const mapbox::vector_tile::feature& feature) const
{
    if (!m_filter.contains(zoom)) {
        // qDebug() << "filter at zoom not found ";
        return -1u; // not found
    }

    for (const auto& pair : m_filter.at(zoom)) {
        if (pair.second == nullptr)
            return pair.first; // no filter is here -> we assume that every feature with layername and zoom is valid
        if (pair.second->matches(feature)) {
            return pair.first; // first filter that matches returns the style index
        }
    }

    // qDebug() << "filter and zoom exist but no match found";
    return -1u; // not found
}

} // namespace nucleus::vector_layer
