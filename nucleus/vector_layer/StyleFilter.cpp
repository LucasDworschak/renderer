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

void StyleFilter::add_filter(uint32_t style_index, uint32_t layer_index, std::shared_ptr<StyleExpressionBase> filter, glm::uvec2 zoom_range)
{
    for (unsigned i = zoom_range.x; i < zoom_range.y + 1; i++) {
        if (!m_filter.contains(i)) {
            m_filter[i] = std::vector<std::tuple<uint32_t, uint32_t, std::shared_ptr<StyleExpressionBase>>>();
        }
        m_filter[i].push_back(std::make_tuple(style_index, layer_index, filter));
    }
}

std::vector<std::pair<uint32_t, uint32_t>> StyleFilter::indices(unsigned zoom, const mapbox::vector_tile::feature& feature) const
{
    if (!m_filter.contains(zoom)) {
        // qDebug() << "filter at zoom not found ";
        return {}; // not found
    }

    const auto type = feature.getType();
    const auto properties = feature.getProperties();

    auto styles = std::vector<std::pair<uint32_t, uint32_t>>();

    for (const auto& values : m_filter.at(zoom)) {
        // qDebug() << (std::get<2>(values) == nullptr);
        if (std::get<2>(values) == nullptr) // no filter is here -> we assume that every feature with layername and zoom is valid
        {
            // if (m_filter.at(zoom).size() > 1) {
            //     qDebug() << "more than one filter is applicable: " << m_filter.at(zoom).size();
            //     qDebug() << std::get<0>(values);
            //     qDebug() << std::get<1>(values);
            // }
            // assert(m_filter.at(zoom).size() == 1); // if there is no filter -> there should only be one valid value
            styles.push_back(std::make_pair(std::get<0>(values), std::get<1>(values)));
            // if (m_filter.at(zoom).size() > 1) {
            //     qDebug() << "why not working?: " << m_filter.at(zoom).size();
            // }
        } else if (std::get<2>(values)->matches(type, properties)) {
            styles.push_back(std::make_pair(std::get<0>(values), std::get<1>(values)));
        }
    }

    // qDebug() << "return style!";

    return styles;
}

} // namespace nucleus::vector_layer
