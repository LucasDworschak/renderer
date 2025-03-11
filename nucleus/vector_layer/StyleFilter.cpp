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

void StyleFilter::add_filter(FilterInfo filter_info, uint8_t zoom)
{
    if (!m_filter.contains(zoom)) {
        m_filter[zoom] = std::vector<FilterInfo>();
    }
    m_filter[zoom].push_back(filter_info);
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

    for (const auto& filter_info : m_filter.at(zoom)) {
        if (filter_info.filter == nullptr) // no filter is here -> we assume that every feature with layername and zoom is valid
        {
            styles.push_back(std::make_pair(filter_info.style_index, filter_info.layer_index));

        } else if (filter_info.filter->matches(type, properties)) {
            styles.push_back(std::make_pair(filter_info.style_index, filter_info.layer_index));
        }
    }

    return styles;
}

} // namespace nucleus::vector_layer
