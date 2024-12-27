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

#include <QObject>
#include <QString>

#include "StyleFilter.h"
#include "nucleus/Raster.h"

#include <radix/hasher.h>

namespace nucleus::vector_layer {

struct Layer_Style {
    uint32_t fill_color; // fill-color or line-color
    uint32_t outline_color; // fill-outline-color or
    float outline_width; // line-width
    uint32_t outline_dash; // line-dasharray

    bool operator==(const Layer_Style& other) const = default;
};

template <class T> inline void hash_combine(std::size_t& seed, T const& v) { seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2); }

struct Hasher {
    size_t operator()(const Layer_Style& style) const
    {
        size_t seed = 0;

        hash_combine<uint32_t>(seed, style.fill_color);
        hash_combine<uint32_t>(seed, style.outline_color);
        hash_combine<float>(seed, style.outline_width);
        hash_combine<uint32_t>(seed, style.outline_dash);

        return seed;
    }
};

struct Style_Buffer_Holder {
    std::shared_ptr<const nucleus::Raster<uint32_t>> fill_styles;
    std::shared_ptr<const nucleus::Raster<uint32_t>> line_styles;
};

class Style : public QObject {
    Q_OBJECT
public:
    Style(const QString& filename);

    uint32_t parse_color(std::string value);
    uint32_t parse_dasharray(QJsonArray dash_values);

    uint32_t layer_style_index(std::string layer_name, unsigned zoom, const mapbox::vector_tile::feature& feature) const;

    Style_Buffer_Holder style_buffer() const;

public slots:
    void load();

signals:
    void load_finished();

private:
    Style_Buffer_Holder m_styles;

    // std::unordered_map<std::tuple<std::string, unsigned>, size_t, radix::hasher::for_tuple<std::string, unsigned>> m_layer_zoom_to_style;
    std::unordered_map<std::string, StyleFilter> m_layer_to_style;

    QString m_filename;
};

} // namespace nucleus::vector_layer
