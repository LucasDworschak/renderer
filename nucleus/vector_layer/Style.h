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

struct LayerStyle {
    uint32_t fill_color; // fill-color or line-color
    uint32_t outline_color; // fill-outline-color or
    float outline_width; // line-width
    uint32_t outline_dash; // line-dasharray

    bool operator==(const LayerStyle& other) const = default;
};

template <class T> inline void hash_combine(std::size_t& seed, T const& v) { seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2); }

struct Hasher {
    size_t operator()(const LayerStyle& style) const
    {
        size_t seed = 0;

        hash_combine<uint32_t>(seed, style.fill_color);
        hash_combine<uint32_t>(seed, style.outline_color);
        hash_combine<float>(seed, style.outline_width);
        hash_combine<uint32_t>(seed, style.outline_dash);

        return seed;
    }
};

struct StyleBufferHolder {
    std::shared_ptr<const nucleus::Raster<glm::u32vec4>> fill_styles;
    std::shared_ptr<const nucleus::Raster<glm::u32vec4>> line_styles;
};

class Style : public QObject {
    Q_OBJECT
public:
    Style(const QString& filename);

    uint32_t parse_color(const QJsonValue& value);
    uint32_t parse_dasharray(const QJsonValue& dash_values);

    std::pair<uint32_t, uint32_t> indices(std::string layer_name, std::string type, unsigned zoom, const mapbox::vector_tile::feature& feature) const;

    StyleBufferHolder style_buffer() const;

    static QJsonArray expand(const QJsonArray& layers);

public slots:
    void load();

signals:
    void load_finished(std::shared_ptr<const nucleus::Raster<glm::u32vec4>> fill_styles, std::shared_ptr<const nucleus::Raster<glm::u32vec4>> line_styles);

private:
    StyleBufferHolder m_styles;

    uint8_t parse_opacity(const QJsonValue& value);

    QJsonValue onlyLastStopValue(const QJsonValue& value);

    static bool sub_is_array(QJsonObject obj, QString sub_key);
    static QJsonValue get_match_value(QJsonArray match_array, QString match_key);
    static std::unordered_map<QString, QJsonArray> get_sub_layer(QJsonArray filter);

    std::unordered_map<std::string, StyleFilter> m_layer_to_style;

    QString m_filename;
};

} // namespace nucleus::vector_layer
