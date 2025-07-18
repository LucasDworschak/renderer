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
    uint16_t outline_width; // line-width
    uint16_t outline_dash; // line-dasharray

    bool operator==(const LayerStyle& other) const = default;
};

struct LowestEncounteredZoom {
    unsigned lowest_zoom;
    unsigned max_zoom;
};

struct StyleHasher {
    size_t operator()(const LayerStyle& style) const
    {
        size_t seed = 0;

        radix::hasher::hash_combine<uint32_t>(seed, style.fill_color);
        radix::hasher::hash_combine<uint32_t>(seed, style.outline_color);
        radix::hasher::hash_combine<float>(seed, style.outline_width);
        radix::hasher::hash_combine<uint32_t>(seed, style.outline_dash);

        return seed;
    }

    size_t operator()(const std::pair<std::pair<std::string, int>, std::shared_ptr<StyleExpressionBase>>& pair) const
    {
        size_t seed = 0;

        radix::hasher::hash_combine<std::string>(seed, pair.first.first);
        radix::hasher::hash_combine<int>(seed, pair.first.second);
        if (pair.second != nullptr) // if there is no filter -> this would be null
            radix::hasher::hash_combine<uint32_t>(seed, pair.second->hash());

        return seed;
    }

    bool operator()(const std::pair<std::pair<std::string, int>, std::shared_ptr<StyleExpressionBase>>& lhs,
        const std::pair<std::pair<std::string, int>, std::shared_ptr<StyleExpressionBase>>& rhs) const
    {
        // in order to have consistent insertion order, we are using first the string compare than the hashes of the styleexpressions
        auto str_comp = lhs.first.first.compare(rhs.first.first);
        if (str_comp == 0) {
            if (lhs.first.second != rhs.first.second) {
                return lhs.first.second < rhs.first.second;
            }

            auto h1 = (lhs.second == nullptr) ? 0 : lhs.second->hash();
            auto h2 = (rhs.second == nullptr) ? 0 : rhs.second->hash();

            return h1 < h2;
        }

        return str_comp < 0;
    }

    size_t operator()(const std::pair<std::string, int>& pair) const
    {
        size_t seed = 0;

        radix::hasher::hash_combine<std::string>(seed, pair.first);
        radix::hasher::hash_combine<int>(seed, pair.second);

        return seed;
    }
};

class Style : public QObject {
    Q_OBJECT
public:
    Style(const QString& filename);
    Style(Style&& style);

    void parse_colors(const QJsonValue& value, std::vector<std::pair<uint8_t, uint32_t>>& colors, float& base);
    void parse_opacities(const QJsonValue& value, std::vector<std::pair<uint8_t, uint8_t>>& opacities, float& base);
    void parse_dashes(const QJsonValue& value, std::vector<std::pair<uint8_t, std::pair<uint8_t, uint8_t>>>& dashes, float& base);
    void parse_line_widths(const QJsonValue& value, std::vector<std::pair<uint8_t, uint16_t>>& widths, float& base);
    uint32_t parse_color(const QJsonValue& value);
    uint8_t parse_opacity(const QJsonValue& value);
    std::pair<uint8_t, uint8_t> parse_dash(const QJsonValue& value);
    uint16_t parse_line_width(const QJsonValue& value);

    static uint32_t premultiply_alpha(uint32_t color);

    std::vector<StyleLayerIndex> indices(std::string layer_name,
        int type,
        unsigned zoom,
        const mapbox::vector_tile::feature& feature,
        std::array<int, constants::max_style_expression_keys>* temp_values);

    std::shared_ptr<const nucleus::Raster<glm::u32vec4>> styles() const;
    std::shared_ptr<const nucleus::Raster<glm::u32vec4>> visible_styles() const;

    bool update_visible_styles();

public slots:
    void load();

signals:
    void load_finished(std::shared_ptr<const nucleus::Raster<glm::u32vec4>> styles);

private:
    // sets alpha to 0 on a 32 bit color
    static constexpr uint remove_alpha_mask = 4294967040u;

    float interpolation_factor(uint8_t zoom, float base, uint8_t zoom1, uint8_t zoom2);
    uint32_t interpolate_color(float t, uint32_t color1, uint32_t color2);

    float stringToFloat(const std::string& value);
    float rgb2linear(uint8_t channel);
    uint8_t linear2rgb(float linear);

    // contains the style info that was parsed from the stylesheet
    std::shared_ptr<const nucleus::Raster<glm::u32vec4>> m_styles;
    // only contains that were encountered by the call of the indices function
    std::shared_ptr<nucleus::Raster<glm::u32vec4>> m_visible_styles;

    std::unordered_map<std::pair<std::string, int>, StyleFilter, StyleHasher> m_layer_to_style;

    // styles on the server and in the stylesheet might be a bit different
    // we want to ensure that we fade to alpha 0 if the next lower tile zoom does not contain a blendable style
    // we ensure that a feature that uses style blending (!all_styles_same) has a style at constants::style_zoom_range.y-1
    // key: style_index of highest style
    // value: lowest zoom_level of current style
    std::unordered_map<uint32_t, LowestEncounteredZoom> m_lowest_encountered_zoom;
    std::vector<StyleLayerIndex> m_styles_to_update;

    QString m_filename;
};

} // namespace nucleus::vector_layer
