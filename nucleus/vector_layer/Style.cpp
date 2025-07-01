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

#include "Style.h"

#include <QFile>

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

#include "nucleus/vector_layer/StyleExpander.h"
#include "nucleus/vector_layer/StyleExpression.h"
#include "nucleus/vector_layer/constants.h"

#include <cmath>
#include <regex>

namespace nucleus::vector_layer {

// https://maplibre.org/maplibre-style-spec/

Style::Style(const QString& filename)
    : m_filename(filename)
{
    // make sure that style_bits and style_buffer_size are correctly matched -> changing one means you also need to change the other
    // -1 because one bit is used to signal if it should blend with next style or not
    assert(((constants::style_buffer_size * constants::style_buffer_size) >> (constants::style_bits - 1)) == 1);
    StyleExpression::initialize();
}

Style::Style(Style&& other)
    : m_styles(std::move(other.m_styles))
    , m_layer_to_style(std::move(other.m_layer_to_style))
    , m_filename(std::move(other.m_filename))
{
}

void Style::load()
{
    QFile file(m_filename);
    const auto open = file.open(QIODeviceBase::OpenModeFlag::ReadOnly);
    assert(open);

    if (!open) {
        qDebug() << "could not open style";
        return;
    }

    const auto data = file.readAll();

    if (data.isEmpty()) {
        qDebug() << "no data";
        return;
    }

    std::vector<glm::u32vec4> style_values;

    QJsonDocument doc = QJsonDocument::fromJson(data);
    QJsonArray layers = style_expander::expand(doc.object().value("layers").toArray());

    // pair is <layer_id, filter>
    // layerid: source_layer+_+type (e.g. transportation_line)
    auto layerid_filter_to_layer_indices
        = std::map<std::pair<std::pair<std::string, int>, std::shared_ptr<StyleExpressionBase>>, std::vector<uint32_t>, StyleHasher>();
    // each entry in this vector is a different layer_index
    auto zoom_to_style = std::vector<std::map<uint8_t, LayerStyle>>();

    auto current_style_map = std::map<uint8_t, LayerStyle>();
    auto previous_layer_filter = std::pair<std::pair<std::string, int>, std::shared_ptr<StyleExpressionBase>>();
    for (const QJsonValue& obj : layers) {

        if (obj.toObject().value("type").toString() != "line" && obj.toObject().value("type").toString() != "fill") {
            continue; // not valid
        }

        if (obj.toObject().value("id").toString().startsWith("TXT"))
            continue; // we are not interested in txt

        if (!obj.toObject().contains("paint") || obj.toObject().value("paint").toObject().keys().size() == 0) {
            // no style was given
            // qDebug() << "no style: " << obj.toObject().value("id").toString().toStdString();
            continue;
        }

        // qDebug() << obj.toObject().value("source-layer").toString();
        // qDebug() << obj.toObject().value("id").toString();

        auto paint_object = obj.toObject().value("paint").toObject();
        auto filter_data = obj.toObject().value("filter").toArray();
        const bool is_line = obj.toObject().value("type").toString() == "line";
        const auto layer_name = std::make_pair(obj.toObject().value("source-layer").toString().toStdString(), is_line ? 0 : 1);

        std::vector<std::pair<uint8_t, uint32_t>> fill_colors;
        std::vector<std::pair<uint8_t, uint32_t>> outline_colors;
        std::vector<std::pair<uint8_t, uint16_t>> widths;
        std::vector<std::pair<uint8_t, std::pair<uint8_t, uint8_t>>> dashes;
        std::vector<std::pair<uint8_t, uint8_t>> opacities;

        // "stop" expressions may have a "base" value. this base value is used to manipulate the interpolation between values
        float fill_color_interpolation_base = 1;
        float outline_color_interpolation_base = 1;
        float width_interpolation_base = 1;
        float dash_interpolation_base = 1;
        float opacity_interpolation_base = 1;

        std::shared_ptr<StyleExpressionBase> filter = StyleExpressionBase::create_filter_expression(filter_data);

        const auto current_layer_filter = std::make_pair(layer_name, filter);

        if (previous_layer_filter.first.first.empty()) {
            // first time -> only set the previous_layer_filter
            previous_layer_filter = current_layer_filter;
        } else {
            // only increase layer_index if filter and or name is different
            // by doing this we prevent edge cases like https://github.com/AlpineMapsOrg/renderer/issues/151#issuecomment-2723695519 from overriding styles
            if (current_layer_filter != previous_layer_filter) {
                // add the previous values to the data structures
                layerid_filter_to_layer_indices[previous_layer_filter].push_back(zoom_to_style.size());
                zoom_to_style.push_back(std::move(current_style_map));
                // renew the current style map
                current_style_map = std::map<uint8_t, LayerStyle>();
                previous_layer_filter = current_layer_filter;
            }
        }

        bool invalid = false;

        for (const QString& key : paint_object.keys()) {

            // fill-antialias ?
            if (key == "fill-color") {
                parse_colors(paint_object.value(key), fill_colors, fill_color_interpolation_base);
            } else if (key == "line-color") {
                parse_colors(paint_object.value(key), fill_colors, fill_color_interpolation_base);
            } else if (key == "fill-opacity" || key == "line-opacity") {
                parse_opacities(paint_object.value(key), opacities, opacity_interpolation_base);
            } else if (key == "fill-outline-color") { // TODO visualize (only used for buildings)
                parse_colors(paint_object.value(key), outline_colors, outline_color_interpolation_base);
            } else if (key == "line-width") { //  "line-gap-width"
                parse_line_widths(paint_object.value(key), widths, width_interpolation_base);
            } else if (key == "line-dasharray") {
                parse_dashes(paint_object.value(key).toArray(), dashes, dash_interpolation_base);
            } else if (key == "line-offset") {
                // might be needed
            } else if (key == "fill-translate" || key == "line-translate-anchor") {
                // currentley not supported
            } else if (key == "fill-pattern") {
                // currently not supported -> causes errors when parsed
                invalid = true;
            } else if (key == "icon-color" || key == "circle-color" || key == "circle-radius" || key == "circle-stroke-color" || key == "circle-stroke-width" || key == "text-color"
                || key == "text-halo-color" || key == "text-halo-width" || key == "fill-antialias" || key == "line-gap-width") {
                // not used
            } else {
                qDebug() << "new unhandled style key detected: " << key.toStdString();
            }
        }

        if (invalid)
            continue;

        // determine zoom range defined in style.json (or fall back to constants::style_zoom_range (zoom in style.json only narrows the range)
        glm::uvec2 zoom_range = constants::style_zoom_range;
        if (obj.toObject().contains("minzoom") && uint8_t(obj.toObject().value("minzoom").toInt()) > zoom_range.x)
            zoom_range.x = obj.toObject().value("minzoom").toInt();
        if (obj.toObject().contains("maxzoom") && uint8_t(obj.toObject().value("maxzoom").toInt()) < zoom_range.y)
            zoom_range.y = obj.toObject().value("maxzoom").toInt();

        // determine style for every zoom level within range
        // fill arrays with default (max zoom, 0 value) if they are empty
        if (fill_colors.empty())
            fill_colors.push_back({ 255, 0 });
        if (outline_colors.empty())
            outline_colors.push_back({ 255, 0 });
        if (widths.empty())
            widths.push_back({ 255, 0 });
        if (dashes.empty())
            dashes.push_back({ 255, { 0, 0 } });
        if (opacities.empty())
            opacities.push_back({ 255, 255 });

        // use the first value of each array to determine the prev/current values
        std::pair<uint8_t, uint32_t> fill_colors_previous_value = fill_colors.front();
        std::pair<uint8_t, uint32_t> outline_colors_previous_value = outline_colors.front();
        std::pair<uint8_t, uint16_t> widths_previous_value = widths.front();
        std::pair<uint8_t, std::pair<uint8_t, uint8_t>> dashes_previous_value = dashes.front();
        std::pair<uint8_t, uint8_t> opacities_previous_value = opacities.front();

        std::pair<uint8_t, uint32_t> fill_colors_current_value = fill_colors.front();
        std::pair<uint8_t, uint32_t> outline_colors_current_value = outline_colors.front();
        std::pair<uint8_t, uint16_t> widths_current_value = widths.front();
        std::pair<uint8_t, std::pair<uint8_t, uint8_t>> dashes_current_value = dashes.front();
        std::pair<uint8_t, uint8_t> opacities_current_value = opacities.front();

        uint8_t fill_colors_index = 0;
        uint8_t outline_colors_index = 0;
        uint8_t widths_index = 0;
        uint8_t dashes_index = 0;
        uint8_t opacities_index = 0;

        bool previous_small_line = false;
        LayerStyle previous_style;

        for (unsigned zoom = zoom_range.x; zoom < zoom_range.y + 1; zoom++) {
            // determine if we have to change current / prev values
            while (fill_colors_current_value.first < zoom) {
                fill_colors_previous_value = fill_colors_current_value;
                fill_colors_current_value = fill_colors[++fill_colors_index];
            }
            while (outline_colors_current_value.first < zoom) {
                outline_colors_previous_value = outline_colors_current_value;
                outline_colors_current_value = outline_colors[++outline_colors_index];
            }
            while (widths_current_value.first < zoom) {
                widths_previous_value = widths_current_value;
                widths_current_value = widths[++widths_index];
            }
            while (dashes_current_value.first < zoom) {
                dashes_previous_value = dashes_current_value;
                dashes_current_value = dashes[++dashes_index];
            }
            while (opacities_current_value.first < zoom) {
                opacities_previous_value = opacities_current_value;
                opacities_current_value = opacities[++opacities_index];
            }

            // interpolate between previous and current values
            auto interpolation_factor_fill_color = 1.0;
            auto interpolation_factor_outline_color = 1.0;
            auto interpolation_factor_width = 1.0;
            auto interpolation_factor_dash = 1.0;
            auto interpolation_factor_opacity = 1.0;

            if (fill_colors_previous_value.second != fill_colors_current_value.second)
                interpolation_factor_fill_color = interpolation_factor(zoom, fill_color_interpolation_base, fill_colors_previous_value.first, fill_colors_current_value.first);
            if (outline_colors_previous_value.second != outline_colors_current_value.second)
                interpolation_factor_outline_color = interpolation_factor(zoom, outline_color_interpolation_base, outline_colors_previous_value.first, outline_colors_current_value.first);
            if (widths_previous_value.second != widths_current_value.second)
                interpolation_factor_width = interpolation_factor(zoom, width_interpolation_base, widths_previous_value.first, widths_current_value.first);
            if (dashes_previous_value.second != dashes_current_value.second)
                interpolation_factor_dash = interpolation_factor(zoom, dash_interpolation_base, dashes_previous_value.first, dashes_current_value.first);
            if (opacities_previous_value.second != opacities_current_value.second)
                interpolation_factor_opacity = interpolation_factor(zoom, opacity_interpolation_base, opacities_previous_value.first, opacities_current_value.first);

            uint32_t fill_color = interpolate_color(interpolation_factor_fill_color, fill_colors_previous_value.second, fill_colors_current_value.second);
            uint32_t outline_color = interpolate_color(interpolation_factor_outline_color, outline_colors_previous_value.second, outline_colors_current_value.second);
            uint16_t width = widths_previous_value.second * (1.0 - interpolation_factor_width) + widths_current_value.second * interpolation_factor_width;
            std::pair<uint8_t, uint8_t> dash = { dashes_previous_value.second.first * (1.0 - interpolation_factor_dash) + dashes_current_value.second.first * interpolation_factor_dash,
                dashes_previous_value.second.second * (1.0 - interpolation_factor_dash) + dashes_current_value.second.second * interpolation_factor_dash };
            uint8_t opacity = opacities_previous_value.second * (1.0 - interpolation_factor_opacity) + opacities_current_value.second * interpolation_factor_opacity;

            // we are not interested in very small lines below a certain zoom level
            constexpr float small_line_scale = constants::tile_extent / 256.0;
            bool small_line
                = (is_line && zoom < constants::small_line_zoom_threshold && width < constants::small_line_px * small_line_scale * constants::style_precision);
            if (small_line) {
                // we want to slowly blend the new line in
                opacity = 0u;
                width = 0u;
            }

            // merge opacity with colors
            if (opacity != 255u) {
                // if opacity is set it overrides any opacity from the color

                fill_color &= 4294967040u; // bit mask that zeros out the opacity bits
                fill_color |= opacity;

                outline_color &= 4294967040u; // bit mask that zeros out the opacity bits
                outline_color |= opacity;
            }

            // premultiply alpha
            // done outside of above if because opacity might be declared in fill_color only
            fill_color = premultiply_alpha(fill_color);
            outline_color = premultiply_alpha(outline_color);

            // dashes consist of gap / dash -> we combine them into one single value that is split again on the shader
            const uint16_t merged_dash = (dash.first << 8) | dash.second;

            if (small_line) {
                // this is a small line and we are not sure if we want to save this style for blending
                // we only want to store the style if the next style is visible
                previous_style = { fill_color, outline_color, width, merged_dash };
                previous_small_line = true;
            } else {
                if (previous_small_line) {
                    // this is no small line anymore -> we want to store the previous style for blending
                    current_style_map[zoom - 1] = previous_style;
                    previous_small_line = false;
                }
                // store the current style
                current_style_map[zoom] = { fill_color, outline_color, width, merged_dash };
            }

            //  DEBUG -> style_index to layername
            // auto id = obj.toObject().value("id").toString();
            // qDebug() << last_style_index << id;
        }
    }

    // at the end add the last filled style map
    layerid_filter_to_layer_indices[previous_layer_filter].push_back(uint32_t(zoom_to_style.size()));
    zoom_to_style.push_back(current_style_map);

    for (const auto& [key, layer_indices] : layerid_filter_to_layer_indices) {

        // we now know that we have valid styles -> create a new StyleFilter for this layer
        // if (!m_layer_to_style.contains(key.first))
        //     m_layer_to_style[key.first] = StyleFilter();

        for (const auto& layer_index : layer_indices) {

            const auto style_map = zoom_to_style[layer_index];
            LayerStyle current_style {};
            bool all_styles_same = true;
            // first loop only checks if all styles are the same -> if so, we do not need any blending and only one style is sufficient
            for (const auto& [zoom, style] : style_map) {
                if (current_style.fill_color == 0) // fill current_style in first iteration
                    current_style = LayerStyle(style);

                // check if any style looks different than first style
                if (current_style != style) {
                    all_styles_same = false;
                    break;
                }
            }

            // create styles below min zoom and fade out
            if (!all_styles_same) {
                // we might need to fill from styles from style.json range to style_zoom_range
                // we only need to add the styles, but we DO NOT need to add them to the m_layer_to_style
                // -> according to style.json there is no style for those values, we only need to add them for blending purposes

                uint8_t first_zoom = style_map.begin()->first;
                const auto first_style = style_map.at(first_zoom);

                for (uint8_t zoom = first_zoom - constants::mipmap_levels + 1; zoom < first_zoom; zoom++) {
                    style_values.push_back({ 0, first_style.outline_color, first_style.outline_width, first_style.outline_dash });
                }
            }

            bool first = true;
            for (const auto& [zoom, style] : style_map) {
                if (!all_styles_same || first) {

                    // add a new style every loop iteration if styles are different
                    // or if styles are the same only add it at first iteration
                    first = false;
                    style_values.push_back({ style.fill_color, style.outline_color, style.outline_width, style.outline_dash });
                }
                // add the styles to the data structure where we later can find the relevant style_index
                uint32_t style_index = (style_values.size() - 1u) << 1; // move style index by 1 for the "blend" flag
                m_layer_to_style[key.first].add_filter({ style_index | ((all_styles_same) ? 0u : 1u), layer_index, key.second }, zoom);
            }

            if (!all_styles_same) {
                // we might need to fill from styles from style.json range to style_zoom_range
                // we only need to add the styles, but we DO NOT need to add them to the m_layer_to_style
                // -> according to style.json there is no style for those values, we only need to add them for blending purposes

                const auto last_style = style_values[style_values.size() - 1];
                uint8_t last_zoom = style_map.rbegin()->first;

                for (uint8_t zoom = last_zoom; zoom < constants::style_zoom_range.y; zoom++) {
                    style_values.push_back({ last_style.x, last_style.y, last_style.z, last_style.w });
                }
            }

            // make sure that layer_index also fits into style_bits (we use this in preprocess)
            assert(layer_index < ((1u << constants::style_bits) - 1u));
        }
    }

    // qDebug() << "style_values: " << style_values.size();

    // make sure that the style values are within the buffer size; resize them to this size and create the raster images
    assert(style_values.size() <= constants::style_buffer_size * constants::style_buffer_size);
    style_values.resize(constants::style_buffer_size * constants::style_buffer_size, glm::u32vec4(-1u));

    m_styles = std::make_shared<const nucleus::Raster<glm::u32vec4>>(nucleus::Raster<glm::u32vec4>(constants::style_buffer_size, std::move(style_values)));

    qDebug() << "vectorlayer style loaded";
}

std::vector<std::pair<uint32_t, uint32_t>> Style::indices(std::string layer_name,
    int type,
    unsigned zoom,
    const mapbox::vector_tile::feature& feature,
    std::array<int, constants::max_style_expression_keys>* temp_values) const
{
    const auto layer = std::make_pair(layer_name, type);

    if (!m_layer_to_style.contains(layer)) {
        // qDebug() << "no style for: " << layer_name;
        return {};
    }

    return m_layer_to_style.at(layer).indices(zoom, feature, temp_values);
}

std::shared_ptr<const nucleus::Raster<glm::u32vec4>> Style::styles() const { return m_styles; }

// NOTE std::stof uses the locale to convert strings
// locale might be german and it expects a "," decimal
// we howewer want to force an english "." decimal point
// we could also temporarily change the locale and change it back afterwards
// but changing the locale might cause performance problems or other unexpected problems
// snippet from: https://stackoverflow.com/a/78993592 -> also talks about possible multi thread issues with changing locale
float Style::stringToFloat(const std::string& value)
{
    auto index = value.find(".");
    if (index == std::string::npos) {
        return std::stoi(value);
    }
    int full = std::stoi(value.substr(0, index));
    int decimals = std::stoi(value.substr(index + 1));
    return full + double(decimals / pow(10, value.substr(index + 1).size()));
}

float Style::rgb2linear(uint8_t channel)
{ // https://stackoverflow.com/a/21010385
    float s = channel / 255.0f;

    return s <= 0.04045 ? s / 12.92 : pow((s + 0.055) / 1.055, 2.4);
}
uint8_t Style::linear2rgb(float linear)
{ // https://stackoverflow.com/a/21010385
    float s = linear <= 0.0031308 ? linear * 12.92 : 1.055 * pow(linear, 1.0 / 2.4) - 0.055;
    return (uint8_t)(s * 255);
}

uint32_t Style::interpolate_color(float t, uint32_t color1, uint32_t color2)
{ // https://stackoverflow.com/a/21010385
    uint32_t interpolated = 0;
    for (int i = 0; i < 4; ++i) {
        const auto channel_bits = ((3 - i) * 8);
        float c1 = rgb2linear((color1 >> channel_bits) & 255);
        float c2 = rgb2linear((color2 >> channel_bits) & 255);

        interpolated |= linear2rgb((c2 - c1) * t + c1) << channel_bits;
    }

    return interpolated;
}

// uses https://github.com/maplibre/maplibre-style-spec/blob/main/src/expression/definitions/interpolate.ts -> exponentialInterpolation()
float Style::interpolation_factor(uint8_t zoom, float base, uint8_t zoom1, uint8_t zoom2)
{
    const auto diff = zoom2 - zoom1;
    uint8_t progress = 0;
    if (zoom > zoom1) // make sure that progress is a positive number or keep at 0
        progress = zoom - zoom1;

    if (diff == 0) // both values are the same
        return 0;

    if (base == 1)
        return std::clamp((double(progress) / double(diff)), 0.0, 1.0);
    else
        return std::clamp((pow(base, progress) - 1) / (pow(base, diff) - 1), 0.0, 1.0);
}

// template <typename T>
// T Style::interpolate(uint8_t zoom, float base, std::pair<uint8_t, T> prev, std::pair<uint8_t, T> current)
// {
//     const auto diff = current.first - prev.first;
//     const auto progress = zoom - diff;

//     if (diff == 0) // both values are the same
//         return current.second;

//     float t = 0;
//     if (base == 1)
//         t = float(progress) / float(diff);
//     else
//         t = (pow(base, progress) - 1) / (pow(base, diff) - 1);

//     return prev.second * (1.0 - t) + current.second * t;
// }

// // interpolate with pairs (for dashes)
// std::pair<uint8_t, uint8_t> Style::interpolate(uint8_t zoom, float base, std::pair<uint8_t, std::pair<uint8_t, uint8_t>> prev, std::pair<uint8_t, std::pair<uint8_t, uint8_t>> current)
// {
//     const auto diff = current.first - prev.first;
//     const auto progress = zoom - diff;

//     if (diff == 0) // both values are the same
//         return current.second;

//     float t = 0;
//     if (base == 1)
//         t = float(progress) / float(diff);
//     else
//         t = (pow(base, progress) - 1) / (pow(base, diff) - 1);

//     return { prev.second.first * (1.0 - t) + current.second.first * t, prev.second.second * (1.0 - t) + current.second.second * t };
// }

void Style::parse_colors(const QJsonValue& value, std::vector<std::pair<uint8_t, uint32_t>>& colors, float& base)
{
    assert(colors.size() == 0); // make sure that we only parse attribute once
    if (value.isObject() && value.toObject().contains("stops")) {
        if (value.toObject().contains("base"))
            base = value.toObject().value("base").toDouble(1);

        const auto stops = value.toObject().value("stops").toArray();
        for (qsizetype i = 0; i < stops.size(); i++) {
            const auto stop = stops[i].toArray();
            colors.push_back({ stop[0].toInt(), parse_color(stop[1]) });
        }

        // last value repeats the last value and sets zoom to max
        colors.push_back({ 255, colors.back().second });
    } else {
        // only one value -> add element with highest zoom value
        colors.push_back({ 255, parse_color(value) });
    }
}
void Style::parse_opacities(const QJsonValue& value, std::vector<std::pair<uint8_t, uint8_t>>& opacities, float& base)
{
    assert(opacities.size() == 0); // make sure that we only parse attribute once
    if (value.isObject() && value.toObject().contains("stops")) {
        if (value.toObject().contains("base"))
            base = value.toObject().value("base").toDouble(1);

        const auto stops = value.toObject().value("stops").toArray();
        for (qsizetype i = 0; i < stops.size(); i++) {
            const auto stop = stops[i].toArray();
            assert(stop[0].isDouble()); // value is not a number -> take a closer look
            opacities.push_back({ uint8_t(stop[0].toDouble()), parse_opacity(stop[1]) });
        }

        // last value repeats the last value and sets zoom to max
        opacities.push_back({ 255, opacities.back().second });
    } else {
        // only one value -> add element with highest zoom value
        opacities.push_back({ 255, parse_opacity(value) });
    }
}
void Style::parse_dashes(const QJsonValue& value, std::vector<std::pair<uint8_t, std::pair<uint8_t, uint8_t>>>& dashes, float& base)
{
    assert(dashes.size() == 0); // make sure that we only parse attribute once
    if (value.isObject() && value.toObject().contains("stops")) {
        if (value.toObject().contains("base"))
            base = value.toObject().value("base").toDouble(1);

        const auto stops = value.toObject().value("stops").toArray();
        for (qsizetype i = 0; i < stops.size(); i++) {
            const auto stop = stops[i].toArray();
            dashes.push_back({ stop[0].toInt(), parse_dash(stop[1]) });
        }

        // last value repeats the last value and sets zoom to max
        dashes.push_back({ 255, dashes.back().second });
    } else {
        // only one value -> add element with highest zoom value
        dashes.push_back({ 255, parse_dash(value) });
    }
}
void Style::parse_line_widths(const QJsonValue& value, std::vector<std::pair<uint8_t, uint16_t>>& widths, float& base)
{
    assert(widths.size() == 0); // make sure that we only parse attribute once
    if (value.isObject() && value.toObject().contains("stops")) {
        if (value.toObject().contains("base"))
            base = value.toObject().value("base").toDouble(1);

        const auto stops = value.toObject().value("stops").toArray();
        for (qsizetype i = 0; i < stops.size(); i++) {
            const auto stop = stops[i].toArray();
            widths.push_back({ stop[0].toInt(), parse_line_width(stop[1]) });
        }

        // last value repeats the last value and sets zoom to max
        widths.push_back({ 255, widths.back().second });
    } else {
        // only one value -> add element with highest zoom value
        widths.push_back({ 255, parse_line_width(value) });
    }
}

uint32_t Style::parse_color(const QJsonValue& value)
{
    std::string colorValue;
    if (value.isString()) {
        colorValue = value.toString().toStdString();
    } else {
        qDebug() << "cannot parse color value: " << value;
        assert(false);
        return 0ul;
    }

    if (colorValue.starts_with("#")) {
        if (colorValue.length() == 4) // transform #9CF to #99CCFF
            colorValue = "#" + std::string(2, colorValue[1]) + std::string(2, colorValue[2]) + std::string(2, colorValue[3]);

        if (colorValue.length() == 7)
            return (std::stoul(colorValue.substr(1), nullptr, 16) << 8) | 255;
        else if (colorValue.length() == 9)
            return std::stoul(colorValue.substr(1), nullptr, 16);
        else {
            qDebug() << "cannot parse hex color: " << colorValue;
            return 0ul;
        }
    } else if (colorValue.starts_with("rgb")) {
        // parses rgb(int,int,int) and rgba(int,int,int,float)
        // ints in range [0-255]; float in range [0-1]
        uint32_t out = 0u;

        auto startPos = colorValue.find("(");
        auto tmp = colorValue.substr(startPos + 1, colorValue.size() - startPos - 2);
        int count = 0;
        auto pos = tmp.find(',');
        while (pos != std::string::npos) {
            count++;
            out = out << 8;

            // find start of next digit or use the full colorValue for the rest
            auto nextPos = tmp.find(',');
            auto nextVal = tmp;
            if (nextPos != std::string::npos)
                nextVal = tmp.substr(0, nextPos);

            // count < 4 necessary since the alpha colorValue might be 1 -> and has to be multiplied by 255
            if (nextVal.find(".") == std::string::npos && count < 4) {
                // integer
                out |= std::stoul(nextVal);
            } else {
                // decimal
                out |= uint32_t(stringToFloat(nextVal) * 255.f);
            }
            // remove the digit we just parsed
            tmp = tmp.substr(nextPos + 1);

            pos = nextPos;
        }

        if (count == 3) // only rgb was given -> add full transparancy
            out = (out << 8) | 255;

        return out;
    } else if (colorValue.starts_with("hsl")) {
        const std::regex regex("hsla?\\((\\d+),\\s?(\\d+)%,\\s?(\\d+)%(?:,\\s?(\\d+.?\\d*))?\\)");
        std::smatch matches;
        if (!std::regex_match(colorValue, matches, regex)) {
            qDebug() << "could not match hsl regex" << colorValue;
            assert(false);
            return 0ul;
        }

        const float h = std::stoi(matches[1]) % 360;
        const float s = std::stoi(matches[2]) / 100.0;
        const float l = std::stoi(matches[3]) / 100.0;

        uint8_t a = 255;
        if (matches.size() == 5 && matches[4].length() > 0) {
            a = stringToFloat(matches[4]) * 255.0;
        }

        // formula from https://www.rapidtables.com/convert/color/hsl-to-rgb.html
        const float c = (1.0 - abs(2.0 * l - 1.0)) * s;
        const float x = c * (1.0 - abs(fmod(h / 60.0, 2) - 1.0));
        const float m = l - c / 2.0;

        glm::uvec3 tmp_rgb;
        if (h < 60)
            tmp_rgb = glm::uvec3(round((c + m) * 255.0), round((x + m) * 255.0), round((0.0 + m) * 255.0));
        else if (h < 120)
            tmp_rgb = glm::uvec3(round((x + m) * 255.0), round((c + m) * 255.0), round((0.0 + m) * 255.0));
        else if (h < 180)
            tmp_rgb = glm::uvec3(round((0.0 + m) * 255.0), round((c + m) * 255.0), round((x + m) * 255.0));
        else if (h < 240)
            tmp_rgb = glm::uvec3(round((0.0 + m) * 255.0), round((x + m) * 255.0), round((c + m) * 255.0));
        else if (h < 300)
            tmp_rgb = glm::uvec3(round((x + m) * 255.0), round((0.0 + m) * 255.0), round((c + m) * 255.0));
        else
            tmp_rgb = glm::uvec3(round((c + m) * 255.0), round((0.0 + m) * 255.0), round((x + m) * 255.0));

        return tmp_rgb.x << 24 | tmp_rgb.y << 16 | tmp_rgb.z << 8 | a;
    } else {
        qDebug() << "cannot parse color: " << colorValue;
        return 0ul;
    }
}

uint32_t Style::premultiply_alpha(uint32_t color)
{
    uint8_t a = color & 255;
    float opacity = float(a) / 255.0;
    uint8_t r = ((color >> 24) & 255) * opacity;
    uint8_t g = ((color >> 16) & 255) * opacity;
    uint8_t b = ((color >> 8) & 255) * opacity;

    return ((r << 24) | (g << 16) | (b << 8) | a);
}

uint8_t Style::parse_opacity(const QJsonValue& value)
{

    if (value.isDouble() && value.toDouble() <= 1.0) {
        return value.toDouble() * 255;
    }

    qDebug() << "unhandled opacity value" << value;
    assert(false);
    return 255;
}

/*
 * Note so far no useful documentation was found for dash-array and how they are structured/constructed
 * -> therefore verify that we are using this correctly
 * Currently: we only support 2 values for dash array -> dash and gap
 * we assume that index 0 is the dash size and index 1 is the gap size.
 * https://docs.mapbox.com/android/maps/api/10.2.0/mapbox-maps-android/com.mapbox.maps.plugin.annotation.generated/-polyline-annotation-manager/line-dasharray.html
 * declares that those values are multiplied by line width to the actual size
 * we store the ratio between both values and the sum of both values in one uint32_t value
 * this currently wastes a bit of space -> two 8bit values should suffice here
 * but we are also not quite clear about all the possible values
 * -> THEREFORE TODO veryfy the assumptions of this method
 */
std::pair<uint8_t, uint8_t> Style::parse_dash(const QJsonValue&)
{

    // TODO there are also values with more than two values
    // currently only 2 values are allowed here
    // assert(dash_values.size() == 2);

    // TODO implement

    return { 0, 0 };

    // double sum = dash_values[0].toDouble() + dash_values[1].toDouble();

    // uint16_t dash_gap_ratio = dash_values[0].toDouble() / sum * 65535.f;

    // return dash_gap_ratio << 16 | uint16_t(sum);
}

uint16_t Style::parse_line_width(const QJsonValue& value)
{
    if (value.isDouble()) {
        float thickness = value.toDouble() * constants::line_width_multiplier;
        if (thickness > constants::max_line_width)
            thickness = constants::max_line_width;
        return uint16_t(thickness * constants::style_precision);
    }

    qDebug() << "unhandled line width value" << value;
    assert(false);
    return 0;
}

} // namespace nucleus::vector_layer
