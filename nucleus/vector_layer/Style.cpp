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

namespace {
static inline uint8_t clamp_u8(int v) { return static_cast<uint8_t>(std::min(std::max(v, 0), 255)); }

static inline float srgb_to_linear(float cs)
{
    // cs in [0,1]
    return (cs <= 0.04045f) ? (cs / 12.92f) : std::pow((cs + 0.055f) / 1.055f, 2.4f);
}

} // namespace

namespace nucleus::vector_layer {

// https://maplibre.org/maplibre-style-spec/

Style::Style(const QString& filename)
    : m_filename(filename)
{
    // make sure that style_bits and style_buffer_size are correctly matched -> changing one means you also need to change the other
    assert((((constants::style_buffer_size * constants::style_buffer_size) / (constants::style_zoom_range.y + 1)) >> (constants::style_bits)) <= 1);
    StyleExpression::initialize();
}

Style::Style(Style&& other)
    : m_styles(std::move(other.m_styles))
    , m_visible_styles(std::move(other.m_visible_styles))
    , m_layer_to_style(std::move(other.m_layer_to_style))
    , m_lowest_encountered_zoom(std::move(other.m_lowest_encountered_zoom))
    , m_styles_to_update(std::move(other.m_styles_to_update))
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


    QJsonDocument doc = QJsonDocument::fromJson(data);
    QJsonArray layers = style_expander::expand(doc.object().value("layers").toArray());

    // key pair is <layer_id, filter>
    // layerid: source_layer+_+type (e.g. transportation_line)
    // example: transportation_line + filter(bridge, class:secondary) -> vector<(zoom_to_style indices)>
    auto layerid_filter_to_layer_indices
        = std::map<std::pair<std::pair<std::string, int>, std::shared_ptr<StyleExpressionBase>>, std::vector<uint32_t>, StyleHasher>();
    // each entry in this vector is a different layer_index
    // vector<map<zoom, style>>
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
        auto layout_object = obj.toObject().value("layout").toObject();
        auto filter_data = obj.toObject().value("filter").toArray();
        const bool is_line = obj.toObject().value("type").toString() == "line";
        const auto layer_name = std::make_pair(obj.toObject().value("source-layer").toString().toStdString(), is_line ? 0 : 1);

        std::vector<std::pair<uint8_t, uint32_t>> fill_colors;
        std::vector<std::pair<uint8_t, uint16_t>> widths;
        std::vector<std::pair<uint8_t, std::pair<uint8_t, float>>> dashes;
        std::vector<std::pair<uint8_t, uint8_t>> opacities;

        // "stop" expressions may have a "base" value. this base value is used to manipulate the interpolation between values
        float fill_color_interpolation_base = 1;
        float width_interpolation_base = 1;
        float dash_interpolation_base = 1;
        float opacity_interpolation_base = 1;

        std::shared_ptr<StyleExpressionBase> filter = StyleExpressionBase::create_filter_expression(filter_data);

        const auto current_layer_filter = std::make_pair(layer_name, filter);

        // TODO is this still necessary? I think this is deprecated and was only used previously for merging similar styles/filters
        if (previous_layer_filter.first.first.empty()) {
            // first time -> only set the previous_layer_filter
            previous_layer_filter = current_layer_filter;
        } else {
            // only increase layer_index if filter and or name is different
            // by doing this we prevent edge cases like https://github.com/AlpineMapsOrg/renderer/issues/151#issuecomment-2723695519 from overriding styles
            if (current_layer_filter != previous_layer_filter) {
                if (!current_style_map.empty()) {
                    // add the previous values to the data structures
                    layerid_filter_to_layer_indices[previous_layer_filter].push_back(zoom_to_style.size());
                    zoom_to_style.push_back(std::move(current_style_map));
                }
                // renew the current style map
                current_style_map = std::map<uint8_t, LayerStyle>();
                previous_layer_filter = current_layer_filter;
            }
        }

        bool invalid = false;

        // default line-cap is butt (https://maplibre.org/maplibre-style-spec/layers/#line-cap)
        // additionally we also currently do not support square line-caps (only one style uses this, and this style probably isn't worth the shader rewrite)
        bool round_line_cap = layout_object.contains("line-cap") && layout_object.value("line-cap") == "round";

        for (const QString& key : paint_object.keys()) {

            // fill-antialias ?
            if (key == "fill-color") {
                parse_colors(paint_object.value(key), fill_colors, fill_color_interpolation_base);
            } else if (key == "line-color") {
                parse_colors(paint_object.value(key), fill_colors, fill_color_interpolation_base);
            } else if (key == "fill-opacity" || key == "line-opacity") {
                parse_opacities(paint_object.value(key), opacities, opacity_interpolation_base);
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
            } else if (key == "fill-outline-color" || key == "icon-color" || key == "circle-color" || key == "circle-radius" || key == "circle-stroke-color"
                || key == "circle-stroke-width" || key == "text-color" || key == "text-halo-color" || key == "text-halo-width" || key == "fill-antialias"
                || key == "line-gap-width") {
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
        if (widths.empty())
            widths.push_back({ 255, 0 });
        if (dashes.empty())
            dashes.push_back({ 255, { 1 * constants::style_precision, 1 } });
        if (opacities.empty())
            opacities.push_back({ 255, 255 });

        // use the first value of each array to determine the prev/current values
        std::pair<uint8_t, uint32_t> fill_colors_previous_value = fill_colors.front();
        std::pair<uint8_t, uint16_t> widths_previous_value = widths.front();
        std::pair<uint8_t, std::pair<uint8_t, float>> dashes_previous_value = dashes.front();
        std::pair<uint8_t, uint8_t> opacities_previous_value = opacities.front();

        std::pair<uint8_t, uint32_t> fill_colors_current_value = fill_colors.front();
        std::pair<uint8_t, uint16_t> widths_current_value = widths.front();
        std::pair<uint8_t, std::pair<uint8_t, float>> dashes_current_value = dashes.front();
        std::pair<uint8_t, uint8_t> opacities_current_value = opacities.front();

        uint8_t fill_colors_index = 0;
        uint8_t widths_index = 0;
        uint8_t dashes_index = 0;
        uint8_t opacities_index = 0;

        for (unsigned zoom = zoom_range.x; zoom < zoom_range.y + 1; zoom++) {
            // determine if we have to change current / prev values
            while (fill_colors_current_value.first < zoom) {
                fill_colors_previous_value = fill_colors_current_value;
                fill_colors_current_value = fill_colors[++fill_colors_index];
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
            auto interpolation_factor_width = 1.0;
            auto interpolation_factor_dash = 1.0;
            auto interpolation_factor_opacity = 1.0;

            if (fill_colors_previous_value.second != fill_colors_current_value.second)
                interpolation_factor_fill_color = interpolation_factor(zoom, fill_color_interpolation_base, fill_colors_previous_value.first, fill_colors_current_value.first);
            if (widths_previous_value.second != widths_current_value.second)
                interpolation_factor_width = interpolation_factor(zoom, width_interpolation_base, widths_previous_value.first, widths_current_value.first);
            if (dashes_previous_value.second != dashes_current_value.second)
                interpolation_factor_dash = interpolation_factor(zoom, dash_interpolation_base, dashes_previous_value.first, dashes_current_value.first);
            if (opacities_previous_value.second != opacities_current_value.second)
                interpolation_factor_opacity = interpolation_factor(zoom, opacity_interpolation_base, opacities_previous_value.first, opacities_current_value.first);

            uint32_t fill_color = interpolate_color(interpolation_factor_fill_color, fill_colors_previous_value.second, fill_colors_current_value.second);
            uint16_t width = widths_previous_value.second * (1.0 - interpolation_factor_width) + widths_current_value.second * interpolation_factor_width;
            std::pair<uint8_t, float> dash = { dashes_previous_value.second.first * (1.0 - interpolation_factor_dash)
                    + dashes_current_value.second.first * interpolation_factor_dash,
                dashes_previous_value.second.second * (1.0 - interpolation_factor_dash) + dashes_current_value.second.second * interpolation_factor_dash };
            uint8_t opacity = opacities_previous_value.second * (1.0 - interpolation_factor_opacity) + opacities_current_value.second * interpolation_factor_opacity;

            // merge opacity with colors
            if (opacity != 255u) {
                // if opacity is set it overrides any opacity from the color

                fill_color &= remove_alpha_mask; // bit mask that zeros out the opacity bits
                fill_color |= opacity;
            }

            // premultiply alpha + gamma_decode
            // done outside of above if because opacity might be declared in fill_color only
            fill_color = premultiply_alpha(gamma_decode(fill_color));

            // make sure that dash_sum with style_precision is not bigger than available bits
            uint8_t dash_sum = 255;
            const uint temp_dash_sum = dash.second * (width / constants::style_precision);
            if (temp_dash_sum <= 255.0) {
                dash_sum = std::max(1u, temp_dash_sum); // guarantee that dash_sum is not 0
            } else {
                qDebug() << "Style: dash sum needs more bits (" << temp_dash_sum << " should be < 255)";
                assert(false);
            }

            // store the current style
            current_style_map[zoom] = { fill_color, width, dash.first, dash_sum, round_line_cap };

            //  DEBUG -> style_index to layername
            // auto id = obj.toObject().value("id").toString();
            // qDebug() << last_style_index << id;
        }
    }

    // at the end add the last filled style map
    if (current_style_map.size() > 0) {
        layerid_filter_to_layer_indices[previous_layer_filter].push_back(uint32_t(zoom_to_style.size()));
        zoom_to_style.push_back(current_style_map);
    }

    std::vector<std::vector<glm::u32vec2>> temp_styles;
    temp_styles.reserve(255);

    for (const auto& [key, layer_indices] : layerid_filter_to_layer_indices) {

        for (const auto& layer_index : layer_indices) {

            const auto style_map = zoom_to_style[layer_index];

            // create styles below min zoom and fade out
            // we might need to fill from styles from style.json range to style_zoom_range start
            // we only need to add the styles, but we DO NOT need to add them to the m_layer_to_style
            // -> according to style.json there is no style for those values, we only need to add them for blending purposes
            const uint first_zoom = style_map.begin()->first;
            const auto first_style = style_map.at(first_zoom);

            // const uint32_t style_index = style_values.size() / num_zooms_per_style;
            const uint32_t style_index = temp_styles.size();
            auto& current_style = temp_styles.emplace_back();
            current_style.resize(constants::num_zooms_per_style, glm::u32vec2(0u));

            for (size_t i = 0; i < first_zoom; i++) {
                // duplicate first style with alpha 0
                // since we premultiply alpha -> alpha: 0 = color: 0,0,0
                current_style[i] = { 0u, first_style.buffer_alignment().y };
            }

            for (const auto& [zoom, style] : style_map) {

                // add a new style every loop iteration
                current_style[zoom] = { style.buffer_alignment() };

                // add the styles to the data structure where we later can find the relevant style_index
                m_layer_to_style[key.first].add_filter({ { style_index, layer_index }, key.second }, zoom);
            }

            // we might need to fill from styles from style.json range to style_zoom_range end
            // we only need to add the styles, but we DO NOT need to add them to the m_layer_to_style
            // -> according to style.json there is no style for those values, we only need to add them for blending purposes
            const uint last_zoom = style_map.rbegin()->first;
            const auto last_style = current_style[last_zoom];

            for (uint i = last_zoom + 1; i < constants::num_zooms_per_style; i++) {
                // duplicate last style with alpha 0 (if we are not yet at maxzoom)
                // since we premultiply alpha -> alpha: 0 = color: 0,0,0
                current_style[i] = { 0, last_style.y };
                if (i == last_zoom + 1)
                    m_layer_to_style[key.first].add_filter({ { style_index, layer_index }, key.second }, i);
            }
            // set it to the highest value
            m_lowest_encountered_zoom[layer_index] = constants::num_zooms_per_style;

            // qDebug() << "style_index" << style_index;
        }
    }

    // put the styles to the correct position in the stylebuffer (that will be a 2d texture)
    auto style_values = create_style_buffer_data(temp_styles);

    // visible styles set the color to 0 for now -> only set to higher value if encountered from server
    // this allows to fade out styles that were not encountered on lower zoom levels
    auto visible_style_values = std::vector<glm::u32vec2>(style_values);
    for (auto& v : visible_style_values) {
        // set color alpha to 0 for fadeout
        // since we premultiply alpha -> alpha: 0 = color: 0,0,0
        v.x = 0;
    }

    // make sure that the style values are within the buffer size; resize them to this size and create the raster images
    assert(style_values.size() <= constants::style_buffer_size * constants::style_buffer_size);
    visible_style_values.resize(constants::style_buffer_size * constants::style_buffer_size, glm::u32vec2(-1u));

#ifdef ALP_ENABLE_DEBUG_VECTOR_TILES
    // add debug styles

    std::vector<LayerStyle> debug_styles {
        { 0x000000FF, 0, 1 * constants::style_precision, 1, false }, // black poly
        { 0xFFFFFFFF, 0, 1 * constants::style_precision, 1, false }, // white poly
        { 0x000000FF, uint(0.1 * constants::style_precision), 1 * constants::style_precision, 1, false }, // black line width 0.1
        { 0x000000FF, uint(0.5 * constants::style_precision), 1 * constants::style_precision, 1, false }, // black line width 0.5
        { 0x000000FF, 1 * constants::style_precision, 1 * constants::style_precision, 1, false }, // black line width 1
        { 0x000000FF, 5 * constants::style_precision, 1 * constants::style_precision, 1, false }, // black line width 5
        { 0x000000FF, 10 * constants::style_precision, 1 * constants::style_precision, 1, false }, // black line width 10
        { 0x000000FF, 15 * constants::style_precision, 1 * constants::style_precision, 1, false }, // black line width 15
        { 0x000000FF, 25 * constants::style_precision, 1 * constants::style_precision, 1, false }, // black line width 25
        { 0x000000FF, 50 * constants::style_precision, 1 * constants::style_precision, 1, false }, // black line width 50
    };

    const auto start_index = uint(style_values.size() / style_zoom_multiplier) - (debug_styles.size());

    // make sure that we do not override valid styles
    assert(style_values[start_index * style_zoom_multiplier] == glm::u32vec2(-1u));
    assert(visible_style_values[start_index * style_zoom_multiplier] == glm::u32vec2(-1u));

    for (uint i = 0; i < debug_styles.size(); i++) {
        const auto style = debug_styles[i].buffer_alignment();
        for (uint j = 0; j < constants::style_zoom_range.y + 1; j++) {

            style_values[(start_index + i) * style_zoom_multiplier + j] = style;
            visible_style_values[(start_index + i) * style_zoom_multiplier + j] = style;
        }
    }

    qDebug() << "debug style start: " << start_index;

#endif

    m_visible_styles
        = std::make_shared<nucleus::Raster<glm::u32vec2>>(nucleus::Raster<glm::u32vec2>(constants::style_buffer_size, std::move(visible_style_values)));
    m_styles = std::make_shared<const nucleus::Raster<glm::u32vec2>>(nucleus::Raster<glm::u32vec2>(constants::style_buffer_size, std::move(style_values)));

    // qDebug() << "vectorlayer style loaded";
}

std::vector<StyleLayerIndex> Style::indices(std::string layer_name,
    int type,
    unsigned zoom,
    const mapbox::vector_tile::feature& feature,
    std::array<int, constants::max_style_expression_keys>* temp_values)
{
    const auto layer = std::make_pair(layer_name, type);

    if (!m_layer_to_style.contains(layer)) {
        // qDebug() << "no style for: " << layer_name;
        return {};
    }

    const auto indices = m_layer_to_style.at(layer).indices(zoom, feature, temp_values);

    return indices;
}

std::vector<StyleLayerIndex> Style::simplify_styles(
    std::vector<StyleLayerIndex>* style_and_layer_indices, const uint zoom_level, const std::vector<glm::u32vec2>& style_buffer)
{
    // we get multiple styles that may have full opacity and the same width
    // creating render calls for both does not make sense -> we only want to draw the top layer
    // this function simplifys all the styles so that only the styles which actually have a change to be rendered will remain

    // TODO this sort should happen at creation of the vector not here
    // order the styles so that we look at layer in descending order
    std::sort(
        style_and_layer_indices->begin(), style_and_layer_indices->end(), [](StyleLayerIndex a, StyleLayerIndex b) { return a.layer_index > b.layer_index; });
    std::vector<StyleLayerIndex> out_styles;
    int accummulative_opacity = 0;
    float width = 0.0;

    for (const auto& indices : *style_and_layer_indices) {

        const auto buffer_index = Style::style_buffer_index(indices.style_index, std::min(zoom_level, zoom_level - 1u));
        const auto style_data_lower = style_buffer[buffer_index];
        const auto style_data_higher = style_buffer[buffer_index + 1];

        const float lower_width = Style::style_width(style_data_lower);
        const float lower_opacity = style_data_lower.x & 255;
        const float higher_width = Style::style_width(style_data_higher);
        const bool uses_dashes = Style::uses_dashes(style_data_higher);
        const float higher_opacity = style_data_higher.x & 255;

        // by mixing the lower and higher style -> we get a value that better represents a real world example
        // this is neccessary for e.g. 1 landcover style that stops at z12 and another that starts at z13
        // -> we need to render both, because rendering only one at z13 would falsely represent a fade to alpha 0 between z12 and z13
        // by using z12.5 for the current opacity and width, we can make sure that any can be countered by the fading in the opposite direction
        const float current_width = (lower_width + higher_width) / 2.0;
        const int current_opacity = (lower_opacity + higher_opacity) / 2.0;

        if (current_opacity == 0)
            continue; // we dont care about 0 opacity geometry

        if (width < current_width) {
            // reset opacity
            accummulative_opacity = 0;
            width = current_width;
        }

        if (accummulative_opacity < 255) {
            if (!uses_dashes) // dashes do not count for accummulative opacity
                accummulative_opacity += current_opacity;
            out_styles.push_back(indices);
        }
    }

    return out_styles;
}

std::vector<glm::u32vec2> Style::create_style_buffer_data(const std::vector<std::vector<glm::u32vec2>>& styles)
{
    // with this encoding we ensure the following:
    // - shader needs a lower and a higher style (but only those two indices for a fragment -> it doesnt care about the other zooms)
    //    -> we therefore duplicate the style and write it in the following column. but we ignore the zoom 0
    //    -> additionally we repeat the last zoom -> this minimizes special cases on the shader while only slightly increasing buffer needs
    // - it is possible that the next needed style will be the next index (and most definitely will have the same zoom if we look at the same fragment)
    //    -> we therefore want to place it on the same row, as near as possible to the previous style
    // - styles are looked at in descending order -> they need to be placed in the same descending order in the buffer
    //
    // example:
    // style_n lower = column 1 [row 0 - max]
    // style_n higher = column 2 [row -1 - (max-1)]
    // style_(n-1) lower = column 3 [row 0 - max]
    // style_(n-1) higher = column 4 [row -1 - (max-1)]

    std::vector<glm::u32vec2> out;
    out.resize(constants::style_buffer_size * constants::style_buffer_size, glm::u32vec2(-1u));

    // we currently only support up to 384 styles
    // -> 19 zooms per style -> styles duplicated for easier higher zoom fetching = max 6 rows * 64 columns
    // possible future improvements:
    // - we could theoretically have 7 rows if the last row is split (e.g. zoom 0-9 and 10-18 on 2 different column pairs)
    //    -> problem that fetching the index is a bit more complicated -> probably worse glsl performance
    // - ignore the first 3 zooms -> we only have 16 zooms and can support up to 8 rows (without wastage) = 512 styles
    assert(styles.size() < 385);

    int column = 0;
    int start_row = 0;
    for (size_t i = 0; i < styles.size(); i++) { // TODO !!!! remove again
        // for (size_t i = styles.size(); i-- > 0;) {
        for (size_t j = 0; j < constants::num_zooms_per_style; j++) {

            out[(start_row + j) * constants::style_buffer_size + column] = styles[i][j];
            if (j < constants::num_zooms_per_style - 1)
                out[(start_row + j) * constants::style_buffer_size + column + 1] = styles[i][j + 1];
            else
                // duplicate the last style
                out[(start_row + j) * constants::style_buffer_size + column + 1] = styles[i][j];
        }

        column += 2;
        if (column >= constants::style_buffer_size) {
            column = 0;
            start_row += constants::num_zooms_per_style;
        }
    }

    return out;
}

// we register the styles after we simplified them
// this way, the dynamic blending of visible styles is more true to what really is needed
void Style::register_used_styles(const uint zoom_level, const std::vector<StyleLayerIndex>& indices)
{
    for (const auto& index : indices) {

        if (zoom_level < m_lowest_encountered_zoom[index.layer_index]) {
            m_lowest_encountered_zoom[index.layer_index] = zoom_level;
            m_styles_to_update.push_back(index);
        }
    }
}

std::shared_ptr<const nucleus::Raster<glm::u32vec2>> Style::styles() const { return m_styles; }
std::shared_ptr<const nucleus::Raster<glm::u32vec2>> Style::visible_styles() const { return m_visible_styles; }

bool Style::update_visible_styles()
{
    if (m_styles_to_update.size() == 0)
        return false; // no update needed

    const auto& style_buffer = m_styles->buffer();
    auto& visible_style_buffer = m_visible_styles->buffer();
    for (const auto& indices : m_styles_to_update) {

        const auto encountered_zoom = m_lowest_encountered_zoom[indices.layer_index];

        // go one step below the lowest encountered zoom -> we need to update the higher zoom here
        const auto start_index = style_buffer_index(indices.style_index, std::min(encountered_zoom, encountered_zoom - 1));
        const auto updateable_styles = constants::num_zooms_per_style - encountered_zoom;

        // update the higher zoom only
        visible_style_buffer[start_index + 1] = style_buffer[start_index + 1];

        // for the rest, update lower and higher zoom
        for (size_t i = 1; i < updateable_styles; i++) {
            const auto index = start_index + (i * constants::style_buffer_size); // -> next zoom = next row of the buffer
            visible_style_buffer[index] = style_buffer[index];
            visible_style_buffer[index + 1] = style_buffer[index + 1];
        }
    }
    m_styles_to_update.clear();

    return true;
}

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

uint32_t Style::style_buffer_index(const uint32_t style_index, const uint zoom_level)
{
    // NOTE: (style_index << 1) necessary since we want to only address every second row (essentially we multiply the index by 2)
    const auto style_buffer_col = (style_index << 1) & (constants::style_buffer_size - 1u);
    const auto style_buffer_row = ((style_index >> (constants::bits_per_buffer_row - 1u)) * constants::num_zooms_per_style) + zoom_level;

    return style_buffer_col + (style_buffer_row * constants::style_buffer_size);
}

float Style::style_width(const glm::u32vec2& style)
{

    // if (style.y == -1u)
    //     qDebug() << "line_width is -1";
    return float(style.y >> 17) / float(constants::style_precision);
}

std::pair<float, float> Style::style_dashes(const glm::u32vec2& style)
{
    const auto dash_gap_ratio = (style.y >> 9) & ((1u << (17 - 9)) - 1u);
    const auto dash_sum = (style.y >> 1) & ((1u << (9 - 1)) - 1u);

    return std::make_pair(float(dash_gap_ratio) / float(constants::style_precision), float(dash_sum));
}

bool Style::uses_dashes(const glm::u32vec2& style) { return float((style.y & ((1u << 17) - 1u)) >> 9) < constants::style_precision; }

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
void Style::parse_dashes(const QJsonValue& value, std::vector<std::pair<uint8_t, std::pair<uint8_t, float>>>& dashes, float& base)
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

uint32_t Style::gamma_decode(uint32_t colour)
{
    // Unpack (0xRRGGBBAA)
    uint8_t r = static_cast<uint8_t>((colour >> 24) & 0xFF);
    uint8_t g = static_cast<uint8_t>((colour >> 16) & 0xFF);
    uint8_t b = static_cast<uint8_t>((colour >> 8) & 0xFF);
    uint8_t a = static_cast<uint8_t>(colour & 0xFF); // unchanged

    // Normalize to [0,1]
    double rf = r / 255.0;
    double gf = g / 255.0;
    double bf = b / 255.0;

    // sRGB gamma decode -> linear
    rf = srgb_to_linear(rf);
    gf = srgb_to_linear(gf);
    bf = srgb_to_linear(bf);

    // Convert back to 8-bit (still representing linear light)
    // Round to nearest and clamp
    uint8_t r_lin = clamp_u8(static_cast<int>(std::lround(rf * 255.0)));
    uint8_t g_lin = clamp_u8(static_cast<int>(std::lround(gf * 255.0)));
    uint8_t b_lin = clamp_u8(static_cast<int>(std::lround(bf * 255.0)));

    // Repack as 0xRRGGBBAA
    return (static_cast<uint32_t>(r_lin) << 24) | (static_cast<uint32_t>(g_lin) << 16) | (static_cast<uint32_t>(b_lin) << 8) | static_cast<uint32_t>(a);
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
std::pair<uint8_t, float> Style::parse_dash(const QJsonValue& dash_values)
{
    // currently only 2 values are allowed here
    // TODO qwant and osm-bright style with id "boundary-land-level-4" uses 4 values
    // allows the following dash pattern: - . - . -
    assert(dash_values.isArray());
    const auto dash_array = dash_values.toArray();
    if (dash_array.size() < 2) {
        return { 1 * constants::style_precision, 1 };
    }
    assert(dash_array.size() >= 2);

    const auto dashes = dash_array[0].toDouble();
    const auto gaps = dash_array[1].toDouble();

    const auto sum = (dashes + gaps);
    assert(sum > 0); // sum=0 is very bad -> but also shouldnt happen

    // ratio between dashes and gaps
    const uint8_t dash_gap_ratio = (dashes / sum) * constants::style_precision;

    return { dash_gap_ratio, sum * constants::dash_multiplier };
}

uint16_t Style::parse_line_width(const QJsonValue& value)
{
    if (value.isDouble()) {
        float thickness = value.toDouble() * constants::line_width_multiplier;
        if (thickness > constants::max_line_width)
            thickness = constants::max_line_width;
        return uint16_t(round(thickness * constants::style_precision));
    }

    qDebug() << "unhandled line width value" << value;
    assert(false);
    return 0;
}

} // namespace nucleus::vector_layer
