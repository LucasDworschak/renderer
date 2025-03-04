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

#include <regex>

namespace nucleus::vector_layer {

// https://maplibre.org/maplibre-style-spec/

Style::Style(const QString& filename)
    : m_filename(filename)
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

    std::unordered_map<LayerStyle, uint32_t, Hasher> style_to_index;

    QJsonDocument doc = QJsonDocument::fromJson(data);
    QJsonArray layers = style_expander::expand(doc.object().value("layers").toArray());

    uint32_t layer_index = 0;

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

        auto paintObject = obj.toObject().value("paint").toObject();
        auto filterData = obj.toObject().value("filter").toArray();
        const auto layer_name = obj.toObject().value("source-layer").toString().toStdString() + "_" + obj.toObject().value("type").toString().toStdString();

        std::vector<std::pair<uint8_t, uint32_t>> fill_colors;
        std::vector<std::pair<uint8_t, uint32_t>> outline_colors;
        std::vector<std::pair<uint8_t, uint16_t>> widths;
        std::vector<std::pair<uint8_t, std::pair<uint8_t, uint8_t>>> dashes;
        std::vector<std::pair<uint8_t, uint8_t>> opacities;

        std::shared_ptr<StyleExpressionBase> filter = StyleExpressionBase::create_filter_expression(filterData);

        bool invalid = false;

        // uint8_t opacity = 255u;

        for (const QString& key : paintObject.keys()) {

            // fill-antialias ?
            if (key == "fill-color") {
                parse_colors(paintObject.value(key), fill_colors);
            } else if (key == "line-color") {
                parse_colors(paintObject.value(key), fill_colors);
            } else if (key == "fill-opacity" || key == "line-opacity") {
                parse_opacities(paintObject.value(key), opacities);
            } else if (key == "fill-outline-color") { // TODO visualize (only used for buildings)
                parse_colors(paintObject.value(key), outline_colors);
            } else if (key == "line-width") { //  "line-gap-width"
                parse_line_widths(paintObject.value(key), widths);
            } else if (key == "line-dasharray") {
                parse_dashes(paintObject.value(key).toArray(), dashes);
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

        // determine zoom range defined in style.json (or fall back to [8-20] range
        glm::uvec2 zoom_range(8u, 20u);
        if (obj.toObject().contains("minzoom"))
            zoom_range.x = obj.toObject().value("minzoom").toInt();
        if (obj.toObject().contains("maxzoom"))
            zoom_range.y = obj.toObject().value("maxzoom").toInt();

        { // determine style for every zoom level within range
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

            for (unsigned zoom = zoom_range.x; zoom < zoom_range.y + 1; zoom++) {

                // determine if we have to change current / prev values
                if (fill_colors_current_value.first < zoom) {
                    fill_colors_previous_value = fill_colors_current_value;
                    fill_colors_current_value = fill_colors[++fill_colors_index];
                }
                if (outline_colors_current_value.first < zoom) {
                    outline_colors_previous_value = outline_colors_current_value;
                    outline_colors_current_value = outline_colors[++outline_colors_index];
                }
                if (widths_current_value.first < zoom) {
                    widths_previous_value = widths_current_value;
                    widths_current_value = widths[++widths_index];
                }
                if (dashes_current_value.first < zoom) {
                    dashes_previous_value = dashes_current_value;
                    dashes_current_value = dashes[++dashes_index];
                }
                if (opacities_current_value.first < zoom) {
                    opacities_previous_value = opacities_current_value;
                    opacities_current_value = opacities[++opacities_index];
                }

                // interpolate between previous and current values
                auto fill_color = interpolate<uint32_t>(zoom, fill_colors_previous_value, fill_colors_current_value);
                auto outline_color = interpolate<uint32_t>(zoom, outline_colors_previous_value, outline_colors_current_value);
                auto width = interpolate<uint16_t>(zoom, widths_previous_value, widths_current_value);
                std::pair<uint8_t, uint8_t> dash = interpolate(zoom, dashes_previous_value, dashes_current_value);
                auto opacity = interpolate<uint8_t>(zoom, opacities_previous_value, opacities_current_value);

                // merge opacity with colors
                if (opacity != 255u) {
                    // if opacity is set it overrides any opacity from the color

                    fill_color &= 4294967040u; // bit mask that zeros out the opacity bits
                    fill_color |= opacity;

                    outline_color &= 4294967040u; // bit mask that zeros out the opacity bits
                    outline_color |= opacity;
                }

                const uint16_t merged_dash = (dash.first << 8) | dash.second;

                LayerStyle s { fill_color, outline_color, width, merged_dash };

                uint32_t style_index = -1u;
                // insert styles if they aren't inserted yet
                // and get the index of the style we want to use
                if (!style_to_index.contains(s)) {
                    style_index = style_values.size();
                    style_to_index[s] = style_index;

                    // add the style to the raster
                    //  // float to uint32 representation // TODO this does not work for release build
                    style_values.push_back({ s.fill_color, s.outline_color, s.outline_width, s.outline_dash });
                } else {
                    style_index = style_to_index[s];
                }

                if (!m_layer_to_style.contains(layer_name))
                    m_layer_to_style[layer_name] = StyleFilter();

                // TODO somehow store both this value and the next style_index value so that we can interpolate between both

                m_layer_to_style[layer_name].add_filter(style_index, layer_index, filter, zoom);

                // auto id = obj.toObject().value("id").toString(); // DEBUG -> what layers with what index are used
                // qDebug() << style_index << id;
            }
        }

        layer_index++;
    }

    // make sure that the style values are within the buffer size; resize them to this size and create the raster images
    assert(style_values.size() <= constants::style_buffer_size * constants::style_buffer_size);
    style_values.resize(constants::style_buffer_size * constants::style_buffer_size, glm::u32vec4(-1u));

    m_styles = std::make_shared<const nucleus::Raster<glm::u32vec4>>(nucleus::Raster<glm::u32vec4>(constants::style_buffer_size, std::move(style_values)));

    emit load_finished(m_styles);
}

std::vector<std::pair<uint32_t, uint32_t>> Style::indices(
    std::string layer_name, std::string type, unsigned zoom, const mapbox::vector_tile::feature& feature) const // + properties // type=Polygon/Point/... // class? subclass
{
    const auto layer = layer_name + "_" + type;

    if (!m_layer_to_style.contains(layer)) {
        // qDebug() << "no style for: " << layer_name;
        return {};
    }

    return m_layer_to_style.at(layer).indices(zoom, feature);
}

std::shared_ptr<const nucleus::Raster<glm::u32vec4>> Style::styles() const { return m_styles; }

// NOTE std::stof uses the locale to convert strings
// locale might be german and it expects a "," decimal
// we howewer want to force an english "." decimal point
// we could also temporarily change the locale and change it back afterwards
// but changing the locale might cause performance problems or other unexpected problems
// snippet from: https://stackoverflow.com/a/78993592 -> also talks about possible multi thread issues with changing locale
float stringToFloat(const std::string& value)
{
    auto index = value.find(".");
    if (index == std::string::npos) {
        return std::stoi(value);
    }
    int full = std::stoi(value.substr(0, index));
    int decimals = std::stoi(value.substr(index + 1));
    return full + double(decimals / pow(10, value.substr(index + 1).size()));
}

// uses https://github.com/maplibre/maplibre-style-spec/blob/main/src/expression/definitions/interpolate.ts -> exponentialInterpolation()
template <typename T>
T Style::interpolate(uint8_t zoom, std::pair<uint8_t, T> prev, std::pair<uint8_t, T> current)
{
    const auto diff = current.first - prev.first;

    if (diff == 0) // both values are the same
        return current.second;

    const auto progress = zoom - diff;

    float t = float(progress) / float(diff);

    return prev.second * (1.0 - t) + current.second * t;
}

// interpolate with pairs (for dashes)
std::pair<uint8_t, uint8_t> Style::interpolate(uint8_t zoom, std::pair<uint8_t, std::pair<uint8_t, uint8_t>> prev, std::pair<uint8_t, std::pair<uint8_t, uint8_t>> current)
{
    const auto diff = current.first - prev.first;
    const auto progress = zoom - diff;

    float t = float(progress) / float(diff);

    return { prev.second.first * (1.0 - t) + current.second.first * t, prev.second.second * (1.0 - t) + current.second.second * t };
}

void Style::parse_colors(const QJsonValue& value, std::vector<std::pair<uint8_t, uint32_t>>& colors)
{
    assert(colors.size() == 0); // make sure that we only parse attribute once
    if (value.isObject() && value.toObject().contains("stops")) {
        auto v = value.toObject().value("stops").toArray().last().toArray().last();
        // TODO iterate over all values
        colors.push_back({ 255, parse_color(v) });
    } else {
        // only one value -> add element with highest zoom value
        colors.push_back({ 255, parse_color(value) });
    }
}
void Style::parse_opacities(const QJsonValue& value, std::vector<std::pair<uint8_t, uint8_t>>& opacities)
{
    assert(opacities.size() == 0); // make sure that we only parse attribute once
    if (value.isObject() && value.toObject().contains("stops")) {
        auto v = value.toObject().value("stops").toArray().last().toArray().last();
        // TODO iterate over all values
        opacities.push_back({ 255, parse_opacity(v) });
    } else {
        // only one value -> add element with highest zoom value
        opacities.push_back({ 255, parse_opacity(value) });
    }
}
void Style::parse_dashes(const QJsonValue& value, std::vector<std::pair<uint8_t, std::pair<uint8_t, uint8_t>>>& dashes)
{
    assert(dashes.size() == 0); // make sure that we only parse attribute once
    if (value.isObject() && value.toObject().contains("stops")) {
        auto v = value.toObject().value("stops").toArray().last().toArray().last();
        // TODO iterate over all values
        dashes.push_back({ 255, parse_dash(v) });
    } else {
        // only one value -> add element with highest zoom value
        dashes.push_back({ 255, parse_dash(value) });
    }
}
void Style::parse_line_widths(const QJsonValue& value, std::vector<std::pair<uint8_t, uint16_t>>& widths)
{
    assert(widths.size() == 0); // make sure that we only parse attribute once
    if (value.isObject() && value.toObject().contains("stops")) {
        auto v = value.toObject().value("stops").toArray().last().toArray().last();
        // TODO iterate over all values
        widths.push_back({ 255, parse_line_width(v) });
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
        return uint16_t(value.toDouble() * constants::style_precision);
    }

    qDebug() << "unhandled line width value" << value;
    assert(false);
    return 0;
}

} // namespace nucleus::vector_layer
