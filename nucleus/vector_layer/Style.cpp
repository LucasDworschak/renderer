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

#include "nucleus/vector_layer/StyleExpression.h"
#include "nucleus/vector_layer/constants.h"

namespace nucleus::vector_layer {

Style::Style(const QString& filename)
    : m_filename(filename)
{
}

void Style::load()
{
    QFile file(m_filename);
    const auto open = file.open(QIODeviceBase::OpenModeFlag::ReadOnly);
    assert(open);

    const auto data = file.readAll();

    if (data.isEmpty()) {
        qDebug() << "no data";
        return;
    }

    std::vector<uint32_t> fill_style_values;
    std::vector<uint32_t> line_style_values;

    std::unordered_map<LayerStyle, uint32_t, Hasher> style_to_fill_index;
    std::unordered_map<LayerStyle, uint32_t, Hasher> style_to_line_index;

    QJsonDocument doc = QJsonDocument::fromJson(data);
    QJsonArray layers = doc.object().value("layers").toArray();

    for (const QJsonValue& obj : layers) {

        bool fill = true;

        if (obj.toObject().value("type").toString() == "line") {
            fill = false;
        } else if (obj.toObject().value("type").toString() == "fill") {
            // valid but nothing to do
        } else {
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

        auto paint = obj.toObject().value("paint").toObject();
        auto filterData = obj.toObject().value("filter").toArray();

        LayerStyle s;
        std::shared_ptr<StyleExpressionBase> filter = StyleExpressionBase::create_filter_expression(filterData);

        for (const QString& key : paint.keys()) {

            // TODO support fill-opacity
            // fill-antialias ?
            if (key == "fill-color") {
                s.fill_color = parse_color(obj.toObject().value("paint").toObject().value(key).toString().toStdString());
            } else if (key == "line-color") {
                s.fill_color = parse_color(obj.toObject().value("paint").toObject().value(key).toString().toStdString());
            } else if (key == "fill-outline-color") {
                s.outline_color = parse_color(obj.toObject().value("paint").toObject().value(key).toString().toStdString());
            } else if (key == "fill-outline-color") {
                s.outline_color = parse_color(obj.toObject().value("paint").toObject().value(key).toString().toStdString());
            } else if (key == "line-width") {
                s.outline_width = obj.toObject().value("paint").toObject().value(key).toDouble();
            } else if (key == "line-dasharray") {
                s.outline_dash = parse_dasharray(obj.toObject().value("paint").toObject().value(key).toArray());
            } else if (key == "line-offset") {
                // might be needed
            } else if (key == "icon-color" || key == "circle-color" || key == "circle-radius" || key == "circle-stroke-color" || key == "circle-stroke-width" || key == "text-color"
                || key == "text-halo-color" || key == "text-halo-width") {
                // not used
            } else {
                qDebug() << "new key detected: " << key.toStdString();
            }
        }

        uint32_t style_index = -1u;

        // insert styles if they aren't inserted yet
        // and get the index of the style we want to use
        if (fill) {
            if (!style_to_fill_index.contains(s)) {
                style_index = fill_style_values.size() / constants::style_data_size;
                style_to_fill_index[s] = style_index;

                // add the style to the raster
                fill_style_values.push_back(s.fill_color);
                fill_style_values.push_back(s.outline_color);
                fill_style_values.push_back(*reinterpret_cast<uint32_t*>(&s.outline_width)); // float to uint32 representation
                fill_style_values.push_back(s.outline_dash);
            } else {
                style_index = style_to_fill_index[s];
            }
        } else {
            if (!style_to_line_index.contains(s)) {
                style_index = line_style_values.size() / constants::style_data_size;
                style_to_line_index[s] = style_index;

                // add the style to the raster
                line_style_values.push_back(s.fill_color);
                line_style_values.push_back(s.outline_color);
                line_style_values.push_back(*reinterpret_cast<uint32_t*>(&s.outline_width)); // float to uint32 representation
                line_style_values.push_back(s.outline_dash);
            } else {
                style_index = style_to_line_index[s];
            }
        }

        glm::uvec2 zoom_range(0u, 19u);
        if (obj.toObject().contains("minzoom"))
            zoom_range.x = obj.toObject().value("minzoom").toInt();
        if (obj.toObject().contains("maxzoom"))
            zoom_range.y = obj.toObject().value("maxzoom").toInt();

        const auto layer_name = obj.toObject().value("source-layer").toString().toStdString();
        if (!m_layer_to_style.contains(layer_name))
            m_layer_to_style[layer_name] = StyleFilter();
        m_layer_to_style[layer_name].add_filter(style_index, filter, zoom_range);
    }

    // make sure that the style values are within the buffer size; resize them to this size and create the raster images
    assert(fill_style_values.size() <= constants::style_buffer_size * constants::style_buffer_size);
    assert(line_style_values.size() <= constants::style_buffer_size * constants::style_buffer_size);
    fill_style_values.resize(constants::style_buffer_size * constants::style_buffer_size, -1u);
    line_style_values.resize(constants::style_buffer_size * constants::style_buffer_size, -1u);

    m_styles.fill_styles = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::style_buffer_size, std::move(fill_style_values)));
    m_styles.line_styles = std::make_shared<const nucleus::Raster<uint32_t>>(nucleus::Raster<uint32_t>(constants::style_buffer_size, std::move(line_style_values)));

    emit load_finished();
}

uint32_t Style::layer_style_index(std::string layer_name, unsigned zoom, const mapbox::vector_tile::feature& feature) const // + properties // type=Polygon/Point/... // class? subclass
{
    if (!m_layer_to_style.contains(layer_name)) {
        qDebug() << "no style for: " << layer_name;
        return -1u;
    }

    return m_layer_to_style.at(layer_name).style_index(zoom, feature);
}

StyleBufferHolder Style::style_buffer() const { return m_styles; }

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
uint32_t Style::parse_dasharray(QJsonArray dash_values)
{
    // TODO there are also values with more than two values
    // currently only 2 values are allowed here
    // assert(dash_values.size() == 2);

    double sum = dash_values[0].toDouble() + dash_values[1].toDouble();

    uint16_t dash_gap_ratio = dash_values[0].toDouble() / sum * 65535.f;

    return dash_gap_ratio << 16 | uint16_t(sum);
}

uint32_t Style::parse_color(std::string value)
{ // TODO support hsl
    if (value.starts_with("#")) {
        if (value.length() == 7)
            return (std::stoul(value.substr(1), nullptr, 16) << 8) | 255;
        else if (value.length() == 9)
            return std::stoul(value.substr(1), nullptr, 16);
        else {
            qDebug() << "cannot parse color: " << value; // TODO change to qdebug (or similar)
            return 0ul;
        }
    } else if (value.starts_with("rgb")) {
        // parses rgb(int,int,int) and rgba(int,int,int,float)
        // ints in range [0-255]; float in range [0-1]
        uint32_t out = 0u;

        auto startPos = value.find("(");
        auto tmp = value.substr(startPos + 1, value.size() - startPos - 2);
        int count = 0;
        auto pos = tmp.find(',');
        while (pos != std::string::npos) {
            count++;
            out = out << 8;

            // find start of next digit or use the full value for the rest
            auto nextPos = tmp.find(',');
            auto nextVal = tmp;
            if (nextPos != std::string::npos)
                nextVal = tmp.substr(0, nextPos);

            // count < 4 necessary since the alpha value might be 1 -> and has to be multiplied by 255
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
    } else {
        qDebug() << "cannot parse color: " << value;
        return 0ul;
    }
}

} // namespace nucleus::vector_layer
