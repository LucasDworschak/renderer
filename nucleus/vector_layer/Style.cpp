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

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QNetworkAccessManager>
#include <QNetworkReply>
#include <QtVersionChecks>

#include <iostream>

// #include <mapbox/vector_tile.hpp>
#include "nucleus/vector_tile/util.h"

namespace nucleus::vector_layer {

Style::Style(const QString& url)
    : m_url(url)
    , m_network_manager(new QNetworkAccessManager(this))
{
}

void Style::load()
{
    QNetworkRequest request(m_url);
    request.setTransferTimeout(m_transfer_timeout);
    request.setAttribute(QNetworkRequest::CacheLoadControlAttribute, QNetworkRequest::PreferCache);
#if QT_VERSION >= QT_VERSION_CHECK(6, 5, 0)
    request.setAttribute(QNetworkRequest::UseCredentialsAttribute, false);
#endif

    QNetworkReply* reply = m_network_manager->get(request);
    connect(reply, &QNetworkReply::finished, [reply, this]() {
        const auto error = reply->error();
        if (error == QNetworkReply::NoError) {
            auto tile = std::make_shared<QByteArray>(reply->readAll());
            parse_load(tile);
        } else if (error == QNetworkReply::ContentNotFoundError) {
            auto tile = std::make_shared<QByteArray>();
            parse_load(tile);
        } else {
            //            qDebug() << reply->url() << ": " << error;
            auto tile = std::make_shared<QByteArray>();
            parse_load(tile);
        }
        reply->deleteLater();
    });
}

unsigned int Style::transfer_timeout() const { return m_transfer_timeout; }

void Style::set_transfer_timeout(unsigned int new_transfer_timeout)
{
    assert(new_transfer_timeout < unsigned(std::numeric_limits<int>::max()));
    m_transfer_timeout = new_transfer_timeout;
}

size_t Style::layer_style_index(std::string layer_name, unsigned zoom) const
{
    const auto key = std::make_tuple(layer_name, zoom);

    // assert(m_layer_zoom_to_style.contains(key)); // no valid style found
    if (!m_layer_zoom_to_style.contains(key)) {
        // std::cout << "no style for: " << layer_name << " " << zoom << std::endl;
        return -1ul;
    }

    return m_layer_zoom_to_style.at(key);
}

void Style::parse_load(std::shared_ptr<QByteArray> data)
{

    if (data->isEmpty()) {
        std::cout << "" << std::endl;
        return;
    }

    const auto style_hasher = Hasher();

    QJsonDocument doc = QJsonDocument::fromJson(*data);
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
            // std::cout << "no style: " << obj.toObject().value("id").toString().toStdString() << std::endl;
            continue;
        }

        auto paint = obj.toObject().value("paint").toObject();

        Layer_Style s;

        for (const QString& key : paint.keys()) {

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
                std::cout << "new key detected: " << key.toStdString() << std::endl;
            }
        }

        unsigned min_zoom = 0u;
        unsigned max_zoom = 19u;

        if (obj.toObject().contains("minzoom"))
            min_zoom = obj.toObject().value("minzoom").toInt();
        if (obj.toObject().contains("maxzoom"))
            max_zoom = obj.toObject().value("maxzoom").toInt();

        const auto style_hash = style_hasher(s);

        // insert styles if they aren't inserted yet
        if (fill) {
            if (!m_fill_styles.contains(s)) {
                m_fill_styles.insert(s);
            }
        } else {
            if (!m_line_styles.contains(s)) {
                m_fill_styles.insert(s);
            }
        }

        const auto layer_name = obj.toObject().value("source-layer").toString().toStdString();
        // map layer name and zoom level to correct style hash
        for (unsigned i = min_zoom; i < max_zoom; i++) {
            m_layer_zoom_to_style[std::make_tuple(layer_name, i)] = style_hash;
        }
    }

    emit load_finished();
}

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
{
    if (value.starts_with("#")) {
        if (value.length() == 7)
            return (std::stoul(value.substr(1), nullptr, 16) << 8) | 255;
        else if (value.length() == 9)
            return std::stoul(value.substr(1), nullptr, 16);
        else {
            std::cout << "cannot parse color: " << value << std::endl; // TODO change to qdebug (or similar)
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
        std::cout << "cannot parse color: " << value << std::endl;
        return 0ul;
    }
}

} // namespace nucleus::vector_layer
