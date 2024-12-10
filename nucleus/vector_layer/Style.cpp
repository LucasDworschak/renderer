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

            // TODO here
            // convert string colors to int colors
            if (key == "fill-color") {
                s.fill_color = 0;
            } else if (key == "line-color") {
                s.fill_color = 0;
            } else if (key == "fill-outline-color") {
                s.outline_color = 0;
            } else if (key == "fill-outline-color") {
                s.outline_color = 0;
            } else if (key == "line-width") {
                s.outline_width = 0;
            } else if (key == "line-dasharray") {
                s.outline_dash = 0;
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
        const auto layer_name = obj.toObject().value("id").toString();

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

        // map layer name and zoom level to correct style hash
        for (unsigned i = min_zoom; i < max_zoom; i++) {
            m_layer_zoom_to_style[std::make_tuple(layer_name, i)] = style_hash;
        }
    }

    emit load_finished();
}

} // namespace nucleus::vector_layer
