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

#include <unordered_set>

#include "nucleus/tile/constants.h"

class QNetworkAccessManager;

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

class Style : public QObject {
    Q_OBJECT
public:
    Style(const QString& url);

    [[nodiscard]] unsigned int transfer_timeout() const;
    void set_transfer_timeout(unsigned int new_transfer_timeout);

    uint32_t parse_color(std::string value);
    uint32_t parse_dasharray(QJsonArray dash_values);

public slots:
    void load();

signals:
    void load_finished();

private:
    unsigned m_transfer_timeout = tile::constants::default_network_timeout;

    std::unordered_set<Layer_Style, Hasher> m_fill_styles;
    std::unordered_set<Layer_Style, Hasher> m_line_styles;

    std::map<std::tuple<QString, unsigned>, size_t> m_layer_zoom_to_style;

    QString m_url;
    std::shared_ptr<QNetworkAccessManager> m_network_manager;

    void parse_load(std::shared_ptr<QByteArray> data);
};

} // namespace nucleus::vector_layer
