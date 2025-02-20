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

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

namespace nucleus::vector_layer::style_expander {
namespace nucleus::vector_layer::style_expander::details {

    bool sub_is_array(QJsonObject obj, QString sub_key);
    QJsonValue get_match_value(QJsonArray match_array, QString match_key);
    std::unordered_map<QString, std::vector<QJsonArray>> get_sub_layer(QJsonArray filter);
    QJsonArray rejoin_filter(std::vector<QJsonArray> filters);

} // namespace nucleus::vector_layer::style_expander::details

QJsonArray expand(const QJsonArray& layers);
} // namespace nucleus::vector_layer::style_expander
