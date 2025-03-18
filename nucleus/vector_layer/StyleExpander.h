/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2025 Lucas Dworschak
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

#include <radix/hasher.h>

namespace nucleus::vector_layer::style_expander::details {

struct PairHasher {
    size_t operator()(const std::pair<QJsonValue, QJsonValue>& pair) const
    {
        size_t seed = 0;

        radix::hasher::hash_combine<size_t>(seed, qHash(pair.first));
        radix::hasher::hash_combine<size_t>(seed, qHash(pair.second));

        return seed;
    }
};

struct SubLayerInfo {
    // key is a hash that is generated from filter_criterium/filter_criterium value (e.g. class, landuse)
    // the key of the filter matches the key to a specific entry in the criterium_values map
    std::vector<size_t> filter_order;
    std::unordered_map<size_t, std::vector<QJsonArray>> filter;
    // the value is a qHash of the criterium_value (for match it is a combination of criterium and value)
    QMap<size_t, QSet<size_t>> criterium_values;
};

QJsonValue evaluate_get(QJsonValue expression);
bool sub_is_array(QJsonObject obj, QString sub_key);
QJsonValue get_match_value(QJsonArray match_array, QSet<size_t> criterium_hashes);
QJsonValue get_case_value(QJsonArray match_array, QSet<size_t> criterium_hashes);
SubLayerInfo generate_expanded_filters(const QJsonObject& paint, const QJsonArray& filter);
QJsonArray rejoin_filter(std::vector<QJsonArray> filters);

} // namespace nucleus::vector_layer::style_expander::details

namespace nucleus::vector_layer::style_expander {
QJsonArray expand(const QJsonArray& layers);
QJsonArray expand2(const QJsonArray& layers);
} // namespace nucleus::vector_layer::style_expander
