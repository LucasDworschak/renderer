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

#include "StyleExpander.h"

#include <unordered_set>

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

namespace nucleus::vector_layer::style_expander::details {

bool sub_is_array(QJsonObject obj, QString sub_key)
{
    if (!obj.contains(sub_key))
        return false; // sub_key not present
    return obj[sub_key].isArray();
}

QJsonValue get_match_value(QJsonArray match_array, QString match_key)
{
    if (match_array.contains(match_key)) {
        // look at every second value and find the key that matches
        for (qsizetype i = 2; i < match_array.size() - 1; i += 2) {
            if (match_array[i].isString()) {
                auto el = match_array[i].toString();
                if (el == match_key) {
                    // we found the key -> return the next value
                    return match_array[i + 1];
                }
            }
        }
    }

    // key not found -> return default value
    return match_array.last();
}

std::unordered_map<QString, std::vector<QJsonArray>> get_sub_layer(QJsonArray filter)
{
    std::unordered_map<QString, std::vector<QJsonArray>> sub_layer;

    std::vector<QJsonArray> sub_filter;

    int in_filter_index = -1;

    if (filter.contains("all")) {
        for (qsizetype i = 1; i < filter.size(); i++) {
            assert(filter[i].isArray()); // all sub filters should also be arrays

            // add filter to sub filters
            sub_filter.push_back(filter[i].toArray());

            if (sub_filter[i - 1].contains("in")) {
                assert(in_filter_index == -1); // 2 in conditions -> function needs to be rewritten do deal with it
                in_filter_index = i - 1;
            }
        }
    } else {
        // no all -> only save the current instance as the sub filter
        sub_filter.push_back(filter);

        if (filter.contains("in")) {
            in_filter_index = 0;
        }
    }

    if (in_filter_index == -1) {
        // no filter splitting necessary here
        sub_layer["all"] = sub_filter;
        return sub_layer;
    }

    const auto criterium = sub_filter[in_filter_index][1].toString();
    // get all layers that are stated in the in expression
    for (qsizetype i = 2; i < sub_filter[in_filter_index].size(); i++) {
        const auto layer_name = sub_filter[in_filter_index][i].toString();

        // replace the "in" filter with an "==" filter
        auto new_filter = std::vector<QJsonArray>(sub_filter);
        new_filter[in_filter_index] = QJsonArray { "==", criterium, layer_name };

        sub_layer[layer_name] = new_filter;
    }

    return sub_layer;
}

QJsonArray rejoin_filter(std::vector<QJsonArray> filters)
{
    assert(filters.size() > 0);

    if (filters.size() == 1) {
        // just one element -> we can just use it without changing anything
        return filters[0];
    }

    // we have to wrap everything in an all expression
    QJsonArray out_filter { "all" };
    for (size_t i = 0; i < filters.size(); i++) {
        out_filter.append(filters[i]);
    }
    return out_filter;
}

} // namespace nucleus::vector_layer::style_expander::details

namespace nucleus::vector_layer::style_expander {
QJsonArray expand(const QJsonArray& layers)
{
    QJsonArray out_layer;

    for (const auto& layer : layers) {

        const auto paint = layer.toObject().value("paint").toObject();

        // if any of the following values is true we will epxand the layer to individual ones
        // else we will simply copy the existing layer
        bool fill_color_match = details::sub_is_array(paint, "fill-color");
        bool fill_opacity_match = details::sub_is_array(paint, "fill-opacity");
        bool line_color_match = details::sub_is_array(paint, "line-color");
        bool line_width_match = details::sub_is_array(paint, "line-width");
        bool line_opacity_match = details::sub_is_array(paint, "line-opacity");

        bool needs_to_expand = fill_color_match || fill_opacity_match || line_color_match || line_width_match || line_opacity_match;

        if (!needs_to_expand) {
            out_layer.append(layer);
            continue;
        }

        // we need to expand the layer
        // qDebug() << layer.toObject().value("id").toString();
        const auto sublayer = details::get_sub_layer(layer.toObject().value("filter").toArray());

        // store criterium_key, criterium_value of all the paint expressions
        std::unordered_map<QString, std::unordered_set<QString>> paint_criterium_values;

        for (const auto& key : paint.keys()) {
            if (paint[key].isArray()) {
                QJsonArray paint_value_array = paint[key].toArray();
                if (paint_value_array[0] == "match") {
                    // asserts make sure that we have something like: ["get", value]
                    assert(paint_value_array[1].isArray());
                    assert(paint_value_array[1].toArray().size() == 2);
                    assert(paint_value_array[1].toArray()[0].isString());
                    assert(paint_value_array[1].toArray()[0].toString() == "get");
                    assert(paint_value_array[1].toArray()[1].isString());

                    QString criterium_key = paint_value_array[1].toArray()[1].toString();

                    for (qsizetype i = 2; i < paint_value_array.size() - 1; i += 2) {
                        if (paint_value_array[i].isString()) {
                            // add criterium_key, and criterium_value to list
                            paint_criterium_values[criterium_key].insert(paint_value_array[i].toString());
                        }
                    }

                } else if (paint_value_array[0] == "case") {
                    // TODO implement
                }
            }
        }

        for (const auto& [layer_name, new_filter] : sublayer) {
            QJsonObject new_paint;
            for (const auto& key : paint.keys()) {
                if (paint[key].isArray() && paint[key].toArray()[0] == "match") {
                    // choose the value that matches and use it for the paint option
                    new_paint[key] = QJsonValue(details::get_match_value(paint[key].toArray(), layer_name));
                } else if (paint[key].isArray() && paint[key].toArray()[0] == "case") {

                    //

                    // TODO implement
                    // new_paint[key] = QJsonValue(get_match_value(paint[key].toArray(), layer_name));
                    new_paint[key] = paint[key];
                } else {
                    // write the previous value
                    new_paint[key] = paint[key];
                }
            }

            // copy layer so that we can change values specific for one sublayer
            QJsonObject new_layer = QJsonObject(layer.toObject());
            new_layer["filter"] = details::rejoin_filter(new_filter);
            new_layer["paint"] = new_paint;

            out_layer.append(new_layer);
        }
    }

    return out_layer;
}

} // namespace nucleus::vector_layer::style_expander
