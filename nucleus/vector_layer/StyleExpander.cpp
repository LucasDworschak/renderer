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

#include "StyleExpander.h"

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

QJsonValue evaluate_get(QJsonValue expression)
{
    if (expression.isArray()) {
        const auto arr = expression.toArray();
        if (arr.size() >= 2 && arr[0].isString() && arr[0].toString() == "get") {
            return arr[1];
        }
    }

    // not a get expression -> return it without changing anything
    return expression;
}

QJsonValue get_match_value(QJsonArray match_array, QSet<size_t> criterium_hashes)
{
    assert(match_array.size() > 4); // we expect that a match expression contains at least 5 values ("match", criterium key, criterium value, value defaultValue)

    // get the criterium value for the criterium key (e.g. get "farmland" for key "subclass")
    // this criterium value must be found in the match_array otherwise we return the default value
    const auto criterium = evaluate_get(match_array[1]).toString();

    // look at every second value and find the criterium that matches
    for (qsizetype i = 2; i < match_array.size() - 1; i += 2) {
        if (match_array[i].isString()) {
            const auto current_criterium_hash = qHash(criterium + match_array[i].toString());
            if (criterium_hashes.contains(current_criterium_hash)) {
                // we found it -> return the next value
                return match_array[i + 1];
            }
        }
    }

    // key not found -> return default value
    return match_array.last();
}

QJsonValue get_case_value(QJsonArray case_array, QSet<size_t> criterium_hashes)
{
    assert(case_array.size() > 3); // we expect that a case expression contains at least 4 values ("case", criterium value, value, default value)

    // look at every second value and find the criterium that matches
    for (qsizetype i = 1; i < case_array.size() - 1; i += 2) {
        assert(case_array[i].isArray());
        assert(case_array[i].toArray().size() == 3);

        const auto current_criterium_hash = qHash(evaluate_get(case_array[i][1]).toString() + case_array[i][2].toString());

        if (criterium_hashes.contains(current_criterium_hash)) {
            // we found it -> return the next value
            return case_array[i + 1];
        }
    }

    // key not found -> return default value
    return case_array.last();
}

/*
 * basic explanation/example:
 * we have something along the following filter/paint expressions:
 *   filter
 *      subclass in [a,b,c]
 *   paint
 *       case
 *           bridge=1
 *           default value
 *       case
 *           intermittent=1
 *           default value
 *       case
 *           subclass=a
 *           subclass=b
 *           default value
 *
 * at the end we want to expand it into the following conditions:
 *  === new filters ===
 *  subclass a (+ default values)
 *  subclass b (+ default values)
 *  subclass c (+ default values)
 *  subclass a + bridge 1
 *  subclass b + bridge 1
 *  subclass c + bridge 1
 *  subclass a + intermittent 1
 *  subclass b + intermittent 1
 *  subclass c + intermittent 1
 *  subclass a + bridge 1 + intermittent 1
 *  subclass b + bridge 1 + intermittent 1
 *  subclass c + bridge 1 + intermittent 1
 *
 *  Note a feature with (subclass c + bridge 1 + intermittent 1)
 *  would also get the styles for:
 *      subclass c (+ default values)
 *      subclass c + bridge 1
 *      subclass c + intermittent 1
 *      subclass c + bridge 1 + intermittent 1
 *  but essentially only the last style (subclass c + bridge 1 + intermittent 1) will be used
 *
 *
 *
 * structure of expressions:
 *  match
 *   criterium key
 *   criterium value1
 *   value1
 *   criterium value2
 *   value2
 *   criterium value3
 *   value3
 *   defaultValue
 *
 *  case
 *   criterium value1
 *   value1
 *   criterium value2
 *   value2
 *   defaultValue
 */

SubLayerInfo generate_expanded_filters(const QJsonObject& paint, const QJsonArray& filter)
{
    SubLayerInfo info;

    ////////////////////////////////////////////
    // get a sub_filter vector (that separates filters if there is an all statement)
    // also stores the location in the vector where an "in" expression is located at
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
        // no all -> only save the current instance as the sub filter (only if there is a filter)
        if (!filter.isEmpty())
            sub_filter.push_back(filter);

        if (filter.contains("in")) {
            in_filter_index = 0;
        }
    }

    ////////////////////////////////////////////
    // create sub filters by splitting the "in" expression into individual "==" filters
    // or if not applicable, create one "default" filter
    const auto hasher = PairHasher();
    QJsonValue filter_criterium;
    if (in_filter_index != -1) {
        filter_criterium = evaluate_get(sub_filter[in_filter_index][1]);
        // get all layers that are stated in the in expression
        for (qsizetype i = 2; i < sub_filter[in_filter_index].size(); i++) {
            const auto criterium_value = sub_filter[in_filter_index][i];

            // replace the "in" filter with an "==" filter
            auto new_filter = std::vector<QJsonArray>(sub_filter);
            new_filter[in_filter_index] = QJsonArray { "==", filter_criterium, criterium_value };
            // new_filter.push_back(QJsonArray { "==", filter_criterium, criterium_value });

            const auto key = hasher({ filter_criterium, criterium_value });
            info.filter_order.push_back(key);
            info.filter[key] = new_filter;
            if (!info.criterium_values.contains(key)) {
                info.criterium_values[key] = QSet<size_t>();
            }

            info.criterium_values[key].insert(qHash(filter_criterium.toString() + criterium_value.toString()));
        }
    } else {
        // no in filter, just add the normal filter
        const auto key = qHash("default");
        info.filter_order.push_back(key);
        info.filter[key] = sub_filter;
        info.criterium_values[key] = QSet<size_t>(); // empty
    }
    const auto filter_criterium_hash = qHash(QJsonArray { "==", filter_criterium });

    ////////////////////////////////////////////
    // check "match" and "case" expressions and find out if we have to generate additional filters
    QMap<size_t, QSet<QJsonArray>> paint_filters;
    for (const auto& key : paint.keys()) {
        if (paint[key].isArray()) {
            QJsonArray paint_value_array = paint[key].toArray();
            if (paint_value_array[0] == "match") {
                auto criterium = evaluate_get(paint_value_array[1]);

                for (qsizetype i = 2; i < paint_value_array.size() - 1; i += 2) {
                    const auto filter = QJsonArray { "==", criterium, paint_value_array[i] };

                    const auto filter_hash = qHash(QJsonArray { filter[0], filter[1] });

                    if (filter_hash == filter_criterium_hash) {
                        // we already have such a filter defined through an "in" expression
                        // -> we only need to add the criterium values
                        const auto key = hasher({ filter_criterium, filter[2] });
                        info.criterium_values[key].insert(qHash(filter_criterium.toString() + filter[2].toString()));
                    } else {
                        // no such filter was found -> we add it to the paint_filters -> and create one filter with and one without this condition (for default values)
                        paint_filters[filter_hash].insert(filter);
                    }
                }
            } else if (paint_value_array[0] == "case") {

                for (qsizetype i = 1; i < paint_value_array.size() - 1; i += 2) {
                    assert(paint_value_array[i].isArray());
                    const auto filter = paint_value_array[i].toArray();
                    assert(filter.size() == 3); // filter does not look like we expected
                    const auto simple_filter = QJsonArray { filter[0], evaluate_get(filter[1]), filter[2] };

                    const auto filter_hash = qHash(QJsonArray { simple_filter[0], simple_filter[1] });

                    if (filter_hash == filter_criterium_hash) {
                        // we already have such a filter defined through an "in" expression
                        // -> we only need to add the criterium values
                        const auto key = hasher({ filter_criterium, filter[2] });
                        info.criterium_values[key].insert(qHash(filter_criterium.toString() + filter[2].toString()));
                    } else {
                        // no such filter was found -> we add it to the paint_filters -> and create one filter with and one without this condition (for default values)
                        paint_filters[filter_hash].insert(simple_filter);
                    }
                }
            }
        }
    }

    ////////////////////////////////////////////
    // for every criterium in paint_filters we have to split previous filters (one without and one with the paint filter)
    // this will allow "case" and "match" expressions (that weren't previous in the filter) to work
    for (const auto& paint_filter_criterium_hash : paint_filters.keys()) {

        // we save current iteration filter -> since we will add new filters to it
        // once we have added all filters for this criterium we redo the current_iteration_filter so that we can combine multiple criteriums
        auto current_iteration_filter = info.filter;

        for (const auto& paint_filter : paint_filters[paint_filter_criterium_hash]) {
            for (auto& [filter_hash, filter] : current_iteration_filter) {
                // generate new key from previous hash and current criterium/value hash
                size_t key = filter_hash;
                radix::hasher::hash_combine<size_t>(key, qHash(paint_filter));

                // create a new filter with the current criterium/value
                auto new_filter = std::vector<QJsonArray>(filter);
                new_filter.push_back(paint_filter);

                // add filter and criterium/value to SubLayerInfo
                info.filter_order.push_back(key);
                info.filter[key] = new_filter;
                if (!info.criterium_values.contains(key))
                    info.criterium_values[key] = QSet<size_t>();
                info.criterium_values[key].insert(qHash(paint_filter));
            }
        }
    }

    return info;
}

QJsonArray rejoin_filter(std::vector<QJsonArray> filters)
{
    if (filters.size() == 0) {
        return {};
    }

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

        ///////////////////////////////////
        // we need to expand the layer

        // qDebug() << layer.toObject().value("id");
        const auto sub_layer_info = details::generate_expanded_filters(paint, layer.toObject().value("filter").toArray());

        // loop over all the filters we just generated
        // use the sub_layer_key to get the criteriums for this specific filter
        for (const auto& sub_layer_key : sub_layer_info.filter_order) {
            QJsonObject new_paint;
            const auto criteriums = sub_layer_info.criterium_values.value(sub_layer_key);

            for (const auto& key : paint.keys()) {
                if (paint[key].isArray() && paint[key].toArray()[0] == "match") {
                    // choose the value that matches and use it for the paint option
                    new_paint[key] = details::get_match_value(paint[key].toArray(), criteriums);
                } else if (paint[key].isArray() && paint[key].toArray()[0] == "case") {
                    new_paint[key] = details::get_case_value(paint[key].toArray(), criteriums);
                } else {
                    // write the previous value
                    new_paint[key] = paint[key];
                }
            }

            // copy layer so that we can change values specific for one sublayer
            QJsonObject new_layer = QJsonObject(layer.toObject());
            new_layer["filter"] = details::rejoin_filter(sub_layer_info.filter.at(sub_layer_key));
            new_layer["paint"] = new_paint;

            out_layer.append(new_layer);
        }
    }

    return out_layer;
}

} // namespace nucleus::vector_layer::style_expander
