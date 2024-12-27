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

#include "StyleExpression.h"

#include <mapbox/vector_tile.hpp>

#include "nucleus/vector_tile/util.h"

namespace nucleus::vector_layer {

std::unique_ptr<StyleExpressionBase> StyleExpressionBase::create_filter_expression(QJsonArray data)
{
    if (data.isEmpty())
        return {}; // no filter was supplied

    if (StyleExpressionCollection::valid(data[0].toString())) {

        if (data.size() > 2 || data[0].toString() == "!") // we have a collection
            return std::make_unique<StyleExpressionCollection>(data);
        else // it says it is an expression, but in reality it only contains one element to evaluate
            if (data[2].isArray())
                return create_filter_expression(data[2].toArray());
            else {
                qDebug() << data[2] << " is not an array";
                // assert(false); // shouldnt happen -> sub expression is not an array
            }
    } else if (StyleExpression::valid(data[0].toString())) {
        // we have an expression
        return std::make_unique<StyleExpression>(data);
    }

    qDebug() << data[0].toString() << " is not a valid operator";

    // assert(false);
    return {};
}

StyleExpression::StyleExpression(QJsonArray data)
    : StyleExpressionBase()
{
    m_comparator = data[0].toString();

    // key
    if (data[1].isArray()) {
        if (data[1].toArray().size() == 1) {
            m_key = data[1].toArray()[0].toString();
        } else {
            qDebug() << "single expression with array that is not correctly sized: " << data[1];
        }
    } else {
        m_key = data[1].toString();
    }

    // values
    if (m_comparator == "in") {
        // multiple values
        for (qsizetype i = 2; i < data.size(); ++i) {
            m_values.push_back(data[i].toString());
        }
    } else {
        // only one more value in data that contains the value to compare to.
        // Note: it might be possible that there is an additional entry for collator for "locale-dependent comparison operations"
        // but we want to make sure for now that our stylesheet doesnt contain this (maybe add later if needed)
        if (data.size() != 3) {
            qDebug() << "data size too large: " << data.size();
            assert(data.size() == 3);
        }

        // even if the value is a number -> convert it into a string -> we are converting back later if needed
        // maybe not too ideal, but probably better if the json stores "long long" as a string and we wonder later why jsonValue.toInt returns 0 for long long numbers
        m_values.push_back(data[2].toVariant().toString());
    }
}

bool compare_equal(std::string) { return true; }
bool compare_equal(int64_t) { return true; }

bool StyleExpression::matches(const mapbox::vector_tile::feature& feature)
{
    mapbox::feature::value value = extract_value(feature);
    // qDebug() << "match: " << m_key;

    if (m_comparator == "in") {
        auto string_value = std::visit(vector_tile::util::string_print_visitor, value); // we are only using the value to see if it is in a list -> string is sufficient

        for (const auto& v : m_values) {
            if (v == string_value) {
                return true;
            }
        }
        return false; // string wasnt in the list of acceptable values
    } else {

        return std::visit([this](auto&& v) { return vector_tile::util::compare_visitor(v, m_values[0], m_comparator); }, value);
    }
}

mapbox::feature::value StyleExpression::extract_value(const mapbox::vector_tile::feature& feature)
{
    if (m_key == "geometry-type" || m_key == "$type") {
        if (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING)
            return "LineString";
        else if (feature.getType() == mapbox::vector_tile::GeomType::POINT)
            return "Point";
        else if (feature.getType() == mapbox::vector_tile::GeomType::POINT)
            return "Polygon";
        else
            return "UNKNOWN";
    }

    // TODO class and subclass have to be supplied by the vectortile server as simple properties
    // alternatively we have to redo the stylesheet to use values our vector tile server can easily serve (there should be utilities to try to convert the style)
    if (feature.getProperties().contains(m_key.toStdString()))
        return feature.getProperties()[m_key.toStdString()];

    return mapbox::feature::null_value; // could not determine what to extract (or feature does not contain key)
}

bool StyleExpression::valid(QString value)
{
    if (value == "==" || value == "!=" || value == "<" || value == ">" || value == "<=" || value == ">=" || value == "in")
        return true;
    return false;
}

StyleExpressionCollection::StyleExpressionCollection(QJsonArray data)
    : StyleExpressionBase()
{
    m_all = data[0] == "all";
    m_negate = data[0] == "!";

    m_subFilters = std::vector<std::unique_ptr<StyleExpressionBase>>();

    for (qsizetype i = 1; i < data.size(); ++i) {
        if (data[i].isArray())
            m_subFilters.push_back(create_filter_expression(data[i].toArray()));
        else
            qDebug() << "not an array??";
    }
}

bool StyleExpressionCollection::matches(const mapbox::vector_tile::feature& feature)
{
    if (m_negate) {
        // negate should only handle one subfilter -> get the match and negate it
        return !m_subFilters[0]->matches(feature);
    }

    for (size_t i = 0; i < m_subFilters.size(); ++i) {
        bool result = m_subFilters[i]->matches(feature);
        if (m_all && !result) {
            return false; // the current match was false but all had to be true
        } else if (!m_all && result) {
            return true; // the current result was true, and at least one had to be true (no need to test the others)
        }
    }

    // all options have been checked
    if (m_all)
        return true; // all had been true

    // not one option was true
    return false;
}

bool StyleExpressionCollection::valid(QString value)
{
    if (value == "all" || value == "any" || value == "!")
        return true;
    return false;
}

} // namespace nucleus::vector_layer
