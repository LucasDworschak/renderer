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

QJsonValue StyleExpressionBase::extract_literal(QJsonValue expression)
{
    assert(expression.isArray());
    assert(expression.toArray()[0].isString());
    assert(expression.toArray()[0].toString() == "literal");
    assert(expression.toArray().size() == 2);

    return expression.toArray()[1];
}

std::unique_ptr<StyleExpressionBase> StyleExpressionBase::create_filter_expression(QJsonArray data)
{
    if (data.isEmpty())
        return {}; // no filter was supplied

    if (StyleExpressionCollection::valid(data[0].toString())) {

        if (data.size() > 2 || data[0].toString() == "!") // we have a collection
            return std::make_unique<StyleExpressionCollection>(data);
        else if (data.size() == 2) {
            // it says it is an expression, but in reality it only contains one element to evaluate
            if (data[1].isArray())
                return create_filter_expression(data[1].toArray());
            else {
                qDebug() << data[1] << " is not an array";
                // assert(false); // shouldnt happen -> sub expression is not an array
            }
        }
    } else if (StyleExpression::valid(data[0].toString())) {
        // we have an expression
        return std::make_unique<StyleExpression>(data);
    }

    qDebug() << data[0].toString() << " is not a valid operator ...";

    // assert(false);
    return {};
}

StyleExpression::StyleExpression(QJsonArray data)
    : StyleExpressionBase()
{
    m_comparator = data[0].toString().toStdString();

    // key
    if (data[1].isArray()) {
        if (data[1].toArray().size() == 1) {
            m_key = data[1].toArray()[0].toString().toStdString();
        } else {
            qDebug() << "single expression with array that is not correctly sized: " << data[1];
        }

    } else {
        m_key = data[1].toString().toStdString();
    }

    // values
    if (m_comparator == "in" || m_comparator == "!in") {
        // multiple values
        if (data.size() == 3 && data[2].isArray()) {
            // we most likely have a "literal" sub array
            const auto values = extract_literal(data[2]).toArray();

            for (qsizetype i = 0; i < values.size(); ++i) {
                assert(values[i].isString());
                m_values.push_back(values[i].toString().toStdString());
            }

        } else {
            for (qsizetype i = 2; i < data.size(); ++i) {
                if (data[i].isDouble()) {
                    m_values.push_back(std::to_string(data[i].toDouble()));
                } else if (data[i].isString()) {
                    m_values.push_back(data[i].toString().toStdString());
                } else {
                    qDebug() << "invalid data: " << data[i];
                    assert(false);
                }
            }
        }

    } else if (m_comparator == "has" || m_comparator == "!has") {
        // no values for has expressions -> we only check if key is present
        if (data.size() > 2) {
            qDebug() << "\"has\" expression with too many arguments:" << data.size();
            assert(false);
        }
    } else {
        // only one more value in data that contains the value to compare to.
        // Note: it might be possible that there is an additional entry for collator for "locale-dependent comparison operations"
        // but we want to make sure for now that our stylesheet doesnt contain this (maybe add later if needed)
        if (data.size() != 3) {
            qDebug() << "data size too large: " << data.size();
            qDebug() << data[0] << data[1];
            assert(data.size() == 3);
        }

        // even if the value is a number -> convert it into a string -> we are converting back later if needed
        // maybe not too ideal, but probably better if the json stores "long long" as a string and we wonder later why jsonValue.toInt returns 0 for long long numbers
        m_values.push_back(data[2].toVariant().toString().toStdString());
    }
}

bool compare_equal(std::string) { return true; }
bool compare_equal(int64_t) { return true; }

bool StyleExpression::matches(const mapbox::vector_tile::GeomType& type, const mapbox::feature::properties_type& properties)
{
    // qDebug() << "match: " << m_key;
    mapbox::feature::value value = extract_value(type, properties);

    if (m_comparator == "in" || m_comparator == "!in") {
        const auto string_value = std::visit(vector_tile::util::string_print_visitor, value).toStdString(); // we are only using the value to see if it is in a list -> string is sufficient
        const bool find_in = !m_comparator.starts_with("!"); // if the comparator is "!in" we want to negate the values

        for (const auto& v : m_values) {
            if (v == string_value) {
                return find_in;
            }
        }
        return !find_in; // string wasnt in the list of acceptable values

    } else if (m_comparator == "has" || m_comparator == "!has") {
        // if the value is a null value -> the value was not in the list
        bool has_not_value = std::visit([](auto&& v) { return typeid(v) == typeid(mapbox::feature::null_value); }, value);
        return m_comparator.starts_with("!") ? has_not_value : !has_not_value;
    } else {
        return std::visit([this](auto&& v) { return vector_tile::util::compare_visitor(v, m_values[0], m_comparator); }, value);
    }
}

mapbox::feature::value StyleExpression::extract_value(const mapbox::vector_tile::GeomType& type, const mapbox::feature::properties_type& properties)
{
    if (m_key == "geometry-type" || m_key == "$type") {
        if (type == mapbox::vector_tile::GeomType::LINESTRING)
            return "LineString";
        else if (type == mapbox::vector_tile::GeomType::POINT)
            return "Point";
        else if (type == mapbox::vector_tile::GeomType::POLYGON)
            return "Polygon";
        else
            return "UNKNOWN";
    }

    if (properties.contains(m_key))
        return properties.at(m_key);

    return mapbox::feature::null_value; // could not determine what to extract (or feature does not contain key)
}

bool StyleExpression::valid(QString value)
{
    if (value == "==" || value == "!=" || value == "<" || value == ">" || value == "<=" || value == ">=" || value == "in" || value == "!in" || value == "has" || value == "!has")
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
        if (data[i].isArray()) {
            if (data[i].toArray().isEmpty()) {
                qDebug() << "empty array";
                assert(false);
            }
            m_subFilters.push_back(create_filter_expression(data[i].toArray()));
        } else {
            qDebug() << "not an array??";
            assert(false);
        }
    }
}

bool StyleExpressionCollection::matches(const mapbox::vector_tile::GeomType& type, const mapbox::feature::properties_type& properties)
{
    if (m_negate) {
        assert(m_subFilters.size() == 1);
        // negate should only handle one subfilter -> get the match and negate it
        return !m_subFilters[0]->matches(type, properties);
    }

    for (size_t i = 0; i < m_subFilters.size(); ++i) {
        assert(m_subFilters[i] != nullptr);
        bool result = m_subFilters[i]->matches(type, properties);
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
