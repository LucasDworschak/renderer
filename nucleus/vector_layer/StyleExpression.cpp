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

#include <radix/hasher.h>

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
                qDebug() << "StyleExpression: " << data[1] << " is not an array";
                assert(false); // shouldnt happen -> sub expression is not an array
            }
        }
    } else if (StyleExpression::valid(data[0].toString())) {
        // we have an expression
        return std::make_unique<StyleExpression>(data);
    }

    qDebug() << "StyleExpression: " << data[0].toString() << " is not a valid operator";

    assert(false);
    return {};
}

size_t StyleExpression::hash()
{
    size_t seed = 0;

    radix::hasher::hash_combine<std::string>(seed, m_key);
    radix::hasher::hash_combine<int>(seed, static_cast<int>(m_comparator));

    for (const auto& v : m_values) {
        radix::hasher::hash_combine<int>(seed, v);
    }

    return seed;
}

size_t StyleExpressionCollection::hash()
{
    size_t seed = 0;

    radix::hasher::hash_combine<bool>(seed, m_negate);
    radix::hasher::hash_combine<bool>(seed, m_all);

    for (const auto& v : m_subFilters) {
        radix::hasher::hash_combine<size_t>(seed, v->hash());
    }

    return seed;
}
void StyleExpression::initialize()
{
    // qDebug() << "initialize style";
    string_value_map = std::unordered_map<std::string, int>();

    // initialize types
    string_value_map["Point"] = 1;
    string_value_map["LineString"] = 2;
    string_value_map["Polygon"] = 3;
    string_value_map["UNKNOWN"] = null_value;

    string_value_map[""] = null_value;

    counter = 4;
}

std::unordered_map<std::string, int> StyleExpression::get_values(const mapbox::vector_tile::feature& feature)
{
    const auto type = feature.getType();
    const auto properties = feature.getProperties();

    auto values = std::unordered_map<std::string, int>();

    if (type == mapbox::vector_tile::GeomType::POINT)
        values["$type"] = 1;
    else if (type == mapbox::vector_tile::GeomType::LINESTRING)
        values["$type"] = 2;
    else if (type == mapbox::vector_tile::GeomType::POLYGON)
        values["$type"] = 3;
    else
        values["$type"] = null_value;

    // TODO parse other properties

    for (const auto& prop : properties) {
        const auto key = prop.first;
        int value = null_value;

        if (std::holds_alternative<std::string>(prop.second)) {
            const auto string_value = std::get<std::string>(prop.second);
            if (string_value_map.contains(string_value)) {
                value = string_value_map[string_value];
            } else {
                // if (key == "class") {
                //     qDebug() << "class does not exist:" << string_value;
                // }
            }

            // else -> values is not in mqster string map -> just save null value
            // we probably can NOT just ignore this value since a "has" expression still must be valid
        } else if (std::holds_alternative<double>(prop.second)) {
            value = int(std::get<double>(prop.second) * float_precision);
        } else if (std::holds_alternative<uint64_t>(prop.second)) {
            value = int(std::get<uint64_t>(prop.second) * float_precision);
        } else if (std::holds_alternative<int64_t>(prop.second)) {
            value = int(std::get<int64_t>(prop.second) * float_precision);
        } else if (std::holds_alternative<bool>(prop.second)) {
            // some bools are defined as numbers in stylesheet -> multiply value by precision to be consistent
            value = int(std::get<bool>(prop.second)) * float_precision;
        }

        values[key] = value;
    }

    return values;
}

StyleExpression::StyleExpression(QJsonArray data)
    : StyleExpressionBase()
{
    const auto comparator = data[0].toString();
    m_negate = false;
    m_comparator_in = false;
    m_comparator_has = false;

    // key
    if (data[1].isArray()) {
        if (data[1].toArray().size() == 1) {
            m_key = data[1].toArray()[0].toString().toStdString();
        } else {
            qDebug() << "StyleExpression: single expression with array that is not correctly sized: " << data[1];
        }
    } else {
        m_key = data[1].toString().toStdString();
    }

    // values
    if (comparator == "in" || comparator == "!in") {
        m_comparator_in = true;
        m_negate = comparator.startsWith("!");
        // m_type = DataType::String;

        // multiple values
        if (data.size() == 3 && data[2].isArray()) {
            // we most likely have a "literal" sub array
            const auto values = extract_literal(data[2]).toArray();

            for (qsizetype i = 0; i < values.size(); ++i) {
                assert(values[i].isString());

                const auto string_value = values[i].toString().toStdString();
                if (string_value_map.contains(string_value)) {
                    const int int_value = string_value_map[string_value];
                    m_values.push_back(int_value);
                } else {
                    const int int_value = counter++;
                    string_value_map[string_value] = int_value;
                    m_values.push_back(int_value);
                }
            }

        } else {
            for (qsizetype i = 2; i < data.size(); ++i) {
                if (data[i].isDouble()) {
                    m_values.push_back(int(data[i].toDouble() * float_precision));
                } else if (data[i].isString()) {

                    const auto string_value = data[i].toString().toStdString();

                    if (string_value_map.contains(string_value)) {
                        const int int_value = string_value_map[string_value];
                        m_values.push_back(int_value);
                    } else {
                        const int int_value = counter++;
                        string_value_map[string_value] = int_value;
                        m_values.push_back(int_value);
                    }

                } else {
                    qDebug() << "StyleExpression: invalid data: " << data[i];
                    assert(false);
                }
            }
        }
    } else if (comparator == "has" || comparator == "!has") {
        m_comparator_has = true;
        m_negate = comparator.startsWith("!");
        // no values for has expressions -> we only check if key is present
        if (data.size() > 2) {
            qDebug() << "StyleExpression: \"has\" expression with too many arguments:" << data.size();
            assert(false);
        }
    } else {

        if (comparator == "==")
            m_comparator = Comparator::Equals;
        else if (comparator == "!=")
            m_comparator = Comparator::NotEquals;
        else if (comparator == "<")
            m_comparator = Comparator::less;
        else if (comparator == "<=")
            m_comparator = Comparator::lessThanEqual;
        else if (comparator == ">")
            m_comparator = Comparator::greater;
        else if (comparator == ">=")
            m_comparator = Comparator::greaterThanEqual;
        else {
            qDebug() << "StyleExpression: unknown comparator";
            assert(false);
        }

        // only one more value in data that contains the value to compare to.
        // Note: it might be possible that there is an additional entry for collator for "locale-dependent comparison operations"
        // but we want to make sure for now that our stylesheet doesnt contain this (maybe add later if needed)
        if (data.size() != 3) {
            qDebug() << "StyleExpression: data size too large: " << data.size();
            qDebug() << data[0] << data[1];
            assert(data.size() == 3);
        }

        if (data[2].type() == QJsonValue::Type::String) {
            const auto string_value = data[2].toString().toStdString();
            if (string_value_map.contains(string_value)) {
                const int int_value = string_value_map[string_value];
                m_values.push_back(int_value);
            } else {
                const int int_value = counter++;
                string_value_map[string_value] = int_value;
                m_values.push_back(int_value);
            }
        } else if (data[2].type() == QJsonValue::Type::Double) {
            m_values.push_back(int(data[2].toDouble() * float_precision));
        } else if (data[2].type() == QJsonValue::Type::Bool) {
            // some bools are defined as numbers in stylesheet -> multiply value by precision to be consistent
            m_values.push_back(int(data[2].toBool()) * float_precision);
        } else {
            qDebug() << "StyleExpression: style value type not defined! (should not happen)";
            assert(false);
        }
    }
}

bool StyleExpression::matches(const std::unordered_map<std::string, int>& value_map)
{
    int value = null_value;
    if (value_map.contains(m_key))
        value = value_map.at(m_key);

    if (m_comparator_in) {
        // a list -> string is sufficient
        const bool find_in = !m_negate; // if the comparator is "!in" we want to negate the values

        for (const auto& v : m_values) {
            if (v == value) {
                return find_in;
            }
        }
        return !find_in; // string wasnt in the list of acceptable values
    } else if (m_comparator_has) {
        // if the value is a null value -> the value was not in the list
        bool has_not_value = value == null_value;
        return m_negate ? has_not_value : !has_not_value;
    } else {

        if (m_comparator == Comparator::Equals)
            return value == m_values[0];
        else if (m_comparator == Comparator::NotEquals)
            return value != m_values[0];
        else if (m_comparator == Comparator::less)
            return value < m_values[0];
        else if (m_comparator == Comparator::lessThanEqual)
            return value <= m_values[0];
        else if (m_comparator == Comparator::greater)
            return value > m_values[0];
        else if (m_comparator == Comparator::greaterThanEqual)
            return value >= m_values[0];
        return false;
    }
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
                qDebug() << "StyleExpression: empty array";
                assert(false);
            }
            m_subFilters.push_back(create_filter_expression(data[i].toArray()));
        } else {
            qDebug() << "StyleExpression: not an array";
            assert(false);
        }
    }
}

bool StyleExpressionCollection::matches(const std::unordered_map<std::string, int>& value_map)
{
    if (m_negate) {
        assert(m_subFilters.size() == 1);
        // negate should only handle one subfilter -> get the match and negate it
        return !m_subFilters[0]->matches(value_map);
    }

    for (size_t i = 0; i < m_subFilters.size(); ++i) {
        assert(m_subFilters[i] != nullptr);
        bool result = m_subFilters[i]->matches(value_map);
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
