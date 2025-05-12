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
                qDebug() << data[1] << " is not an array";
                assert(false); // shouldnt happen -> sub expression is not an array
            }
        }
    } else if (StyleExpression::valid(data[0].toString())) {
        // we have an expression
        return std::make_unique<StyleExpression>(data);
    }

    qDebug() << data[0].toString() << " is not a valid operator ...";

    assert(false);
    return {};
}

size_t StyleExpression::hash()
{
    size_t seed = 0;

    radix::hasher::hash_combine<std::string>(seed, m_key);
    radix::hasher::hash_combine<int>(seed, static_cast<int>(m_comparator));

    for (const auto& v : m_values) {
        if (m_type == DataType::String) {
            radix::hasher::hash_combine<std::string>(seed, std::get<std::string>(v));
        } else if (m_type == DataType::Number) {

            if (std::holds_alternative<double>(v)) {
                radix::hasher::hash_combine<double>(seed, std::get<double>(v));
            } else if (std::holds_alternative<uint64_t>(v)) {
                radix::hasher::hash_combine<uint64_t>(seed, std::get<uint64_t>(v));
            } else if (std::holds_alternative<int64_t>(v)) {
                radix::hasher::hash_combine<int64_t>(seed, std::get<int64_t>(v));
            }

        } else if (m_type == DataType::Bool) {
            radix::hasher::hash_combine<bool>(seed, std::get<bool>(v));
        }
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

StyleExpression::StyleExpression(QJsonArray data)
    : StyleExpressionBase()
{
    const auto comparator = data[0].toString();
    m_negate = false;
    m_comparator_in = false;
    m_comparator_has = false;
    m_type = DataType::Undefined;

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
                m_values.push_back(values[i].toString().toStdString());
                assert(m_type == DataType::Undefined || m_type == DataType::String);
                m_type = DataType::String;
            }

        } else {
            for (qsizetype i = 2; i < data.size(); ++i) {
                if (data[i].isDouble()) {
                    m_values.push_back(data[i].toDouble());
                    assert(m_type == DataType::Undefined || m_type == DataType::Number);
                    m_type = DataType::Number;
                } else if (data[i].isString()) {
                    m_values.push_back(data[i].toString().toStdString());
                    assert(m_type == DataType::Undefined || m_type == DataType::String);
                    m_type = DataType::String;
                } else {
                    qDebug() << "invalid data: " << data[i];
                    assert(false);
                }
            }
        }
    } else if (comparator == "has" || comparator == "!has") {
        m_comparator_has = true;
        m_negate = comparator.startsWith("!");
        // no values for has expressions -> we only check if key is present
        if (data.size() > 2) {
            qDebug() << "\"has\" expression with too many arguments:" << data.size();
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
            qDebug() << "Unknown comparator";
            assert(false);
        }

        // only one more value in data that contains the value to compare to.
        // Note: it might be possible that there is an additional entry for collator for "locale-dependent comparison operations"
        // but we want to make sure for now that our stylesheet doesnt contain this (maybe add later if needed)
        if (data.size() != 3) {
            qDebug() << "data size too large: " << data.size();
            qDebug() << data[0] << data[1];
            assert(data.size() == 3);
        }

        // even if the value is a number -> convert it into a string -> we are converting back later if needed
        // maybe not too ideal, but probably better if the json stores "long long" as a string and we wonder later why jsonValue.toInt returns 0 for long long
        // numbers m_values.push_back(data[2].toVariant().toString().toStdString());

        // enum Type {
        //     Null =  0x0,
        //     Bool = 0x1,
        //     Double = 0x2,
        //     String = 0x3,
        //     Array = 0x4,
        //     Object = 0x5,
        //     Undefined = 0x80
        // };

        switch (data[2].type()) {
        case QJsonValue::Type::String:
            m_values.push_back(data[2].toString().toStdString());
            m_type = DataType::String;
            break;
        case QJsonValue::Type::Double:
            m_values.push_back(data[2].toDouble());
            m_type = DataType::Number;
            break;

        // case QJsonValue::Type::INT:
        //     m_values.push_back( data[2].get_int64();
        //     break;
        // case QJsonValue::Type::UINT:
        //     m_values.push_back( data[2].get_uint64();
        //     break;
        // case QJsonValue::Type::SINT:
        //     m_values.push_back( data[2].get_sint64();
        //     break;
        case QJsonValue::Type::Bool:
            m_values.push_back(data[2].toBool());
            m_type = DataType::Bool;
            break;
        default:
            qDebug() << "style value type not defined! (should not happen)";
            assert(false);
        }
    }
}

bool StyleExpression::compare_values(const mapbox::feature::value& expression_value, const mapbox::feature::value& current_value, const Comparator& comparator)
{

    if (m_type == DataType::String) {

        if (const auto current_string_value = std::get_if<std::string>(&current_value)) {
            const auto expression_string_value = std::get<std::string>(expression_value);

            if (comparator == Comparator::Equals) {
                return expression_string_value == *current_string_value;
            } else if (comparator == Comparator::NotEquals) {
                return expression_string_value != *current_string_value;
            }
        } else {

            // current string is not a string (most likely null) -> we cannot compare -> return false if equals; true if not equals
            return comparator == Comparator::NotEquals;
        }
    } else if (m_type == DataType::Number) {

        double current_number_value = 0.0;
        if (std::holds_alternative<double>(current_value)) {
            current_number_value = std::get<double>(current_value);
        } else if (std::holds_alternative<uint64_t>(current_value)) {
            current_number_value = double(std::get<uint64_t>(current_value));
        } else if (std::holds_alternative<int64_t>(current_value)) {
            current_number_value = double(std::get<int64_t>(current_value));
        } else if (std::holds_alternative<mapbox::feature::null_value_t>(current_value)) {
            current_number_value = 0;
        } else {
            qDebug() << "could not parse number" << current_value.index() << m_key;
            assert(false);
        }

        if (const auto expression_number_value = std::get_if<double>(&expression_value)) {

            if (comparator == Comparator::Equals)
                return current_number_value == *expression_number_value;
            else if (comparator == Comparator::NotEquals)
                return current_number_value != *expression_number_value;
            else if (comparator == Comparator::less)
                return current_number_value < *expression_number_value;
            else if (comparator == Comparator::lessThanEqual)
                return current_number_value <= *expression_number_value;
            else if (comparator == Comparator::greater)
                return current_number_value > *expression_number_value;
            else if (comparator == Comparator::greaterThanEqual)
                return current_number_value >= *expression_number_value;
            return false;

        } else {
            qDebug() << "StyleExpression: expression declared as number, but it isn't -> should not happen";
            assert(false);

            // value types do not match
            return false;
        }
    } else if (m_type == DataType::Bool) {
        if (const auto expression_bool_value = std::get_if<bool>(&expression_value)) {
            const auto current_bool_value = std::get<bool>(current_value);

            if (comparator == Comparator::Equals) {
                return *expression_bool_value == current_bool_value;
            } else if (comparator == Comparator::NotEquals) {
                return *expression_bool_value != current_bool_value;
            }
        } else {
            qDebug() << "StyleExpression: expression declared as bool, but it isn't -> should not happen";
            assert(false);

            // value types do not match
            return false;
        }
    }

    qDebug() << "Unknown type in StyleExpression";
    assert(false); // shouldnt happen
    return false;
}

bool StyleExpression::matches(const mapbox::vector_tile::GeomType& type, const mapbox::feature::properties_type& properties)
{
    // qDebug() << "match: " << m_key;
    mapbox::feature::value value = extract_value(type, properties);

    if (m_comparator_in) {
        // a list -> string is sufficient
        const bool find_in = !m_negate; // if the comparator is "!in" we want to negate the values

        for (const auto& v : m_values) {
            if (compare_values(v, value, Comparator::Equals)) {
                return find_in;
            }
        }
        return !find_in; // string wasnt in the list of acceptable values
    } else if (m_comparator_has) {
        // if the value is a null value -> the value was not in the list
        bool has_not_value = std::visit([](auto&& v) { return typeid(v) == typeid(mapbox::feature::null_value); }, value);
        return m_negate ? has_not_value : !has_not_value;
    } else {

        return compare_values(m_values[0], value, m_comparator);

        // return std::visit([this](auto&& v) { return compare_visitor(v, m_values[0], m_comparator); }, value);
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
