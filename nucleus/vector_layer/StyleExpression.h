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

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

#include <unordered_set>

#include <nucleus/vector_layer/constants.h>

namespace mapbox { // forward declare mapbox classes
namespace vector_tile {
    enum GeomType : std::uint8_t;
    class feature;
}

namespace feature {
    class value;
    // class properties_type;
    using properties_type = std::unordered_map<std::string, value>;
}
} // namespace mapbox

namespace nucleus::vector_layer {

class StyleExpressionBase {
public:
    StyleExpressionBase() { }
    virtual ~StyleExpressionBase() { }

    virtual size_t hash() = 0;

    /**
     * checks if the supplied feature matches the expression that was initialized
     */
    virtual bool matches(const std::array<int, constants::max_style_expression_keys>& values) = 0;

    static std::unique_ptr<StyleExpressionBase> create_filter_expression(QJsonArray data);
    QJsonValue extract_literal(QJsonValue expression);
};

enum class Comparator { Equals = 1, NotEquals = 2, less = 3, lessThanEqual = 4, greater = 5, greaterThanEqual = 6 };

class StyleExpression : public StyleExpressionBase {
public:
    StyleExpression(QJsonArray data);

    size_t hash() override;

    /**
     * checks if the supplied feature matches the expression that was initialized
     */
    bool matches(const std::array<int, constants::max_style_expression_keys>& values) override;

    /**
     * parses the values of the incomming feature and maps them to an array. the index of the array corresponds to the m_key variable, while the value is used
     * to compare to the m_values
     */
    static void get_values(const mapbox::vector_tile::feature& feature, std::array<int, constants::max_style_expression_keys>* values);

    static void initialize();

    /**
     * checks if the supplied argument is in a list of valid arguments
     * -> if it is this class guarantees that it is suitable to store and parse it.
     */
    static bool valid(QString value);

private:
    int m_key;
    Comparator m_comparator;
    std::vector<int> m_values;

    bool m_negate;
    bool m_comparator_in;
    bool m_comparator_has;

    // master map that converts strings into integer values (this way we can compare strings more efficiently)
    inline static std::unordered_map<std::string, int> all_keys;
    inline static std::unordered_map<std::string, int> string_value_map;
    inline static int last_key_index; // saves the index of the last element inserted into the all_keys map
    inline static int counter; // counter that is used for the string value map
    static constexpr int null_value = INT_MAX;
    static constexpr int float_precision = 1000; // multiplier for float values to make them integers
};

class StyleExpressionCollection : public StyleExpressionBase {
public:
    StyleExpressionCollection(QJsonArray data);

    size_t hash() override;

    /**
     * checks if the supplied feature matches the expression of all/any of the expression it holds
     */
    bool matches(const std::array<int, constants::max_style_expression_keys>& values) override;

    /**
     * checks if the supplied argument is in a list of valid arguments
     * -> if it is this class guarantees that it is suitable to store and parse it.
     */
    static bool valid(QString value);

private:
    // m_all==true: must match m_all the Expressions
    // m_all==false: (=any) -> must match at least one Expression
    bool m_all;
    bool m_negate;

    std::vector<std::unique_ptr<StyleExpressionBase>> m_subFilters;
};

} // namespace nucleus::vector_layer
