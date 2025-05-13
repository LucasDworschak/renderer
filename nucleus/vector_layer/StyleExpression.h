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
    virtual bool matches(const std::unordered_map<std::string, int>& value_map) = 0;

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
    bool matches(const std::unordered_map<std::string, int>& value_map) override;

    /**
     * @brief get_values
     * @param feature
     * @return map that converts key to an appropriate integer value. this integer value is universal and can be used to check if two values are the same
     */
    static std::unordered_map<std::string, int> get_values(const mapbox::vector_tile::feature& feature);

    static void initialize();

    /**
     * checks if the supplied argument is in a list of valid arguments
     * -> if it is this class guarantees that it is suitable to store and parse it.
     */
    static bool valid(QString value);

private:
    std::string m_key;
    Comparator m_comparator;
    std::vector<int> m_values;

    bool m_negate;
    bool m_comparator_in;
    bool m_comparator_has;

    // master map that converts strings into integer values (this way we can compare strings more efficiently)
    inline static std::unordered_map<std::string, int> string_value_map;
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
    bool matches(const std::unordered_map<std::string, int>& value_map) override;

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
