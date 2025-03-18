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
    virtual bool matches(const mapbox::vector_tile::GeomType& type, const mapbox::feature::properties_type& properties) = 0;

    static std::unique_ptr<StyleExpressionBase> create_filter_expression(QJsonArray data);
    QJsonValue extract_literal(QJsonValue expression);
};

class StyleExpression : public StyleExpressionBase {
public:
    StyleExpression(QJsonArray data);

    size_t hash() override;

    /**
     * checks if the supplied feature matches the expression that was initialized
     */
    bool matches(const mapbox::vector_tile::GeomType& type, const mapbox::feature::properties_type& properties) override;

    /**
     * checks if the supplied argument is in a list of valid arguments
     * -> if it is this class guarantees that it is suitable to store and parse it.
     */
    static bool valid(QString value);

private:
    std::string m_key;
    std::string m_comparator;
    std::vector<std::string> m_values;

    mapbox::feature::value extract_value(const mapbox::vector_tile::GeomType& type, const mapbox::feature::properties_type& properties);
};

class StyleExpressionCollection : public StyleExpressionBase {
public:
    StyleExpressionCollection(QJsonArray data);

    size_t hash() override;

    /**
     * checks if the supplied feature matches the expression of all/any of the expression it holds
     */
    bool matches(const mapbox::vector_tile::GeomType& type, const mapbox::feature::properties_type& properties) override;

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
