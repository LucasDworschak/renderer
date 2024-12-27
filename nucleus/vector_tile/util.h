/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Jörg Christian Reiher
 * Copyright (C) 2024 Adam Celarek
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
#include <QString>
#include <mapbox/vector_tile.hpp>

namespace nucleus::vector_tile::util {
namespace detail {
    class StringPrintVisitor {
    public:
        QString operator()(std::vector<mapbox::feature::value>) const { return QStringLiteral("vector"); }
        QString operator()(std::unordered_map<std::string, mapbox::feature::value>) const { return QStringLiteral("unordered_map"); }
        QString operator()(mapbox::feature::null_value_t) const { return QStringLiteral("null"); }
        QString operator()(std::nullptr_t) const { return QStringLiteral("nullptr"); }
        QString operator()(uint64_t val) const { return QString::number(val); }
        QString operator()(int64_t val) const { return QString::number(val); }
        QString operator()(double val) const { return QString::number(val); }
        QString operator()(std::string const& val) const { return QString::fromStdString(val); }
        QString operator()(bool val) const { return val ? QStringLiteral("true") : QStringLiteral("false"); }
    };

    class CompareVisitor {
    public:
        bool operator()(std::vector<mapbox::feature::value>, QString, QString) const { return false; }
        bool operator()(std::unordered_map<std::string, mapbox::feature::value>, QString, QString) const { return false; }
        bool operator()(mapbox::feature::null_value_t, QString, QString) const { return false; }
        bool operator()(std::nullptr_t, QString, QString) const { return false; }
        bool operator()(uint64_t val, QString compare, QString comparator) const
        {
            // qDebug() << "compare uint: " << val << comparator << compare.toULongLong();
            if (comparator == "==")
                return val == compare.toULongLong();
            else if (comparator == "!=")
                return val != compare.toULongLong();
            else if (comparator == "<")
                return val < compare.toULongLong();
            else if (comparator == "<=")
                return val <= compare.toULongLong();
            else if (comparator == ">")
                return val > compare.toULongLong();
            else if (comparator == ">=")
                return val >= compare.toULongLong();
            return false;
        }
        bool operator()(int64_t val, QString compare, QString comparator) const
        {
            // qDebug() << "compare int: " << val << comparator << compare.toLongLong();
            if (comparator == "==")
                return val == compare.toLongLong();
            else if (comparator == "!=")
                return val != compare.toLongLong();
            else if (comparator == "<")
                return val < compare.toLongLong();
            else if (comparator == "<=")
                return val <= compare.toLongLong();
            else if (comparator == ">")
                return val > compare.toLongLong();
            else if (comparator == ">=")
                return val >= compare.toLongLong();
            return false;
        }
        bool operator()(double val, QString compare, QString comparator) const
        {
            // qDebug() << "compare double: " << val << comparator << compare.toDouble();
            if (comparator == "==")
                return val == compare.toDouble();
            else if (comparator == "!=")
                return val != compare.toDouble();
            else if (comparator == "<")
                return val < compare.toDouble();
            else if (comparator == "<=")
                return val <= compare.toDouble();
            else if (comparator == ">")
                return val > compare.toDouble();
            else if (comparator == ">=")
                return val >= compare.toDouble();
            return false;
        }
        bool operator()(std::string const& val, QString compare, QString comparator) const
        {
            // qDebug() << "compare string: " << val << comparator << compare.toStdString();
            if (comparator == "==")
                return val == compare.toStdString();
            else if (comparator == "!=")
                return val != compare.toStdString();
            return false;
        }
        bool operator()(bool val, QString compare, QString comparator) const
        {
            // qDebug() << "compare bool: " << val << comparator << compare.toStdString();
            if (comparator == "==")
                return val == (compare.toLower() == "true") ? true : false;
            else if (comparator == "!=")
                return val != (compare.toLower() == "true") ? true : false;

            return false;
        }
    };

} // namespace detail
const static detail::StringPrintVisitor string_print_visitor = {};
const static detail::CompareVisitor compare_visitor = {};

} // namespace nucleus::vector_tile::util
