/*****************************************************************************
 * Alpine Renderer
 * Copyright (C) 2022 Adam Celarek
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

#include <catch2/catch_test_macros.hpp>

#include "nucleus/map_label/Charset.h"

#include <set>

TEST_CASE("nucleus/charset")
{

    // Checks if reading the charset file and creating the charset class works
    SECTION("charset reading")
    {
        nucleus::maplabel::Charset& c = nucleus::maplabel::Charset::get_instance();
        
        std::set<char16_t> all_chars = c.all_chars();

        // is the number of chars correct
        CHECK(all_chars.size() == 123);
        // does the set contain certain chars (e.g. " " -> 32)
        CHECK(all_chars.contains(32));
    }

    // check if inserting new values and detecting those inserts works
    SECTION("char insert/compare")
    {
        std::set<char16_t> all_char_list = {};

        // initialize the all_char_list with values from charset.txt
        {
            nucleus::maplabel::Charset& c = nucleus::maplabel::Charset::get_instance();
            auto cached_chars = c.all_chars();
            all_char_list.insert(cached_chars.begin(), cached_chars.end());
        }

        // insert a new char to the charset class
        // NOTE: this step only happens in VectorTileManager.cpp
        {
            nucleus::maplabel::Charset& c = nucleus::maplabel::Charset::get_instance();
            
            auto cached_chars = c.all_chars();

            // making sure that this test value was not already present in the charset.txt
            CHECK(!cached_chars.contains(28779));
            cached_chars.insert(28779); // insert this test value to the set
            // update the charset with new chars
            c.add_chars(cached_chars);
        }

        // check if new chars have been found
        // NOTE: this step only happens in LabelFactory.cpp
        {
            nucleus::maplabel::Charset& c = nucleus::maplabel::Charset::get_instance();

            CHECK(c.is_update_necessary(all_char_list.size()));
            
            auto diff = c.char_diff(all_char_list);

            CHECK(diff.size() == 1);
            CHECK(diff.contains(28779));
        }
    }
}
