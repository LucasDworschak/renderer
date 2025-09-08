/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Lucas Dworschak
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

#include <QSignalSpy>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <QFile>
#include <QImage>
#include <QString>
#include <glm/glm.hpp>

#include <nucleus/tile/conversion.h>
#include <nucleus/tile/utils.h>
#include <radix/tile.h>

#include "nucleus/Raster.h"

#include "nucleus/utils/rasterizer.h"
#include "nucleus/vector_layer/GpuMultiArrayHelper.h"
#include "nucleus/vector_layer/constants.h"

using namespace nucleus::vector_layer;

const glm::uvec2 id_to_pixel(nucleus::tile::Id id, nucleus::Raster<glm::u32vec2> packed_ids)
{
    const auto hash_to_pixel = [](uint16_t hash) { return glm::uvec2(hash & 255, hash >> 8); };

    auto hash = nucleus::srs::hash_uint16(id);
    auto wanted_packed_tile_id = nucleus::srs::pack(id);

    auto found_packed_tile_id = packed_ids.pixel(hash_to_pixel(hash));
    auto iter = 0u;

    while (found_packed_tile_id != wanted_packed_tile_id && found_packed_tile_id != glm::u32vec2(-1, -1)) {
        hash++;
        found_packed_tile_id = packed_ids.pixel(hash_to_pixel(hash));
        if (iter++ > 50u) {
            break;
        }
    }

    return hash_to_pixel(hash);
};

TEST_CASE("nucleus/gpu_array_helper")
{
    SECTION("Multi Array Helper - buffer sizes")
    {
        CHECK(custom_array_layer_index() == 2); // make sure that the index where custom data starts is still the same (if not change the checks below)

        {
            GpuMultiArrayHelper helper;

            CHECK(helper.buffer_amount() == 2);

            helper.set_tile_limit(2048);

            CHECK(helper.layer_amount(0) == 2048);
            CHECK(helper.layer_amount(1) == 2048);
            CHECK(helper.layer_amount(2) == constants::array_layer_tile_amount[2]);
        }
        {
            GpuMultiArrayHelper helper;

            helper.set_tile_limit(128);

            CHECK(helper.layer_amount(0) == 128);
            CHECK(helper.layer_amount(1) == 128);
            CHECK(helper.layer_amount(2) == 128); // although it is a custom size, the lower value is used instead
        }

    }

    SECTION("Multi Array Helper - tile encoding in same buffer")
    {
        constexpr auto buffer_offset = (constants::array_helper_all_bits - constants::array_helper_buffer_info_bits);
        constexpr auto layer_mask = ((1u << buffer_offset) - 1u);

        GpuMultiArrayHelper helper;

        helper.set_tile_limit(2048);

        const nucleus::tile::Id id0 = { 10, { 10, 10 } };
        const nucleus::tile::Id id1 = { 12, { 12, 12 } };
        const nucleus::tile::Id id2 = { 13, { 13, 13 } };

        auto layer0 = helper.add_tile(id0, 0);
        auto layer1 = helper.add_tile(id1, 0);
        auto layer2 = helper.add_tile(id2, 0);

        auto [packed_ids, layers] = helper.generate_dictionary();

        auto p_layer0 = layers.pixel(id_to_pixel(id0, packed_ids));
        auto p_layer1 = layers.pixel(id_to_pixel(id1, packed_ids));
        auto p_layer2 = layers.pixel(id_to_pixel(id2, packed_ids));

        auto buffer_info0 = (p_layer0.y & ((-1u << buffer_offset))) >> buffer_offset;
        auto buffer_info1 = (p_layer1.y & ((-1u << buffer_offset))) >> buffer_offset;
        auto buffer_info2 = (p_layer2.y & ((-1u << buffer_offset))) >> buffer_offset;

        p_layer0.x = p_layer0.x & layer_mask;
        p_layer1.x = p_layer1.x & layer_mask;
        p_layer2.x = p_layer2.x & layer_mask;
        p_layer0.y = p_layer0.y & layer_mask;
        p_layer1.y = p_layer1.y & layer_mask;
        p_layer2.y = p_layer2.y & layer_mask;

        CHECK(buffer_info0 == 0);
        CHECK(buffer_info1 == 0);
        CHECK(buffer_info2 == 0);

        CHECK(layer0 == p_layer0);
        CHECK(layer1 == p_layer1);
        CHECK(layer2 == p_layer2);
    }

    SECTION("Multi Array Helper - tile encoding in different buffer")
    {
        constexpr auto buffer_offset = (constants::array_helper_all_bits - constants::array_helper_buffer_info_bits);
        constexpr auto layer_mask = ((1u << buffer_offset) - 1u);

        GpuMultiArrayHelper helper;

        helper.set_tile_limit(2048);

        const nucleus::tile::Id id0 = { 6, { 23, 56 } };
        const nucleus::tile::Id id1 = { 8, { 122, 132 } };
        const nucleus::tile::Id id2 = { 10, { 1553, 1663 } };

        auto layer0 = helper.add_tile(id0, 0);
        auto layer1 = helper.add_tile(id1, 1);
        auto layer2 = helper.add_tile(id2, 2);

        auto [packed_ids, layers] = helper.generate_dictionary();

        auto p_layer0 = layers.pixel(id_to_pixel(id0, packed_ids));
        auto p_layer1 = layers.pixel(id_to_pixel(id1, packed_ids));
        auto p_layer2 = layers.pixel(id_to_pixel(id2, packed_ids));

        auto buffer_info0 = (p_layer0.y & ((-1u << buffer_offset))) >> buffer_offset;
        auto buffer_info1 = (p_layer1.y & ((-1u << buffer_offset))) >> buffer_offset;
        auto buffer_info2 = (p_layer2.y & ((-1u << buffer_offset))) >> buffer_offset;

        p_layer0.x = p_layer0.x & layer_mask;
        p_layer1.x = p_layer1.x & layer_mask;
        p_layer2.x = p_layer2.x & layer_mask;
        p_layer0.y = p_layer0.y & layer_mask;
        p_layer1.y = p_layer1.y & layer_mask;
        p_layer2.y = p_layer2.y & layer_mask;

        CHECK(buffer_info0 == 0);
        CHECK(buffer_info1 == 1);
        CHECK(buffer_info2 == 2);

        // make sure that the values are the expected values
        CHECK(glm::u16vec2(0, 0) == p_layer0);
        CHECK(glm::u16vec2(1, 1) == p_layer1);
        CHECK(glm::u16vec2(2, 0) == p_layer2);

        // make sure that the values from add_tile are the correct values
        CHECK(layer0 == p_layer0);
        CHECK(layer1 == p_layer1);
        CHECK(layer2 == p_layer2);
    }
}
