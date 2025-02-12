/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2025 Lucas Dworschak
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

const highp float tile_extent = 512.0;

struct VectorLayerData{
    highp ivec2 a;
    highp ivec2 b;
    highp ivec2 c;
    highp uint style_index;
};

highp uvec3 pack_vectorlayer_data(VectorLayerData data) {
    highp uvec3 packed_data;

    const highp uint sign_mask = (1u << 31);
    const highp uint bitmask_12 = (1u << 12) - 1u;
    const highp uint bitmask_6 = (1u << 6) - 1u;
    const highp uint bitmask_2 = (1u << 2) - 1u;

    // move the values by half extent
    const highp ivec2 half_extent = ivec2(tile_extent / 2.0);
    highp uvec2 a = uvec2(data.a - half_extent);
    highp uvec2 b = uvec2(data.b - half_extent);
    highp uvec2 c = uvec2(data.c - half_extent);

    packed_data.x = a.x << (32u - 13u);
    packed_data.x = packed_data.x | (a.x & sign_mask);
    packed_data.x = packed_data.x | ((a.y & bitmask_12) << (32u - 26u));
    packed_data.x = packed_data.x | ((a.y & sign_mask) >> (13u));

    packed_data.y = b.x << (32u - 13u);
    packed_data.y = packed_data.y | (b.x & sign_mask);
    packed_data.y = packed_data.y | ((b.y & bitmask_12) << (32u - 26u));
    packed_data.y = packed_data.y | ((b.y & sign_mask) >> (13u));

    packed_data.z = c.x << (32u - 13u);
    packed_data.z = packed_data.z | (c.x & sign_mask);
    packed_data.z = packed_data.z | ((c.y & bitmask_12) << (32u - 26u));
    packed_data.z = packed_data.z | ((c.y & sign_mask) >> (13u));

    packed_data.x = packed_data.x | ((data.style_index >> (14u - 6u)) & bitmask_6);
    packed_data.y = packed_data.y | ((data.style_index >> (14u - 12u)) & bitmask_6);
    packed_data.z = packed_data.z | ((data.style_index & bitmask_2) << 4u);

    return packed_data;
}

VectorLayerData unpack_vectorlayer_data(highp uvec3 data) {
    VectorLayerData unpacked_data;

    const highp uint sign_mask_first = (1u << 31);
    const highp uint sign_mask_middle = (1u << (31u - 13u));
    const highp uint bitmask_12_first = ((1u << 12) - 1u) << (32u - 13u);
    const highp uint bitmask_12_middle = ((1u << 12) - 1u) << (32u - 26u);
    const highp uint bitmask_6 = (1u << 6) - 1u;
    const highp uint bitmask_2 = (1u << 2) - 1u;

    const highp int negative_bits = int(((1u << 31) - 1u) << 12u);

    highp int x_first = int((data.x & bitmask_12_first) >> (32u - 13u));
    highp int x_middle = int((data.x & bitmask_12_middle) >> (32u - 26u));
    highp int x_sign_first = int((data.x & sign_mask_first) >> 31);
    highp int x_sign_middle = int((data.x & sign_mask_middle) >> (31 - 13));

    highp int y_first = int((data.y & bitmask_12_first) >> (32u - 13u));
    highp int y_middle = int((data.y & bitmask_12_middle) >> (32u - 26u));
    highp int y_sign_first = int((data.y & sign_mask_first) >> 31);
    highp int y_sign_middle = int((data.y & sign_mask_middle) >> (31 - 13));

    highp int z_first = int((data.z & bitmask_12_first) >> (32u - 13u));
    highp int z_middle = int((data.z & bitmask_12_middle) >> (32u - 26u));
    highp int z_sign_first = int((data.z & sign_mask_first) >> 31);
    highp int z_sign_middle = int((data.z & sign_mask_middle) >> (31 - 13));

    unpacked_data.a.x = x_first | (negative_bits & -x_sign_first);
    unpacked_data.a.y = x_middle | (negative_bits & -x_sign_middle);
    unpacked_data.b.x = y_first | (negative_bits & -y_sign_first);
    unpacked_data.b.y = y_middle | (negative_bits & -y_sign_middle);
    unpacked_data.c.x = z_first | (negative_bits & -z_sign_first);
    unpacked_data.c.y = z_middle | (negative_bits & -z_sign_middle);

    unpacked_data.style_index = (data.x & bitmask_6) << (14u - 6u);
    unpacked_data.style_index = unpacked_data.style_index | (data.y & bitmask_6) << (14u - 12u);
    unpacked_data.style_index = unpacked_data.style_index | ((data.z & (bitmask_2 << 4u)) >> 4u);

    unpacked_data.style_index = (data.x & ((1u << 6) - 1u)) << (14u - 6u);
    unpacked_data.style_index = unpacked_data.style_index | (data.y & (1u << 6) - 1u) << (14u - 12u);
    unpacked_data.style_index = unpacked_data.style_index | ((data.z & (((1u << 2) - 1u) << 4u)) >> 4u);

    // move the values back to the correct coordinates
    const highp ivec2 half_extent = ivec2(tile_extent / 2.0);
    unpacked_data.a = unpacked_data.a + half_extent;
    unpacked_data.b = unpacked_data.b + half_extent;
    unpacked_data.c = unpacked_data.c + half_extent;

    return unpacked_data;
}

