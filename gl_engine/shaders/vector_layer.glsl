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

const float tile_extent = 4096.0;

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
    data.a = data.a - highp ivec2(tile_extent / 2.0);
    data.b = data.b - highp ivec2(tile_extent / 2.0);
    data.c = data.c - highp ivec2(tile_extent / 2.0);

    packed_data.x = highp uint(data.a.x) << (32u - 13u);
    packed_data.x = packed_data.x | (highp uint(data.a.x) & sign_mask);
    packed_data.x = packed_data.x | ((highp uint(data.a.y) & bitmask_12) << (32u - 26u));
    packed_data.x = packed_data.x | ((highp uint(data.a.y) & sign_mask) >> (13u));

    packed_data.y = highp uint(data.b.x) << (32u - 13u);
    packed_data.y = packed_data.y | (highp uint(data.b.x) & sign_mask);
    packed_data.y = packed_data.y | ((highp uint(data.b.y) & bitmask_12) << (32u - 26u));
    packed_data.y = packed_data.y | ((highp uint(data.b.y) & sign_mask) >> (13u));

    packed_data.z = highp uint(data.c.x) << (32u - 13u);
    packed_data.z = packed_data.z | (highp uint(data.c.x) & sign_mask);
    packed_data.z = packed_data.z | ((highp uint(data.c.y) & bitmask_12) << (32u - 26u));
    packed_data.z = packed_data.z | ((highp uint(data.c.y) & sign_mask) >> (13u));

    packed_data.x = packed_data.x | ((highp uint(data.style_index) >> (14u - 6u)) & bitmask_6);
    packed_data.y = packed_data.y | ((highp uint(data.style_index) >> (14u - 12u)) & bitmask_6);
    packed_data.z = packed_data.z | ((highp uint(data.style_index) & bitmask_2) << 4u);

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

    const highp int negative_bits = highp int(((1u << 31) - 1u) << 12u);

    unpacked_data.a.x = highp int((data.x & bitmask_12_first) >> (32u - 13u)) | (negative_bits & -highp int((data.x & sign_mask_first) >> 31));
    unpacked_data.a.y = highp int((data.x & bitmask_12_middle) >> (32u - 26u)) | (negative_bits & -highp int((data.x & sign_mask_middle) >> (31 - 13)));
    unpacked_data.b.x = highp int((data.y & bitmask_12_first) >> (32u - 13u)) | (negative_bits & -highp int((data.y & sign_mask_first) >> 31));
    unpacked_data.b.y = highp int((data.y & bitmask_12_middle) >> (32u - 26u)) | (negative_bits & -highp int((data.y & sign_mask_middle) >> (31 - 13)));
    unpacked_data.c.x = highp int((data.z & bitmask_12_first) >> (32u - 13u)) | (negative_bits & -highp int((data.z & sign_mask_first) >> 31));
    unpacked_data.c.y = highp int((data.z & bitmask_12_middle) >> (32u - 26u)) | (negative_bits & -highp int((data.z & sign_mask_middle) >> (31 - 13)));

    unpacked_data.style_index = (data.x & bitmask_6) << (14u - 6u);
    unpacked_data.style_index = unpacked_data.style_index | (data.y & bitmask_6) << (14u - 12u);
    unpacked_data.style_index = unpacked_data.style_index | ((data.z & (bitmask_2 << 4u)) >> 4u);

    unpacked_data.style_index = (data.x & ((1u << 6) - 1u)) << (14u - 6u);
    unpacked_data.style_index = unpacked_data.style_index | (data.y & (1u << 6) - 1u) << (14u - 12u);
    unpacked_data.style_index = unpacked_data.style_index | ((data.z & (((1u << 2) - 1u) << 4u)) >> 4u);

    // move the values back to the correct coordinates
    unpacked_data.a = unpacked_data.a + highp ivec2(tile_extent / 2.0);
    unpacked_data.b = unpacked_data.b + highp ivec2(tile_extent / 2.0);
    unpacked_data.c = unpacked_data.c + highp ivec2(tile_extent / 2.0);

    return unpacked_data;
}

