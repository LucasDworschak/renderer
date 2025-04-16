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

const highp float tile_extent = 4096.0;

struct VectorLayerData{
    highp ivec2 a;
    highp ivec2 b;
    highp ivec2 c;

    highp uint style_layer;
};

/////////////////////////////////////////////
// constants for data packing/unpacking
const lowp int all_bits = 32;
const lowp int coordinate_bits = 8;
const lowp int coordinate_bits_buffer = 2; // those bits can be used to address coordinates outside cell

const lowp int style_bits = all_bits - (2 * coordinate_bits);

const lowp int coordinate_shift1 = all_bits - coordinate_bits;
const lowp int coordinate_shift2 = all_bits - (2 * coordinate_bits);
const lowp int coordinate_shift3 = all_bits - (3 * coordinate_bits);
const lowp int coordinate_shift4 = all_bits - (4 * coordinate_bits);

const highp uint coordinate_bitmask = (1u << coordinate_bits) - 1u;
const highp int geometry_offset = int(coordinate_bitmask >> coordinate_bits_buffer);
// end constants for data packing/unpacking
/////////////////////////////////////////////

highp uvec2 pack_vectorlayer_data(VectorLayerData data) {

    highp uvec2 packed_data;

    data.a += geometry_offset;
    data.b += geometry_offset;
    data.c += geometry_offset;

    packed_data.x = uint(data.a.x) << coordinate_shift1;
    packed_data.x = packed_data.x | ((uint(data.a.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.x = packed_data.x | ((uint(data.b.x) & coordinate_bitmask) << coordinate_shift3);
    packed_data.x = packed_data.x | ((uint(data.b.y) & coordinate_bitmask) << coordinate_shift4);

    packed_data.y = uint(data.c.x) << coordinate_shift1;
    packed_data.y = packed_data.y | ((uint(data.c.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.y = packed_data.y | data.style_layer;

    return packed_data;
}

VectorLayerData unpack_vectorlayer_data(highp uvec2 packed_data) {

    VectorLayerData unpacked_data;

    unpacked_data.a.x = int((packed_data.x & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    unpacked_data.a.y = int((packed_data.x & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);
    unpacked_data.b.x = int((packed_data.x & (coordinate_bitmask << coordinate_shift3)) >> coordinate_shift3);
    unpacked_data.b.y = int((packed_data.x & (coordinate_bitmask << coordinate_shift4)) >> coordinate_shift4);
    unpacked_data.c.x = int((packed_data.y & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    unpacked_data.c.y = int((packed_data.y & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);

    unpacked_data.style_layer = packed_data.y & ((1u << style_bits) - 1u);

    unpacked_data.a -= geometry_offset;
    unpacked_data.b -= geometry_offset;
    unpacked_data.c -= geometry_offset;

    return unpacked_data;
}

