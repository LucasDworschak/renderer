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
    highp uint style_index;
};

/////////////////////////////////////////////
// constants for data packing/unpacking
const lowp int all_bits = 32;
const lowp int coordinate_bits = 14;

const lowp int style_bits_per_coord = all_bits - (2 * coordinate_bits);
const lowp int all_style_bits = (style_bits_per_coord * 3);

const lowp int coordinate_shift1 = all_bits - coordinate_bits;
const lowp int coordinate_shift2 = all_bits - (2 * coordinate_bits);
const lowp int style_shift1 = all_style_bits - style_bits_per_coord;
const lowp int style_shift2 = all_style_bits - (2 * style_bits_per_coord);

const highp uint coordinate_bitmask = (1u << coordinate_bits) - 1u;
const highp uint style_bitmask = (1u << style_bits_per_coord) - 1u;
// end constants for data packing/unpacking
/////////////////////////////////////////////

highp uvec3 pack_vectorlayer_data(VectorLayerData data) {
    highp uvec3 packed_data;

    // the extent from the tile gives us e.g. 4096
    // in reality the coordinates can be a little over and a little under the extent -> something like [-128, 4096+128]
    // solution: coordinate normalization (move from [-tile_extent - tile_extent] to [0 - 2*tile_extent])
    data.a += ivec2(tile_extent);
    data.b += ivec2(tile_extent);
    data.c += ivec2(tile_extent);

    packed_data.x = uint(data.a.x) << coordinate_shift1;
    packed_data.x = packed_data.x | ((uint(data.a.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.y = uint(data.b.x) << coordinate_shift1;
    packed_data.y = packed_data.y | ((uint(data.b.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.z = uint(data.c.x) << coordinate_shift1;
    packed_data.z = packed_data.z | ((uint(data.c.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.x = packed_data.x | ((data.style_index >> style_shift1) & style_bitmask);
    packed_data.y = packed_data.y | ((data.style_index >> style_shift2) & style_bitmask);
    packed_data.z = packed_data.z | ((data.style_index & style_bitmask));

    return packed_data;
}

VectorLayerData unpack_vectorlayer_data(highp uvec3 packed_data) {
    VectorLayerData unpacked_data;

    unpacked_data.a.x = int((packed_data.x & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    unpacked_data.a.y = int((packed_data.x & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);
    unpacked_data.b.x = int((packed_data.y & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    unpacked_data.b.y = int((packed_data.y & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);
    unpacked_data.c.x = int((packed_data.z & (coordinate_bitmask << coordinate_shift1)) >> coordinate_shift1);
    unpacked_data.c.y = int((packed_data.z & (coordinate_bitmask << coordinate_shift2)) >> coordinate_shift2);

    unpacked_data.style_index = (packed_data.x & style_bitmask) << style_shift1;
    unpacked_data.style_index = unpacked_data.style_index | (packed_data.y & style_bitmask) << style_shift2;
    unpacked_data.style_index = unpacked_data.style_index | (packed_data.z & style_bitmask);

    // move the values back to the correct coordinates
    unpacked_data.a -= ivec2(tile_extent);
    unpacked_data.b -= ivec2(tile_extent);
    unpacked_data.c -= ivec2(tile_extent);

    return unpacked_data;
}

