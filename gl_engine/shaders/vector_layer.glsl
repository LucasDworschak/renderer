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

#line 10020

struct VectorLayerData{
    highp ivec2 a;
    highp ivec2 b;
    highp ivec2 c;

    highp uint style_index;
    bool should_blend;
    bool is_polygon;
};

/////////////////////////////////////////////
// constants for data packing/unpacking

// defined in c++ code
// const lowp int all_bits = 32; // per output channel
// const lowp int coordinate_bits = 8;
// const highp int aa_border = 2;

const lowp int available_style_bits = all_bits - (2 * coordinate_bits_polygons);

const highp uint coordinate_bitmask = (1u << coordinate_bits_polygons) - 1u;
const highp uint coordinate_bitmask_lines = (1u << coordinate_bits_lines) - 1u;
const lowp int remaining_coordinate_bits_lines = coordinate_bits_lines - coordinate_bits_polygons;
const highp uint remaining_coordinate_bitmask_lines = (1u << remaining_coordinate_bits_lines) - 1u;

const lowp int coordinate_shift1 = all_bits - coordinate_bits_polygons;
const lowp int coordinate_shift2 = all_bits - (2 * coordinate_bits_polygons);
const lowp int coordinate_shift3 = all_bits - (3 * coordinate_bits_polygons);
const lowp int coordinate_shift4 = all_bits - (4 * coordinate_bits_polygons);

const highp uint coordinate_bitmask_shift1 = coordinate_bitmask << coordinate_shift1;
const highp uint coordinate_bitmask_shift2 = coordinate_bitmask << coordinate_shift2;
const highp uint coordinate_bitmask_shift3 = coordinate_bitmask << coordinate_shift3;
const highp uint coordinate_bitmask_shift4 = coordinate_bitmask << coordinate_shift4;
const highp uint remaining_coordinates_bitmask_shift = remaining_coordinate_bitmask_lines << remaining_coordinate_bits_lines;

// for packed.y -> we need to divide coordinate_bits_polygons by 2
const lowp int coordinate_shift1_lines = all_bits - int(0.5 * coordinate_bits_polygons);
const lowp int coordinate_shift2_lines = all_bits - int(1.0 * coordinate_bits_polygons);
const lowp int coordinate_shift3_lines = all_bits - int(1.5 * coordinate_bits_polygons);
const lowp int coordinate_shift4_lines = all_bits - int(2.0 * coordinate_bits_polygons);

const highp int cell_width_polygons = int((tile_extent * scale_polygons) / grid_size);
const highp int cell_width_lines = int((tile_extent * scale_lines) / grid_size);

const highp int max_cell_width_polygons = (1 << (coordinate_bits_polygons));
const highp int geometry_offset_polygons = (max_cell_width_polygons - cell_width_polygons) / 2;
const highp int max_cell_width_line = (1 << (coordinate_bits_lines));
const highp int geometry_offset_line = (max_cell_width_line - cell_width_lines) / 2;



// end constants for data packing/unpacking
/////////////////////////////////////////////

highp uvec2 pack_vectorlayer_data(VectorLayerData data)
{
    highp uvec2 packed_data;

    if (data.is_polygon) {
        data.a += geometry_offset_polygons;
        data.b += geometry_offset_polygons;
        data.c += geometry_offset_polygons;
    } else {
        data.a += geometry_offset_line;
        data.b += geometry_offset_line;
    }

    packed_data.x = uint(data.a.x) << coordinate_shift1;
    packed_data.x = packed_data.x | ((uint(data.a.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.x = packed_data.x | ((uint(data.b.x) & coordinate_bitmask) << coordinate_shift3);
    packed_data.x = packed_data.x | ((uint(data.b.y) & coordinate_bitmask) << coordinate_shift4);

    if (data.is_polygon) {
        packed_data.y = uint(data.c.x) << coordinate_shift1;
        packed_data.y = packed_data.y | ((uint(data.c.y) & coordinate_bitmask) << coordinate_shift2);
    } else {
        packed_data.y = ((uint(data.a.x) >> coordinate_bits_polygons) << coordinate_shift1_lines);
        packed_data.y = packed_data.y | ((uint(data.a.y) >> coordinate_bits_polygons) << coordinate_shift2_lines);
        packed_data.y = packed_data.y | ((uint(data.b.x) >> coordinate_bits_polygons) << coordinate_shift3_lines);
        packed_data.y = packed_data.y | ((uint(data.b.y) >> coordinate_bits_polygons) << coordinate_shift4_lines);
    }

    // packed_data.y = packed_data.y | ((data.style_index << 1) | ((data.is_polygon) ? 1u : 0u));
    // alternative only for neceesary for shader testing
    packed_data.y = packed_data.y | ((data.style_index << 2) | (((data.should_blend) ? 1u : 0u) << 1) | ((data.is_polygon) ? 1u : 0u));

    return packed_data;
}

VectorLayerData unpack_data(highp uvec2 packed_data)
{
    VectorLayerData unpacked_data;

    unpacked_data.a.x = int((packed_data.x & (coordinate_bitmask_shift1)) >> coordinate_shift1);
    unpacked_data.a.y = int((packed_data.x & (coordinate_bitmask_shift2)) >> coordinate_shift2);
    unpacked_data.b.x = int((packed_data.x & (coordinate_bitmask_shift3)) >> coordinate_shift3);
    unpacked_data.b.y = int((packed_data.x & (coordinate_bitmask_shift4)) >> coordinate_shift4);

    highp uvec2 c;
    c.x = (packed_data.y & (coordinate_bitmask_shift1)) >> coordinate_shift1;
    c.y = (packed_data.y & (coordinate_bitmask_shift2)) >> coordinate_shift2;

    highp uint style_and_blend = (packed_data.y & ((1u << available_style_bits) - 1u)) >> 1;
    // NOTE: the should_blend bit is only necessary for the shader -> on cpu side we do not need to separate them
    unpacked_data.style_index = style_and_blend >> 1;
    unpacked_data.should_blend = (style_and_blend & 1u) == 1u;

    unpacked_data.is_polygon = (packed_data.y & 1u) == 1u;

    if (unpacked_data.is_polygon) {
        unpacked_data.a -= geometry_offset_polygons;
        unpacked_data.b -= geometry_offset_polygons;
        unpacked_data.c = ivec2(c) - geometry_offset_polygons;
    } else {
        // unpack most significant coordinates of the line and add them to the unpacked lines
        unpacked_data.a.x = unpacked_data.a.x | int(((c.x & (remaining_coordinates_bitmask_shift)) << remaining_coordinate_bits_lines));
        unpacked_data.a.y = unpacked_data.a.y | int(((c.x & remaining_coordinate_bitmask_lines) << coordinate_bits_polygons));
        unpacked_data.b.x = unpacked_data.b.x | int(((c.y & (remaining_coordinates_bitmask_shift)) << remaining_coordinate_bits_lines));
        unpacked_data.b.y = unpacked_data.b.y | int(((c.y & remaining_coordinate_bitmask_lines) << coordinate_bits_polygons));

        unpacked_data.a -= geometry_offset_line;
        unpacked_data.b -= geometry_offset_line;
    }

    return unpacked_data;
}

