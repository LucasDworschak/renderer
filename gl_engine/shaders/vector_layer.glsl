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
    bool should_blend;
    bool is_polygon;
};

/////////////////////////////////////////////
// constants for data packing/unpacking
const lowp int  all_bits = 32; // per output channel
const lowp int  coordinate_bits_triangle = 8;
const lowp int  coordinate_bits_line = 10;
const lowp int  coordinate_bits_cell = 6; // those bits can be used to address coordinates outside cell

const lowp int  available_style_bits_triangle = all_bits - (2 * coordinate_bits_triangle);
const lowp int  available_style_bits_line = all_bits - (2 * coordinate_bits_triangle);

const lowp int  coordinate_shift1_triangle = all_bits - coordinate_bits_triangle;
const lowp int  coordinate_shift2_triangle = all_bits - (2 * coordinate_bits_triangle);
const lowp int  coordinate_shift3_triangle = all_bits - (3 * coordinate_bits_triangle);
const lowp int  coordinate_shift4_triangle = all_bits - (4 * coordinate_bits_triangle);

const lowp int  coordinate_shift1_line = all_bits - coordinate_bits_line;
const lowp int  coordinate_shift2_line = all_bits - (2 * coordinate_bits_line);
const lowp int  coordinate_shift3_line = all_bits - (3 * coordinate_bits_line);

const highp uint coordinate_bitmask_triangle = (1u << coordinate_bits_triangle) - 1u;
const highp uint coordinate_bitmask_line = (1u << coordinate_bits_line) - 1u;

const highp int cell_width = (1 << (coordinate_bits_cell));

const highp int max_cell_width_triangle = (1 << (coordinate_bits_triangle));
const highp int geometry_offset_triangle = (max_cell_width_triangle - cell_width) / 2;
const highp int max_cell_width_line = (1 << (coordinate_bits_line));
const highp int geometry_offset_line = (max_cell_width_line - cell_width) / 2;
// end constants for data packing/unpacking
/////////////////////////////////////////////

highp uvec2 pack_vectorlayer_data(VectorLayerData data)
{
    highp uvec2 packed_data;

    if(data.is_polygon)
    {
        data.a += geometry_offset_triangle;
        data.b += geometry_offset_triangle;
        data.c += geometry_offset_triangle;

        packed_data.x = uint(data.a.x) << coordinate_shift1_triangle;
        packed_data.x = packed_data.x | ((uint(data.a.y) & coordinate_bitmask_triangle) << coordinate_shift2_triangle);

        packed_data.x = packed_data.x | ((uint(data.b.x) & coordinate_bitmask_triangle) << coordinate_shift3_triangle);
        packed_data.x = packed_data.x | ((uint(data.b.y) & coordinate_bitmask_triangle) << coordinate_shift4_triangle);

        packed_data.y = uint(data.c.x) << coordinate_shift1_triangle;
        packed_data.y = packed_data.y | ((uint(data.c.y) & coordinate_bitmask_triangle) << coordinate_shift2_triangle);

        // packed_data.y = packed_data.y | ((data.style_index << 1) | ((data.is_polygon) ? 1u : 0u));
        // alternative only for neceesary for shader testing
        packed_data.y = packed_data.y | ((data.style_index << 2) | (((data.should_blend) ? 1u : 0u) << 1) | ((data.is_polygon) ? 1u : 0u));
    }
    else
    {
        data.a += geometry_offset_line;
        data.b += geometry_offset_line;

        packed_data.x = uint(data.a.x) << coordinate_shift1_line;
        packed_data.x = packed_data.x | ((uint(data.a.y) & coordinate_bitmask_line) << coordinate_shift2_line);

        packed_data.x = packed_data.x | ((uint(data.b.x) & coordinate_bitmask_line) << coordinate_shift3_line);
        packed_data.y = uint(data.b.y) << coordinate_shift1_line;

        // packed_data.y = packed_data.y | (data.style_index << 1); // | ((is_polygon) ? 1u : 0u));
        // alternative only for neceesary for shader testing
        packed_data.y = packed_data.y | ((data.style_index << 2) | (((data.should_blend) ? 1u : 0u) << 1) | ((data.is_polygon) ? 1u : 0u));

        return packed_data;
    }

    return packed_data;
}

VectorLayerData unpack_line_data(highp uvec2 packed_data)
{
    VectorLayerData unpacked_data;

    unpacked_data.c = ivec2(0,0);

    unpacked_data.a.x = int((packed_data.x & (coordinate_bitmask_line << coordinate_shift1_line)) >> coordinate_shift1_line);
    unpacked_data.a.y = int((packed_data.x & (coordinate_bitmask_line << coordinate_shift2_line)) >> coordinate_shift2_line);
    unpacked_data.b.x = int((packed_data.x & (coordinate_bitmask_line << coordinate_shift3_line)) >> coordinate_shift3_line);
    unpacked_data.b.y = int((packed_data.y & (coordinate_bitmask_line << coordinate_shift1_line)) >> coordinate_shift1_line);

    highp uint style_and_blend = (packed_data.y & ((1u << available_style_bits_triangle) - 1u)) >> 1;
    // NOTE: the should_blend bit is only necessary for the shader -> on cpu side we do not need to separate them
    unpacked_data.style_index = style_and_blend >> 1;
    unpacked_data.should_blend = (style_and_blend & 1u) == 1u;

    unpacked_data.is_polygon = (packed_data.y & 1u) == 1u;

    unpacked_data.a -= geometry_offset_line;
    unpacked_data.b -= geometry_offset_line;

    return unpacked_data;
}

VectorLayerData unpack_triangle_data(highp uvec2 packed_data)
{
    VectorLayerData unpacked_data;

    unpacked_data.a.x = int((packed_data.x & (coordinate_bitmask_triangle << coordinate_shift1_triangle)) >> coordinate_shift1_triangle);
    unpacked_data.a.y = int((packed_data.x & (coordinate_bitmask_triangle << coordinate_shift2_triangle)) >> coordinate_shift2_triangle);
    unpacked_data.b.x = int((packed_data.x & (coordinate_bitmask_triangle << coordinate_shift3_triangle)) >> coordinate_shift3_triangle);
    unpacked_data.b.y = int((packed_data.x & (coordinate_bitmask_triangle << coordinate_shift4_triangle)) >> coordinate_shift4_triangle);
    unpacked_data.c.x = int((packed_data.y & (coordinate_bitmask_triangle << coordinate_shift1_triangle)) >> coordinate_shift1_triangle);
    unpacked_data.c.y = int((packed_data.y & (coordinate_bitmask_triangle << coordinate_shift2_triangle)) >> coordinate_shift2_triangle);

    highp uint style_and_blend = (packed_data.y & ((1u << available_style_bits_triangle) - 1u)) >> 1;
    // NOTE: the should_blend bit is only necessary for the shader -> on cpu side we do not need to separate them
    unpacked_data.style_index = style_and_blend >> 1;
    unpacked_data.should_blend = (style_and_blend & 1u) == 1u;

    unpacked_data.is_polygon = (packed_data.y & 1u) == 1u;

    unpacked_data.a -= geometry_offset_triangle;
    unpacked_data.b -= geometry_offset_triangle;
    unpacked_data.c -= geometry_offset_triangle;

    return unpacked_data;
}

VectorLayerData unpack_data(highp uvec2 packed_data)
{

    if ((packed_data.y & 1u) == 1u) {
        return unpack_triangle_data(packed_data);
    } else {
        return unpack_line_data(packed_data);
    }
}

