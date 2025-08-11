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
    highp vec2 a;
    highp vec2 b;
    highp vec2 c;

    highp uint style_index;
    bool is_polygon;
};

/////////////////////////////////////////////
// constants for data packing/unpacking

// defined in c++ code
// const lowp int all_bits = 32; // per output channel
// const lowp int coordinate_bits = 8;

const highp vec2 cell_size = vec2(tile_extent) / vec2(grid_size);


const highp uint coordinate_bitmask = (1u << coordinate_bits_polygons) - 1u;
const highp uint coordinate_bitmask_lines = (1u << coordinate_bits_lines) - 1u;
const lowp int remaining_coordinate_bits_lines = coordinate_bits_lines - coordinate_bits_polygons;
const highp uint remaining_coordinate_bitmask_lines = (1u << remaining_coordinate_bits_lines) - 1u;
const highp uint is_polygon_bitmask = (1u << style_bits);

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
const lowp int coordinate_shift1_lines = all_bits - int(0.5 * float(coordinate_bits_polygons));
const lowp int coordinate_shift2_lines = all_bits - int(1.0 * float(coordinate_bits_polygons));
const lowp int coordinate_shift3_lines = all_bits - int(1.5 * float(coordinate_bits_polygons));
const lowp int coordinate_shift4_lines = all_bits - int(2.0 * float(coordinate_bits_polygons));

const highp int cell_width_polygons = int((tile_extent * scale_polygons) / grid_size.x);
const highp int cell_width_lines = int((tile_extent * scale_lines) / grid_size.x);

const highp int max_cell_width_polygons = (1 << (coordinate_bits_polygons));
const highp int geometry_offset_polygons = (max_cell_width_polygons - cell_width_polygons) / 2;
const highp int max_cell_width_line = (1 << (coordinate_bits_lines));
const highp int geometry_offset_line = (max_cell_width_line - cell_width_lines) / 2;

const highp uint style_bitmask = ((1u << style_bits) - 1u);

// end constants for data packing/unpacking
/////////////////////////////////////////////

highp uvec2 pack_vectorlayer_data(VectorLayerData data)
{
    highp uvec2 packed_data;

    ivec2 a = ivec2(data.a);
    ivec2 b = ivec2(data.b);
    ivec2 c = ivec2(data.c);

    if (data.is_polygon) {
        a += geometry_offset_polygons;
        b += geometry_offset_polygons;
        c += geometry_offset_polygons;
    } else {
        a += geometry_offset_line;
        b += geometry_offset_line;
    }

    packed_data.x = uint(a.x) << coordinate_shift1;
    packed_data.x = packed_data.x | ((uint(a.y) & coordinate_bitmask) << coordinate_shift2);

    packed_data.x = packed_data.x | ((uint(b.x) & coordinate_bitmask) << coordinate_shift3);
    packed_data.x = packed_data.x | ((uint(b.y) & coordinate_bitmask) << coordinate_shift4);

    if (data.is_polygon) {
        packed_data.y = uint(c.x) << coordinate_shift1;
        packed_data.y = packed_data.y | ((uint(c.y) & coordinate_bitmask) << coordinate_shift2);
    } else {
        packed_data.y = ((uint(a.x) >> coordinate_bits_polygons) << coordinate_shift1_lines);
        packed_data.y = packed_data.y | ((uint(a.y) >> coordinate_bits_polygons) << coordinate_shift2_lines);
        packed_data.y = packed_data.y | ((uint(b.x) >> coordinate_bits_polygons) << coordinate_shift3_lines);
        packed_data.y = packed_data.y | ((uint(b.y) >> coordinate_bits_polygons) << coordinate_shift4_lines);
    }

    highp uint is_polygon = (data.is_polygon ? 1u : 0u) << style_bits;

    packed_data.y = packed_data.y | is_polygon | data.style_index;

    return packed_data;
}

highp uint unpack_style_index(highp uvec2 packed_data)
{
    return (packed_data.y & style_bitmask) * uint(max_zoom+1);
}

bool is_polygon(highp uvec2 packed_data)
{
    return (packed_data.y & is_polygon_bitmask) != 0u;
}

VectorLayerData unpack_data(highp uvec2 packed_data, lowp vec2 grid_cell)
{
    VectorLayerData unpacked_data;

    ivec2 a;
    ivec2 b;
    ivec2 c = ivec2(0,0);

    a.x = int((packed_data.x & (coordinate_bitmask_shift1)) >> coordinate_shift1);
    a.y = int((packed_data.x & (coordinate_bitmask_shift2)) >> coordinate_shift2);
    b.x = int((packed_data.x & (coordinate_bitmask_shift3)) >> coordinate_shift3);
    b.y = int((packed_data.x & (coordinate_bitmask_shift4)) >> coordinate_shift4);

    highp uvec2 c_u;
    c_u.x = (packed_data.y & (coordinate_bitmask_shift1)) >> coordinate_shift1;
    c_u.y = (packed_data.y & (coordinate_bitmask_shift2)) >> coordinate_shift2;

    unpacked_data.style_index = unpack_style_index(packed_data);

    unpacked_data.is_polygon = is_polygon(packed_data);

    if (unpacked_data.is_polygon) {
        a -= geometry_offset_polygons;
        b -= geometry_offset_polygons;
        c = ivec2(c_u) - geometry_offset_polygons;
    } else {
        // unpack most significant coordinates of the line and add them to the unpacked lines
        a.x = a.x | int(((c_u.x & (remaining_coordinates_bitmask_shift)) << remaining_coordinate_bits_lines));
        a.y = a.y | int(((c_u.x & remaining_coordinate_bitmask_lines) << coordinate_bits_polygons));
        b.x = b.x | int(((c_u.y & (remaining_coordinates_bitmask_shift)) << remaining_coordinate_bits_lines));
        b.y = b.y | int(((c_u.y & remaining_coordinate_bitmask_lines) << coordinate_bits_polygons));

        a -= geometry_offset_line;
        b -= geometry_offset_line;
    }

    highp float tile_scale = float(scale_lines) * (1.0-float(unpacked_data.is_polygon)) + float(scale_polygons) * float(unpacked_data.is_polygon);

    highp vec2 cell_offset = grid_cell * cell_size * tile_scale;
    highp float division_extent = 1.0 / (tile_extent * tile_scale);

    unpacked_data.a = (vec2(a) + cell_offset) * division_extent;
    unpacked_data.b = (vec2(b) + cell_offset) * division_extent;
    unpacked_data.c = (vec2(c) + cell_offset) * division_extent;

    return unpacked_data;
}


// for the base unittest we are not interested in testing if the polygons scaling and uv scaling works
// we only are interested if the data we send to the gpu is correctly encoded and decoded without anly loss of precision
VectorLayerData normalize_unpack_for_unittest(VectorLayerData unpacked_data, lowp vec2 grid_cell)
{
    highp float tile_scale = float(scale_lines) * (1.0-float(unpacked_data.is_polygon)) + float(scale_polygons) * float(unpacked_data.is_polygon);

    highp vec2 cell_offset = grid_cell * cell_size * tile_scale;
    highp float division_extent = 1.0 / (tile_extent * tile_scale);

    unpacked_data.a = (vec2(unpacked_data.a) / division_extent) - cell_offset;
    unpacked_data.b = (vec2(unpacked_data.b) / division_extent) - cell_offset;
    unpacked_data.c = (vec2(unpacked_data.c) / division_extent) - cell_offset;

    unpacked_data.style_index /= uint(max_zoom+1);

    return unpacked_data;
}



