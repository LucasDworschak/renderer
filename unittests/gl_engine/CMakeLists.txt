#############################################################################
# Alpine Terrain Renderer
# Copyright (C) 2023 Adam Celarek <family name at cg tuwien ac at>
# Copyright (C) 2023 Gerald Kimmersdorfer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#############################################################################

project(alpine-renderer-unittests_gl_engine LANGUAGES CXX)

alp_add_unittest(unittests_gl_engine
    UnittestGLContext.h UnittestGLContext.cpp
    framebuffer.cpp
    uniformbuffer.cpp
    texture.cpp
)

target_sources(unittests_gl_engine
  PRIVATE
    tile_funs.cpp
    # tile_drawing.cpp
)

qt_add_resources(unittests_gl_engine "gl_test_data"
    PREFIX "/gl_test_data"
    BASE ${CMAKE_CURRENT_SOURCE_DIR}/data/
    FILES
    data/vectortile_openmaptile_18_140276_92195_kals.pbf
    data/vectortile_mountain_14_8781_5760.pbf
    data/vectortile_skilift_18_140269_92190.pbf
    data/vectortile_highway.png
    data/vectortile_landstrasse_view.png
    data/vectortile_grossglockner_view.png
    data/vectortile_straight_polygon.png
    data/vectortile_straight_thin_line.png
    data/vectortile_vienna_uniform.png
    data/vectortile_vienna_view.png
    data/vectortile_weichtalhaus_view.png
)

target_link_libraries(unittests_gl_engine PUBLIC gl_engine Qt::Test)

target_compile_definitions(unittests_gl_engine PUBLIC "ALP_GL_TEST_DATA_DIR=\":/gl_test_data/\"")

