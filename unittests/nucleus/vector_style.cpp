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

#include <nucleus/tile/TileLoadService.h>
#include <nucleus/tile/conversion.h>
#include <nucleus/tile/utils.h>
#include <nucleus/utils/bit_coding.h>
#include <radix/tile.h>

#include "nucleus/vector_layer/Style.h"
#include "nucleus/vector_layer/StyleExpression.h"
#include "nucleus/vector_layer/constants.h"
#include "nucleus/vector_tile/util.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

#include "nucleus/utils/rasterizer.h"

#include <mapbox/vector_tile.hpp>

using namespace nucleus::vector_layer;

TEST_CASE("nucleus/vector_style")
{

    SECTION("Style parsing color")
    {
        Style s("");
        CHECK(s.parse_color("rgba(134,179,1,0.5)") == 0x86B3017F);
        CHECK(s.parse_color("rgba(134,179,1,1)") == 0x86B301FF);
        CHECK(s.parse_color("rgb(134,179,1)") == 0x86B301FF);
        CHECK(s.parse_color("rgb(50,0,100)") == 0x320064FF);
        CHECK(s.parse_color("#FF55AA") == 0xFF55AAFF);
        CHECK(s.parse_color("#ABCDEF45") == 0xABCDEF45);
    }

    SECTION("Style parsing basemap")
    {
        Style s(":/vectorlayerstyles/basemap.json");
        s.load();

        QString filepath = QString("%1%2").arg(ALP_TEST_DATA_DIR, "vectortile_13_4412_2893.pbf");
        QFile file(filepath);
        file.open(QIODevice::ReadOnly | QIODevice::Unbuffered);
        QByteArray data = file.readAll();
        const auto zoom = 13;

        const auto d = data.toStdString();
        const mapbox::vector_tile::buffer tile(d);

        std::map<std::string, uint32_t> feature_to_style;

        for (const auto& layer_name : tile.layerNames()) {
            // qDebug() << layer_name;

            if (layer_name == "GIP_L_GIP_144/label")
                continue; // for some reason this is has a line geom, but is actually a text, which we do not currently parse

            const mapbox::vector_tile::layer layer = tile.getLayer(layer_name);
            std::size_t feature_count = layer.featureCount();

            for (std::size_t i = 0; i < feature_count; ++i) {                
                const auto feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

                if (feature.getType() == mapbox::vector_tile::GeomType::POINT || mapbox::vector_tile::GeomType::UNKNOWN)
                    continue; // we are only interested in polygon and line strings

                const auto style = s.layer_style_index(layer_name, zoom, feature);

                const auto symbol_property = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["_symbol"]).toStdString();
                const auto key = layer_name + "__" + symbol_property;
                if (feature_to_style.contains(key)) // make sure that all features with the same key share the style (if not we have to expand how we generate the key)
                    CHECK(feature_to_style[key] == style);
                else
                    feature_to_style[key] = style; // add the style to the map

                CHECK(style != -1u); // make sure that every style is available in the stylesheet
            }
        }

        // check if the color stored int he style buffer points to the correct color in the stylesheet
        const auto fill_style_buffer = s.style_buffer().fill_styles->buffer();
        const auto line_style_buffer = s.style_buffer().line_styles->buffer();

        CHECK(line_style_buffer[feature_to_style["BEV_BEZIRK_L_BEZIRKSGRENZE__0"] * constants::style_data_size] == s.parse_color("#EAE0EF"));
        CHECK(line_style_buffer[feature_to_style["BEV_GEMEINDE_L_GEMEINDEGRENZE__0"] * constants::style_data_size] == s.parse_color("#EAE0EF"));
        CHECK(fill_style_buffer[feature_to_style["GEBAEUDE_F_AGG__0"] * constants::style_data_size] == s.parse_color("#EDCACA"));
        CHECK(fill_style_buffer[feature_to_style["GEBAEUDE_F_AGG__1"] * constants::style_data_size] == s.parse_color("#E6B8B8"));
        CHECK(fill_style_buffer[feature_to_style["GEWAESSER_F_GEWF__1"] * constants::style_data_size] == s.parse_color("#B3D9FF"));
        CHECK(fill_style_buffer[feature_to_style["GEWAESSER_F_GEWF__3"] * constants::style_data_size] == s.parse_color("#B3D9FF"));
        CHECK(line_style_buffer[feature_to_style["GEWAESSER_L_GEWL __4"] * constants::style_data_size] == s.parse_color("#B3D9FF"));
        CHECK(line_style_buffer[feature_to_style["GIP_BAUWERK_L_BRÜCKE__0"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["GIP_BAUWERK_L_BRÜCKE__1"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["GIP_BAUWERK_L_BRÜCKE__3"] * constants::style_data_size] == s.parse_color("#CDAA66"));
        CHECK(line_style_buffer[feature_to_style["GIP_BAUWERK_L_TUNNEL_BRUNNENCLA__0"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["GIP_L_GIP_144__0"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["GIP_L_GIP_144__1"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["GIP_L_GIP_144__3"] * constants::style_data_size] == s.parse_color("#CDAA66"));
        CHECK(line_style_buffer[feature_to_style["GIP_L_GIP_144__4"] * constants::style_data_size] == s.parse_color("#B2B2B2"));
        CHECK(line_style_buffer[feature_to_style["GIP_L_GIP_144__5"] * constants::style_data_size] == s.parse_color("#B2B2B2"));
        CHECK(line_style_buffer[feature_to_style["GIP_L_GIP_144__6"] * constants::style_data_size] == s.parse_color("#B2B2B2"));
        CHECK(fill_style_buffer[feature_to_style["NUTZUNG_L15_12__0"] * constants::style_data_size] == s.parse_color("#EFEBE9"));
        CHECK(fill_style_buffer[feature_to_style["NUTZUNG_L15_12__1"] * constants::style_data_size] == s.parse_color("rgba(136,204,102,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["NUTZUNG_L15_12__2"] * constants::style_data_size] == s.parse_color("rgba(235,255,170,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["NUTZUNG_L15_12__3"] * constants::style_data_size] == s.parse_color("rgba(71,179,18,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["NUTZUNG_L15_12__5"] * constants::style_data_size] == s.parse_color("rgba(255,255,255,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["NUTZUNG_L15_12__6"] * constants::style_data_size] == s.parse_color("rgba(163,255,115,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["NUTZUNG_L15_12__7"] * constants::style_data_size] == s.parse_color("rgba(153,125,77,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["NUTZUNG_L15_12__8"] * constants::style_data_size] == s.parse_color("rgba(102,153,77,0.25)"));

        // // DEBUG show all keys to styles
        // for (const auto& el : feature_to_style) {
        //     // easy to renew check:
        //     std::cout << "CHECK(feature_to_style[\"" << el.first << "\"] == " << el.second << ");" << std::endl;
        //     // simple output to check:
        //     // std::cout << el.first << " " << el.second << std::endl;
        // }
    }
}

TEST_CASE("nucleus/vector_style benchmarks")
{
    BENCHMARK("Style parsing basemap")
    {
        Style s(":/vectorlayerstyles/basemap.json");
        s.load();
    };
    // BENCHMARK("triangulize polygons")
    // {
    //     const std::vector<glm::vec2> polygon_points = { glm::vec2(10.5, 10.5), glm::vec2(30.5, 10.5), glm::vec2(50.5, 50.5), glm::vec2(10.5, 30.5) };
    //     const auto edges = nucleus::utils::rasterizer::generate_neighbour_edges(polygon_points);
    //     nucleus::utils::rasterizer::triangulize(polygon_points, edges);
    // };
}
