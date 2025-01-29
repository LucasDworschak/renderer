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

template <typename Func>
std::map<std::string, uint32_t> parse_tile(
    const Style& s, const mapbox::vector_tile::buffer& tile, int zoom, Func key_generator, std::unordered_set<std::string> skipped_layers, bool check_valid_style)
{
    std::map<std::string, uint32_t> feature_to_style;

    for (const auto& layer_name : tile.layerNames()) {
        // qDebug() << layer_name;

        if (skipped_layers.contains(layer_name))
            continue;

        const mapbox::vector_tile::layer layer = tile.getLayer(layer_name);
        std::size_t feature_count = layer.featureCount();

        for (std::size_t i = 0; i < feature_count; ++i) {
            const auto feature = mapbox::vector_tile::feature(layer.getFeature(i), layer);

            if (feature.getType() == mapbox::vector_tile::GeomType::POINT || mapbox::vector_tile::GeomType::UNKNOWN)
                continue; // we are only interested in polygon and line strings

            const std::string type = (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) ? "line" : "fill";

            const auto style = s.layer_style_index(layer_name, type, zoom, feature);

            const auto key = type + "__" + key_generator(layer_name, feature);

            if (feature_to_style.contains(key)) // make sure that all features with the same key share the style (if not we have to expand how we generate the key)
            {
                if (feature_to_style[key] != style)
                    qDebug() << "key does already exist with diferent value: " << key;
                CHECK(feature_to_style[key] == style);
            } else
                feature_to_style[key] = style; // add the style to the map

            if (check_valid_style)
                CHECK(style != -1u); // make sure that every style is available in the stylesheet
        }
    }

    return feature_to_style;
}

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
        CHECK(s.parse_color("#9CF") == 0x99CCFFFF);
        CHECK(s.parse_color("hsl(36, 6%, 74%)") == 0XC1BDB9FF);
    }

    SECTION("Style parsing openmaptile")
    {

        Style s(":/vectorlayerstyles/openstreetmap.json");
        s.load();

        QString filepath = QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer/vectortile_openmaptile_13_4412_2893.pbf");
        QFile file(filepath);
        file.open(QIODevice::ReadOnly | QIODevice::Unbuffered);
        QByteArray data = file.readAll();
        const auto zoom = 13;

        const auto d = data.toStdString();
        const mapbox::vector_tile::buffer tile(d);

        auto key_generator = [](std::string layer_name, mapbox::vector_tile::feature feature) {
            std::string out = "";
            for (auto& prop : feature.getProperties()) {
                if (prop.first == "class" || prop.first == "subclass" || prop.first == "id" || prop.first == "mvt_id" || prop.first.starts_with("name"))
                    continue;
                out += "__" + std::visit(nucleus::vector_tile::util::string_print_visitor, prop.second).toStdString();
            }

            const auto class_name = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["class"]).toStdString();
            const auto subclass_name = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["subclass"]).toStdString();
            return layer_name + "__" + class_name + "__" + subclass_name + out;
        };

        auto skipped_layers = std::unordered_set<std::string> { "transportation_name", "water_name" };

        auto feature_to_style = parse_tile(s, tile, zoom, key_generator, skipped_layers, false);

        // check if the color stored int he style buffer points to the correct color in the stylesheet
        const auto fill_style_buffer = s.style_buffer().fill_styles->buffer();
        const auto line_style_buffer = s.style_buffer().line_styles->buffer();

        // TODO here somewhere -> exception is thrown ...

        CHECK(fill_style_buffer[feature_to_style["fill__landuse__residential__null"] * constants::style_data_size] == s.parse_color("#e0dfdf"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__motorway__null__paved__1__1__1"] * constants::style_data_size] == s.parse_color("#dc2a67"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__motorway__null__paved__1__bridge__1__1__1"] * constants::style_data_size] == s.parse_color("#dc2a67"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__grass__park"] * constants::style_data_size] == s.parse_color("#c8facc"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__rock__bare_rock"] * constants::style_data_size] == s.parse_color("#eee5dc"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__grass__scrub"] * constants::style_data_size] == s.parse_color("#c8d7ab"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__farmland__farmland"] * constants::style_data_size] == s.parse_color("#eef0d5"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__motorway__null__paved__-1__tunnel__1__1"] * constants::style_data_size] == s.parse_color("#c24e6b"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__grass__grass"] * constants::style_data_size] == s.parse_color("#cdebb0"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__grass__grassland"] * constants::style_data_size] == s.parse_color("#cdebb0"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__grass__meadow"] * constants::style_data_size] == s.parse_color("#cdebb0"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__wetland__wetland"] * constants::style_data_size] == s.parse_color("#add19e"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__wood__forest"] * constants::style_data_size] == s.parse_color("#add19e"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__wood__wood"] * constants::style_data_size] == s.parse_color("#add19e"));
        CHECK(fill_style_buffer[feature_to_style["fill__landcover__farmland__vineyard"] * constants::style_data_size] == s.parse_color("#aedfa3"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__path__path__1__no__no"] * constants::style_data_size] == s.parse_color("#fa8072"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__path__path__1__no__no__unpaved"] * constants::style_data_size] == s.parse_color("#fa8072"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__path__path__2__no__no"] * constants::style_data_size] == s.parse_color("#fa8072"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__path__path__2__no__yes__no__no"] * constants::style_data_size] == s.parse_color("#fa8072"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__path__path__no__no"] * constants::style_data_size] == s.parse_color("#fa8072"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__path__path__no__no__no__unpaved__no"] * constants::style_data_size] == s.parse_color("#fa8072"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__path__path__no__unpaved__no"] * constants::style_data_size] == s.parse_color("#fa8072"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__path__path__no__unpaved__yes"] * constants::style_data_size] == s.parse_color("#fa8072"));
        CHECK(fill_style_buffer[feature_to_style["fill__water__pond__null__1"] * constants::style_data_size] == s.parse_color("rgba(172, 218, 251, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__paved__-1__tunnel"] * constants::style_data_size] == s.parse_color("#fff"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__tunnel__-1"] * constants::style_data_size] == s.parse_color("#fff"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__tertiary__null"] * constants::style_data_size] == s.parse_color("#fff"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__tertiary__null__paved"] * constants::style_data_size] == s.parse_color("#fff"));
        CHECK(fill_style_buffer[feature_to_style["fill__building__null__null"] * constants::style_data_size] == s.parse_color("#d9d0c9"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__1__alley"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__1__no"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__1__unpaved"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__alley"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__no"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__paved"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__paved__no__no__no"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__unpaved"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__service__null__yes"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__track__null"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__track__null__no__yes__unpaved__no"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__track__null__no__yes__unpaved__yes"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__track__null__unpaved"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__track__null__unpaved__1__no__yes__no"] * constants::style_data_size] == s.parse_color("#bbbbbb"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__primary__null__paved__1__1"] * constants::style_data_size] == s.parse_color("#a06b00"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__raceway__null__unpaved"] * constants::style_data_size] == s.parse_color("rgba(254, 190, 200, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__1"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__1__paved"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__no__paved__no__1"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__paved"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__paved__yes__no"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__unpaved"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__yes__paved"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__yes__paved__yes"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__yes__paved__yes__ford"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__yes__yes"] * constants::style_data_size] == s.parse_color("rgba(255, 255, 255, 1)"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__secondary__null__1__paved"] * constants::style_data_size] == s.parse_color("#f7fabf"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__secondary__null__paved"] * constants::style_data_size] == s.parse_color("#f7fabf"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__primary__null__bridge__1"] * constants::style_data_size] == s.parse_color("#000000"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__primary__null__paved__1__bridge"] * constants::style_data_size] == s.parse_color("#000000"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__secondary__null__paved__1__bridge"] * constants::style_data_size] == s.parse_color("#000000"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__bridge__1"] * constants::style_data_size] == s.parse_color("#fff"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__paved__1__bridge"] * constants::style_data_size] == s.parse_color("#fff"));
        CHECK(line_style_buffer[feature_to_style["line__transportation__minor__null__paved__yes__yes__1__bridge"] * constants::style_data_size] == s.parse_color("#fff"));

        // // DEBUG show all keys to styles
        // for (const auto& el : feature_to_style) {

        //     if (el.second == -1u)
        //         continue;

        //     // easy to renew check:
        //     std::cout << el.second << "\tCHECK(fill_style_buffer[feature_to_style[\"" << el.first << "\"] * constants::style_data_size] == s.parse_color(\"\"));" << std::endl;

        //     // simple output to check:
        //     // qDebug() << el.first << " " << el.second;
        // }
    }

    SECTION("Style parsing basemap")
    {
        Style s(":/vectorlayerstyles/basemap.json");
        s.load();

        QString filepath = QString("%1%2").arg(ALP_TEST_DATA_DIR, "vector_layer/vectortile_13_4412_2893.pbf");
        QFile file(filepath);
        file.open(QIODevice::ReadOnly | QIODevice::Unbuffered);
        QByteArray data = file.readAll();
        const auto zoom = 13;

        const auto d = data.toStdString();
        const mapbox::vector_tile::buffer tile(d);

        auto key_generator = [](std::string layer_name, mapbox::vector_tile::feature feature) {
            const auto symbol_property = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["_symbol"]).toStdString();
            return layer_name + "__" + symbol_property;
        };

        auto skipped_layers = std::unordered_set<std::string> { "GIP_L_GIP_144/label" }; // for some reason this is has a line geom, but is actually a text, which we do not currently parse

        auto feature_to_style = parse_tile(s, tile, zoom, key_generator, skipped_layers, true);

        // check if the color stored int he style buffer points to the correct color in the stylesheet
        const auto fill_style_buffer = s.style_buffer().fill_styles->buffer();
        const auto line_style_buffer = s.style_buffer().line_styles->buffer();

        CHECK(line_style_buffer[feature_to_style["line__BEV_BEZIRK_L_BEZIRKSGRENZE__0"] * constants::style_data_size] == s.parse_color("#EAE0EF"));
        CHECK(line_style_buffer[feature_to_style["line__BEV_GEMEINDE_L_GEMEINDEGRENZE__0"] * constants::style_data_size] == s.parse_color("#EAE0EF"));
        CHECK(fill_style_buffer[feature_to_style["fill__GEBAEUDE_F_AGG__0"] * constants::style_data_size] == s.parse_color("#EDCACA"));
        CHECK(fill_style_buffer[feature_to_style["fill__GEBAEUDE_F_AGG__1"] * constants::style_data_size] == s.parse_color("#E6B8B8"));
        CHECK(fill_style_buffer[feature_to_style["fill__GEWAESSER_F_GEWF__1"] * constants::style_data_size] == s.parse_color("#B3D9FF"));
        CHECK(fill_style_buffer[feature_to_style["fill__GEWAESSER_F_GEWF__3"] * constants::style_data_size] == s.parse_color("#B3D9FF"));
        CHECK(line_style_buffer[feature_to_style["line__GEWAESSER_L_GEWL __4"] * constants::style_data_size] == s.parse_color("#B3D9FF"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_BAUWERK_L_BRÜCKE__0"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_BAUWERK_L_BRÜCKE__1"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_BAUWERK_L_BRÜCKE__3"] * constants::style_data_size] == s.parse_color("#CDAA66"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_BAUWERK_L_TUNNEL_BRUNNENCLA__0"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_L_GIP_144__0"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_L_GIP_144__1"] * constants::style_data_size] == s.parse_color("#CD8966"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_L_GIP_144__3"] * constants::style_data_size] == s.parse_color("#CDAA66"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_L_GIP_144__4"] * constants::style_data_size] == s.parse_color("#B2B2B2"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_L_GIP_144__5"] * constants::style_data_size] == s.parse_color("#B2B2B2"));
        CHECK(line_style_buffer[feature_to_style["line__GIP_L_GIP_144__6"] * constants::style_data_size] == s.parse_color("#B2B2B2"));
        CHECK(fill_style_buffer[feature_to_style["fill__NUTZUNG_L15_12__0"] * constants::style_data_size] == s.parse_color("#EFEBE9"));
        CHECK(fill_style_buffer[feature_to_style["fill__NUTZUNG_L15_12__1"] * constants::style_data_size] == s.parse_color("rgba(136,204,102,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["fill__NUTZUNG_L15_12__2"] * constants::style_data_size] == s.parse_color("rgba(235,255,170,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["fill__NUTZUNG_L15_12__3"] * constants::style_data_size] == s.parse_color("rgba(71,179,18,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["fill__NUTZUNG_L15_12__5"] * constants::style_data_size] == s.parse_color("rgba(255,255,255,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["fill__NUTZUNG_L15_12__6"] * constants::style_data_size] == s.parse_color("rgba(163,255,115,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["fill__NUTZUNG_L15_12__7"] * constants::style_data_size] == s.parse_color("rgba(153,125,77,0.25)"));
        CHECK(fill_style_buffer[feature_to_style["fill__NUTZUNG_L15_12__8"] * constants::style_data_size] == s.parse_color("rgba(102,153,77,0.25)"));

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
