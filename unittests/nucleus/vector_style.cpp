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

#include "nucleus/vector_layer/Preprocessor.h"
#include "nucleus/vector_layer/Style.h"
#include "nucleus/vector_layer/StyleExpander.h"
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

    const auto style_buffer = s.styles()->buffer();

    std::array<int, nucleus::vector_layer::constants::max_style_expression_keys> temp_values;

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

            const auto type = (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) ? 0 : 1;
            const std::string type_string = (feature.getType() == mapbox::vector_tile::GeomType::LINESTRING) ? "line" : "fill";
            const auto key = type_string + "__" + key_generator(layer_name, feature);

            // qDebug() << key;

            auto style_indices = s.indices(layer_name, type, zoom, feature, &temp_values);
            style_indices = nucleus::vector_layer::Preprocessor::simplify_styles(&style_indices, style_buffer);

            if (check_valid_style)
                CHECK(style_indices.size() > 0); // make sure that every style is available in the stylesheet

            if (style_indices.size() == 0)
                continue;

            int index = 0;

            for (const auto& style : style_indices) {
                const uint32_t style_index = style.first;
                std::string indexed_key = key + "_" + std::to_string(index);
                if (feature_to_style.contains(indexed_key)) // make sure that all features with the same key share the style (if not we have to expand how we generate the key)
                {
                    if (feature_to_style[indexed_key] != style_index)
                        qDebug() << "key does already exist with diferent value: " << indexed_key;
                    CHECK(feature_to_style[indexed_key] == style_index);
                } else
                    feature_to_style[indexed_key] = style_index; // add the style to the map

                index++;
            }
        }
    }

    return feature_to_style;
}

void create_debug_filter_checks(std::map<std::string, uint32_t> feature_to_style, const std::vector<glm::u32vec4> style_buffer)
{
    std::cout << "CHECK(feature_to_style.size()==" << std::dec << feature_to_style.size() << ");" << std::endl;
    for (const auto& el : feature_to_style) {

        if (el.second == -1u)
            continue;

        // easy to renew check:
        std::cout << "CHECK(style_buffer[feature_to_style.at(\"" << el.first << "\")].x == s.parse_color(\"#" << std::hex << style_buffer[el.second].x << "\"));" << std::endl;
        // std::cout << "CHECK((style_buffer[feature_to_style.at(\"" << el.first << "\")].z  / nucleus::vector_layer::constants::style_precision) == " << std::dec
        // << (float(style_buffer[el.second].z) / nucleus::vector_layer::constants::style_precision) << ");" << std::endl;

        // simple output to check:
        // qDebug() << el.first << " " << el.second;
    }
    std::cout << std::endl << std::endl;
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

    SECTION("Premultiply alpha")
    {
        // full alpha
        CHECK(Style::premultiply_alpha(0xFFFFFFFF) == 0xFFFFFFFF);
        CHECK(Style::premultiply_alpha(0xAABBCCFF) == 0xAABBCCFF);
        CHECK(Style::premultiply_alpha(0x334455FF) == 0x334455FF);

        // half alpha
        CHECK(Style::premultiply_alpha(0xFFFFFF80) == 0x80808080);
        CHECK(Style::premultiply_alpha(0xAABBCC80) == 0x555D6680);
        CHECK(Style::premultiply_alpha(0x33445580) == 0x19222A80);

        CHECK(Style::premultiply_alpha(0xFFFFFF15) == 0x15151515);
        CHECK(Style::premultiply_alpha(0xFFFFFFCC) == 0xCCCCCCCC);
        CHECK(Style::premultiply_alpha(0xAA558815) == 0x0E070B15);
        CHECK(Style::premultiply_alpha(0x991188CC) == 0x7A0D6CCC);
    }

    SECTION("Simple style parsing")
    {
        // define your own corner case styles in the test-style.json
        // and verify how they are parsed
        // you might want to display the style_index to layername qdebug in style.cpp
        Style s(":/test_data/vector_layer/test-style.json");
        s.load();

        const auto style_buffer = s.styles()->buffer();
        const auto line_multipliers = nucleus::vector_layer::constants::line_width_multiplier * nucleus::vector_layer::constants::style_precision;

        // note order in style buffer may not be the order defined in style.json
        // since we are using a key comparer and a map the order should still be the same between compilers/os
        // the layer_index order is however preserved (just not tested in this testcase)

        // "opacity outside of zoom range"
        CHECK(style_buffer[0].x == 0xbbbbbbff); // z 13-3
        CHECK(style_buffer[1].x == 0xbbbbbbff); // z 13-2
        CHECK(style_buffer[2].x == 0xbbbbbbff); // z 13-1
        CHECK(style_buffer[3].x == 0xbbbbbbff); // z 13
        CHECK(style_buffer[4].x == 0xbbbbbbff); // z 14
        CHECK(style_buffer[5].x == 0xbbbbbbff); // z 15
        CHECK(style_buffer[6].x == 0xbbbbbbff); // z 16
        CHECK(style_buffer[7].x == 0xbbbbbbff); // z 17
        CHECK(style_buffer[8].x == 0xbbbbbbff); // z 18

        // reuse style if no blending
        CHECK(style_buffer[9].x == 0xaaaaaaff); // make sure that we are in the right style instruction here (not using other style)
        CHECK(style_buffer[9].z == 8 * line_multipliers); // z 11-3
        CHECK(style_buffer[10].z == 8 * line_multipliers); // z 11-2
        CHECK(style_buffer[11].z == 8 * line_multipliers); // z 11-1
        CHECK(style_buffer[11].x == 0xaaaaaaff); // z 11-1 // color
        CHECK(style_buffer[12].x == 0xaaaaaaff); // z 11 // color
        CHECK(style_buffer[12].z == 8 * line_multipliers); // z 11
        CHECK(style_buffer[13].z == 8 * line_multipliers); // z 12
        CHECK(style_buffer[14].z == 8 * line_multipliers); // z 13
        CHECK(style_buffer[15].z == 9 * line_multipliers); // z 14
        CHECK(style_buffer[16].z == 10 * line_multipliers); // z 15
        CHECK(style_buffer[17].z == 10 * line_multipliers); // z 16
        CHECK(style_buffer[18].z == 10 * line_multipliers); // z 17
        CHECK(style_buffer[19].z == 10 * line_multipliers); // z 18

        CHECK(style_buffer[20].x == -1u); // no data
        CHECK(style_buffer[21].z == -1u); // no data
    }

    SECTION("Style expand openstreetmap")
    {
        QFile file(":/vectorlayerstyles/openstreetmap.json");
        file.open(QIODeviceBase::OpenModeFlag::ReadOnly);

        const auto data = file.readAll();

        if (data.isEmpty()) {
            CHECK(false);
            return;
        }

        QJsonDocument doc = QJsonDocument::fromJson(data);
        QJsonArray layers = doc.object().value("layers").toArray();

        QJsonArray expanded_layers = style_expander::expand(layers);

        CHECK(layers.size() == 208); // makes sure that the input file is still the same
        CHECK(expanded_layers.size() == 312);

        // // DEBUG view what is written in expanded layers
        // QFile out_file("expanded.style.json");
        // out_file.open(QFile::WriteOnly);
        // QJsonDocument out = QJsonDocument(expanded_layers);
        // out_file.write(out.toJson());
    }

    SECTION("Style expand qwant")
    {
        // expanding the following filters
        // landcover-grass (+1)
        // water (+1)
        // bridge-link (+4)
        // bridge-trunk-primary (+2)
        // = adds 8 filters

        QFile file(":/vectorlayerstyles/qwant.json");
        file.open(QIODeviceBase::OpenModeFlag::ReadOnly);

        const auto data = file.readAll();

        if (data.isEmpty()) {
            CHECK(false);
            return;
        }

        QJsonDocument doc = QJsonDocument::fromJson(data);
        QJsonArray layers = doc.object().value("layers").toArray();

        QJsonArray expanded_layers = style_expander::expand(layers);

        CHECK(layers.size() == 115); // makes sure that the input file is still the same
        CHECK(expanded_layers.size() == 123);

        // qDebug() << layers.size();
        // qDebug() << expanded_layers.size();

        // // DEBUG view what is written in expanded layers
        // QFile out_file("qwant-expanded-style.json");
        // out_file.open(QFile::WriteOnly);
        // QJsonDocument out = QJsonDocument(expanded_layers);
        // out_file.write(out.toJson());
    }

    SECTION("Style expand osm-bright")
    {
        // osm-bright style does not expand

        QFile file(":/vectorlayerstyles/osm-bright.json");
        file.open(QIODeviceBase::OpenModeFlag::ReadOnly);

        const auto data = file.readAll();

        if (data.isEmpty()) {
            CHECK(false);
            return;
        }

        QJsonDocument doc = QJsonDocument::fromJson(data);
        QJsonArray layers = doc.object().value("layers").toArray();

        QJsonArray expanded_layers = style_expander::expand(layers);

        CHECK(layers.size() == 128); // makes sure that the input file is still the same
        CHECK(expanded_layers.size() == 128);

        // qDebug() << layers.size();
        // qDebug() << expanded_layers.size();

        // DEBUG view what is written in expanded layers
        // QFile out_file("osm-bright-expanded-style.json");
        // out_file.open(QFile::WriteOnly);
        // QJsonDocument out = QJsonDocument(expanded_layers);
        // out_file.write(out.toJson());
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
            std::vector<std::pair<std::string, std::string>> sorted_props;

            for (auto& prop : feature.getProperties()) {
                if (prop.first == "class" || prop.first == "subclass" || prop.first == "id" || prop.first == "mvt_id" || prop.first.starts_with("name"))
                    continue;
                sorted_props.push_back(std::make_pair(prop.first, std::visit(nucleus::vector_tile::util::string_print_visitor, prop.second).toStdString()));
            }

            std::sort(sorted_props.begin(), sorted_props.end(), [](std::pair<std::string, std::string> a, std::pair<std::string, std::string> b) { return a.first < b.first; });

            for (auto& prop : sorted_props) {
                out += "__" + prop.second;
            }

            const auto class_name = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["class"]).toStdString();
            const auto subclass_name = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["subclass"]).toStdString();
            return layer_name + "__" + class_name + "__" + subclass_name + out;
        };

        auto skipped_layers = std::unordered_set<std::string> { "transportation_name", "water_name" };

        auto feature_to_style = parse_tile(s, tile, zoom, key_generator, skipped_layers, false);

        // check if the color stored int he style buffer points to the correct color in the stylesheet
        const auto style_buffer = s.styles()->buffer();

        CHECK(feature_to_style.size() == 127);
        CHECK(style_buffer[feature_to_style.at("fill__building__null__null_0")].x == s.parse_color("#ded5cfff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__farmland__farmland_0")].x == s.parse_color("#eef0d5ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__farmland__vineyard_0")].x == s.parse_color("#aedfa3ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__grass_0")].x == s.parse_color("#cdebb0ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__grassland_0")].x == s.parse_color("#cdebb0ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__meadow_0")].x == s.parse_color("#cdebb0ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__park_0")].x == s.parse_color("#c8faccff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__scrub_0")].x == s.parse_color("#c8d7abff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__rock__bare_rock_0")].x == s.parse_color("#eee5dcff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wetland__wetland_0")].x == s.parse_color("#add19eff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wood__forest_0")].x == s.parse_color("#add19eff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wood__wood_0")].x == s.parse_color("#add19eff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__commercial__null_0")].x == s.parse_color("#f2dad9ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__industrial__null_0")].x == s.parse_color("#ebdbe8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__pitch__null_0")].x == s.parse_color("#aae0cbff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__residential__null_0")].x == s.parse_color("#e0dfdfff"));
        CHECK(style_buffer[feature_to_style.at("fill__transportation__bridge__null__bridge__1_0")].x == s.parse_color("#b8b8b8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__transportation__pier__null_0")].x == s.parse_color("#f6f1e5ff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__lake__null__0_0")].x == s.parse_color("#aad3dfff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__pond__null__0_0")].x == s.parse_color("#aad3dfff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__pond__null__1_0")].x == s.parse_color("#91b8d4d8"));
        CHECK(style_buffer[feature_to_style.at("fill__water__river__null__0_0")].x == s.parse_color("#aad3dfff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__swimming_pool__null__0_0")].x == s.parse_color("#aad3dfff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1__paved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__no__paved__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__no__paved__1_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__yes__paved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__paved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__tunnel__-1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__tunnel__-1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__unpaved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__bridge__yes__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__ford__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__ford__yes__paved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__paved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes__paved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__1__paved__1_0")].x == s.parse_color("#e892a2ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__1__paved__1_1")].x == s.parse_color("#dc2a67ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__paved__1_0")].x == s.parse_color("#e892a2ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__paved__1_1")].x == s.parse_color("#dc2a67ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__1__paved__1_0")].x == s.parse_color("#e892a2ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__1__paved__1_1")].x == s.parse_color("#dc2a67ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__paved__1_0")].x == s.parse_color("#e892a2ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__paved__1_1")].x == s.parse_color("#000000ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__tunnel__-1__1__paved__1_0")].x == s.parse_color("#f1bcc6ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__tunnel__-1__1__paved__1_1")].x == s.parse_color("#c24e6bff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__1_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__1__unpaved_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__2_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__no__no__unpaved_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__unpaved_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__yes__no__2_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__unpaved_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__yes_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__unpaved_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__yes__no__unpaved_0")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null_0")].x == s.parse_color("#fcd6a4ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null_1")].x == s.parse_color("#a07400ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__1__1__paved_0")].x == s.parse_color("#fcd6a4ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__1__1__paved_1")].x == s.parse_color("#a06b00ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1_0")].x == s.parse_color("#fcd6a4ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1_1")].x == s.parse_color("#000000ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1__paved_0")].x == s.parse_color("#fcd6a4ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1__paved_1")].x == s.parse_color("#000000ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__paved_0")].x == s.parse_color("#fcd6a4ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__paved_1")].x == s.parse_color("#a07400ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__raceway__null__unpaved_0")].x == s.parse_color("#febec8ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__1__paved_0")].x == s.parse_color("#f7fabfff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__1__paved_1")].x == s.parse_color("#707d05ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__bridge__1__paved_0")].x == s.parse_color("#f7fabfff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__bridge__1__paved_1")].x == s.parse_color("#c3bdbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__paved_0")].x == s.parse_color("#f7fabfff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__paved_1")].x == s.parse_color("#707d05ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__alley_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__alley_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__unpaved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__alley_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__alley_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__1_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__no__no__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__no__no__paved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__paved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__unpaved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__yes_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__yes_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null_1")].x == s.parse_color("#8f8f8fff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null__paved_1")].x == s.parse_color("#8f8f8fff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__1__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__1__unpaved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__unpaved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__unpaved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__yes__yes__no__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__yes__yes__no__unpaved_1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__river__null__0_0")].x == s.parse_color("#a0c8f0ff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__stream__null__0_0")].x == s.parse_color("#a0c8f0ff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__stream__null__1_0")].x == s.parse_color("#a0c8f0ff"));

        // DEBUG show all keys to styles
        // create_debug_filter_checks(feature_to_style, style_buffer);
    }

    SECTION("Style parsing qwant")
    {

        Style s(":/vectorlayerstyles/qwant.json");

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
            std::vector<std::pair<std::string, std::string>> sorted_props;

            for (auto& prop : feature.getProperties()) {
                if (prop.first == "class" || prop.first == "subclass" || prop.first == "id" || prop.first == "mvt_id" || prop.first.starts_with("name"))
                    continue;
                sorted_props.push_back(std::make_pair(prop.first, std::visit(nucleus::vector_tile::util::string_print_visitor, prop.second).toStdString()));
            }

            std::sort(sorted_props.begin(), sorted_props.end(), [](std::pair<std::string, std::string> a, std::pair<std::string, std::string> b) { return a.first < b.first; });

            for (auto& prop : sorted_props) {
                out += "__" + prop.second;
            }

            const auto class_name = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["class"]).toStdString();
            const auto subclass_name = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["subclass"]).toStdString();
            return layer_name + "__" + class_name + "__" + subclass_name + out;
        };

        auto skipped_layers = std::unordered_set<std::string> { "transportation_name", "water_name" };

        auto feature_to_style = parse_tile(s, tile, zoom, key_generator, skipped_layers, false);

        // check if the color stored int he style buffer points to the correct color in the stylesheet
        const auto style_buffer = s.styles()->buffer();

        CHECK(feature_to_style.size() == 125);
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__grass_0")].x == s.parse_color("#e0f2d3ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__grassland_0")].x == s.parse_color("#e0f2d3ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__meadow_0")].x == s.parse_color("#e0f2d3ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__park_0")].x == s.parse_color("#e0f2d3ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__scrub_0")].x == s.parse_color("#e0f2d3ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wood__forest_0")].x == s.parse_color("#cae4beff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wood__wood_0")].x == s.parse_color("#cae4beff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__commercial__null_0")].x == s.parse_color("#56524156"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__industrial__null_0")].x == s.parse_color("#56524156"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__pitch__null_0")].x == s.parse_color("#151e1533"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__residential__null_0")].x == s.parse_color("#5554545c"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__zoo__null_0")].x == s.parse_color("#0a02040c"));
        CHECK(style_buffer[feature_to_style.at("fill__transportation__bridge__null__bridge__1_0")].x == s.parse_color("#d4d4d4e5"));
        CHECK(style_buffer[feature_to_style.at("fill__transportation__pier__null_0")].x == s.parse_color("#f8f4f0ff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__lake__null__0_0")].x == s.parse_color("#bbe0fcff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__pond__null__0_0")].x == s.parse_color("#bbe0fcff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__pond__null__1_0")].x == s.parse_color("#bbe0fcff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__river__null__0_0")].x == s.parse_color("#bbe0fcff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__swimming_pool__null__0_0")].x == s.parse_color("#bbe0fcff"));
        CHECK(style_buffer[feature_to_style.at("line__boundary__null__null__6__0__0_0")].x == s.parse_color("#9e9cabff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__no__paved__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__no__paved__1_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__yes__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__tunnel__-1_0")].x == s.parse_color("#cfcdcaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__tunnel__-1__paved_0")].x == s.parse_color("#cfcdcaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__unpaved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__bridge__yes__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__bridge__yes__1__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__ford__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__ford__yes__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__1__paved__1_0")].x == s.parse_color("#ffcc88ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__paved__1_0")].x == s.parse_color("#ffcc88ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__1__paved__1_0")].x == s.parse_color("#ffcc88ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__paved__1_0")].x == s.parse_color("#ffcc88ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__tunnel__-1__1__paved__1_0")].x == s.parse_color("#ffdaa6ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__tunnel__-1__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__1_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__1__unpaved_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__2_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__no__no__unpaved_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__unpaved_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__yes__no__2_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__unpaved_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__yes_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__unpaved_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__yes__no__unpaved_0")].x == s.parse_color("#00000011"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null_0")].x == s.parse_color("#fdeab2ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__1__1__paved_0")].x == s.parse_color("#fdeab2ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__1__1__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1_1")].x == s.parse_color("#eba76bff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1__paved_1")].x == s.parse_color("#eba76bff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__paved_0")].x == s.parse_color("#fdeab2ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__1__paved_0")].x == s.parse_color("#fef1ccff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__1__paved_1")].x == s.parse_color("#fcdc7fff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__bridge__1__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__bridge__1__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__paved_0")].x == s.parse_color("#fef1ccff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__paved_1")].x == s.parse_color("#fcdc7fff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__alley_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__alley_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__unpaved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__alley_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__alley_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__1_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__no__no__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__no__no__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__paved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__unpaved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__yes_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__yes_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null_0")].x == s.parse_color("#fef1ccff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null_1")].x == s.parse_color("#fcdc7fff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null__paved_0")].x == s.parse_color("#fef1ccff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null__paved_1")].x == s.parse_color("#fcdc7fff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__1__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__1__unpaved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__unpaved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__unpaved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__yes__yes__no__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__yes__yes__no__unpaved_1")].x == s.parse_color("#edededff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__river__null__0_0")].x == s.parse_color("#a0c8f0ff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__stream__null__0_0")].x == s.parse_color("#a0c8f0ff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__stream__null__1_0")].x == s.parse_color("#a0c8f0ff"));

        // DEBUG show all keys to styles
        // create_debug_filter_checks(feature_to_style, style_buffer);
    }

    SECTION("Style parsing osm-bright")
    {

        Style s(":/vectorlayerstyles/osm-bright.json");

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
            std::vector<std::pair<std::string, std::string>> sorted_props;

            for (auto& prop : feature.getProperties()) {
                if (prop.first == "class" || prop.first == "subclass" || prop.first == "id" || prop.first == "mvt_id" || prop.first.starts_with("name"))
                    continue;
                sorted_props.push_back(std::make_pair(prop.first, std::visit(nucleus::vector_tile::util::string_print_visitor, prop.second).toStdString()));
            }

            std::sort(sorted_props.begin(), sorted_props.end(), [](std::pair<std::string, std::string> a, std::pair<std::string, std::string> b) { return a.first < b.first; });

            for (auto& prop : sorted_props) {
                out += "__" + prop.second;
            }

            const auto class_name = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["class"]).toStdString();
            const auto subclass_name = std::visit(nucleus::vector_tile::util::string_print_visitor, feature.getProperties()["subclass"]).toStdString();
            return layer_name + "__" + class_name + "__" + subclass_name + out;
        };

        auto skipped_layers = std::unordered_set<std::string> { "transportation_name", "water_name" };

        auto feature_to_style = parse_tile(s, tile, zoom, key_generator, skipped_layers, false);

        // check if the color stored int he style buffer points to the correct color in the stylesheet
        const auto style_buffer = s.styles()->buffer();

        CHECK(feature_to_style.size() == 94);
        // CHECK(style_buffer[feature_to_style.at("fill__building__null__null_0")].x == s.parse_color("#00000000"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__grass_0")].x == s.parse_color("#d8e8c8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__grassland_0")].x == s.parse_color("#d8e8c8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__meadow_0")].x == s.parse_color("#d8e8c8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__park_0")].x == s.parse_color("#d8e8c8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__scrub_0")].x == s.parse_color("#d8e8c8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wood__forest_0")].x == s.parse_color("#0a100619"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wood__wood_0")].x == s.parse_color("#0a100619"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__commercial__null_0")].x == s.parse_color("#372d2d3a"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__industrial__null_0")].x == s.parse_color("#56524156"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__residential__null_0")].x == s.parse_color("#5452515c"));
        CHECK(style_buffer[feature_to_style.at("fill__transportation__bridge__null__bridge__1_0")].x == s.parse_color("#cbcbcbe5"));
        CHECK(style_buffer[feature_to_style.at("fill__transportation__pier__null_0")].x == s.parse_color("#f8f4f0ff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__lake__null__0_0")].x == s.parse_color("#bfd9f2ff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__pond__null__0_0")].x == s.parse_color("#bfd9f2ff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__pond__null__1_0")].x == s.parse_color("#8597a8b2"));
        CHECK(style_buffer[feature_to_style.at("fill__water__river__null__0_0")].x == s.parse_color("#bfd9f2ff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__swimming_pool__null__0_0")].x == s.parse_color("#bfd9f2ff"));
        CHECK(style_buffer[feature_to_style.at("line__boundary__null__null__6__0__0_0")].x == s.parse_color("#9e9cabff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__no__paved__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__tunnel__-1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__tunnel__-1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__bridge__yes__1__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__ford__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__1__paved__1_0")].x == s.parse_color("#ffcc88ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__paved__1_0")].x == s.parse_color("#ffcc88ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__1__paved__1_0")].x == s.parse_color("#ffcc88ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__paved__1_0")].x == s.parse_color("#ffcc88ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__tunnel__-1__1__paved__1_0")].x == s.parse_color("#ffdaa6ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__tunnel__-1__1__paved__1_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__1_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__1__unpaved_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__2_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__no__no__unpaved_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__unpaved_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__yes__no__2_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__unpaved_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__yes_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__unpaved_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__yes__no__unpaved_0")].x == s.parse_color("#ccbbaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__1__1__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__1__1__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1_1")].x == s.parse_color("#eba76bff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1__paved_1")].x == s.parse_color("#eba76bff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__1__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__1__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__bridge__1__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__bridge__1__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__alley_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__alley_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__1_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__no__no__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__paved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__yes_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null__paved_0")].x == s.parse_color("#ffeeaaff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null__paved_1")].x == s.parse_color("#e9ac77ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__1__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__yes__yes__no__unpaved_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__river__null__0_0")].x == s.parse_color("#a0c8f0ff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__stream__null__0_0")].x == s.parse_color("#a0c8f0ff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__stream__null__1_0")].x == s.parse_color("#a0c8f0ff"));

        // DEBUG show all keys to styles
        // create_debug_filter_checks(feature_to_style, style_buffer);
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
        const auto style_buffer = s.styles()->buffer();

        CHECK(feature_to_style.size() == 36);
        CHECK(style_buffer[feature_to_style.at("fill__GEBAEUDE_F_AGG__0_0")].x == s.parse_color("#edcacaff"));
        CHECK(style_buffer[feature_to_style.at("fill__GEBAEUDE_F_AGG__1_0")].x == s.parse_color("#e6b8b8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__GEWAESSER_F_GEWF__1_0")].x == s.parse_color("#b3d9ffff"));
        CHECK(style_buffer[feature_to_style.at("fill__GEWAESSER_F_GEWF__3_0")].x == s.parse_color("#b3d9ffff"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__0_0")].x == s.parse_color("#efebe9ff"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__1_0")].x == s.parse_color("#2132193f"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__2_0")].x == s.parse_color("#3a3f2a3f"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__3_0")].x == s.parse_color("#112c043f"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__5_0")].x == s.parse_color("#3f3f3f3f"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__6_0")].x == s.parse_color("#283f1c3f"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__7_0")].x == s.parse_color("#251e133f"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__8_0")].x == s.parse_color("#1925133f"));
        CHECK(style_buffer[feature_to_style.at("line__BEV_BEZIRK_L_BEZIRKSGRENZE__0_0")].x == s.parse_color("#b094a0ff"));
        CHECK(style_buffer[feature_to_style.at("line__BEV_BEZIRK_L_BEZIRKSGRENZE__0_1")].x == s.parse_color("#eae0efff"));
        CHECK(style_buffer[feature_to_style.at("line__BEV_GEMEINDE_L_GEMEINDEGRENZE__0_0")].x == s.parse_color("#b094a0ff"));
        CHECK(style_buffer[feature_to_style.at("line__BEV_GEMEINDE_L_GEMEINDEGRENZE__0_1")].x == s.parse_color("#eae0efff"));
        CHECK(style_buffer[feature_to_style.at("line__GEWAESSER_L_GEWL __4_0")].x == s.parse_color("#b3d9ffff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__0_0")].x == s.parse_color("#ffd37fff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__0_1")].x == s.parse_color("#cd8966ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__1_0")].x == s.parse_color("#ffd37fff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__1_1")].x == s.parse_color("#cd8966ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__3_0")].x == s.parse_color("#ffff99ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__3_1")].x == s.parse_color("#cdaa66ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_TUNNEL_BRUNNENCLA__0_0")].x == s.parse_color("#feefd8ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_TUNNEL_BRUNNENCLA__0_1")].x == s.parse_color("#cd8966ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__0_0")].x == s.parse_color("#ffd37fff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__0_1")].x == s.parse_color("#cd8966ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__1_0")].x == s.parse_color("#ffd37fff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__1_1")].x == s.parse_color("#cd8966ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__3_0")].x == s.parse_color("#ffff99ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__3_1")].x == s.parse_color("#cdaa66ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__4_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__4_1")].x == s.parse_color("#b2b2b2ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__5_0")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__5_1")].x == s.parse_color("#b2b2b2ff"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__6_0")].x == s.parse_color("#b2b2b2ff"));

        // DEBUG show all keys to styles
        // create_debug_filter_checks(feature_to_style, style_buffer);
    }
}

TEST_CASE("nucleus/vector_style benchmarks")
{
    BENCHMARK("Style parsing basemap")
    {
        Style s(":/vectorlayerstyles/basemap.json");
        s.load();
    };
    BENCHMARK("Style parsing openmaptile")
    {
        Style s(":/vectorlayerstyles/openstreetmap.json");
        s.load();
    };
}
