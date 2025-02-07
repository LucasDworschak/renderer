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

#include <iomanip>

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
            const auto key = type + "__" + key_generator(layer_name, feature);

            const auto style = s.indices(layer_name, type, zoom, feature);

            if (check_valid_style)
                CHECK(style.first != -1u); // make sure that every style is available in the stylesheet

            if (style.first == -1u)
                continue;

            if (feature_to_style.contains(key)) // make sure that all features with the same key share the style (if not we have to expand how we generate the key)
            {
                if (feature_to_style[key] != style.first)
                    qDebug() << "key does already exist with diferent value: " << key;
                CHECK(feature_to_style[key] == style.first);
            } else
                feature_to_style[key] = style.first; // add the style to the map
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

    SECTION("Style expand")
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

        QJsonArray expanded_layers = Style::expand(layers);

        CHECK(layers.size() == 208); // makes sure that the input file is still the same
        CHECK(expanded_layers.size() == 312);

        // // DEBUG view what is written in expanded layers
        // QFile out_file("expanded.style.json");
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

        CHECK(style_buffer[feature_to_style.at("fill__building__null__null")].x == s.parse_color("#d9d0c9ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__farmland__farmland")].x == s.parse_color("#eef0d5ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__farmland__vineyard")].x == s.parse_color("#aedfa3ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__grass")].x == s.parse_color("#cdebb0ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__grassland")].x == s.parse_color("#cdebb0ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__meadow")].x == s.parse_color("#cdebb0ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__park")].x == s.parse_color("#c8faccff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__grass__scrub")].x == s.parse_color("#c8d7abff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__rock__bare_rock")].x == s.parse_color("#eee5dcff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wetland__wetland")].x == s.parse_color("#add19eff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wood__forest")].x == s.parse_color("#add19eff"));
        CHECK(style_buffer[feature_to_style.at("fill__landcover__wood__wood")].x == s.parse_color("#add19eff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__commercial__null")].x == s.parse_color("#f2dad9ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__industrial__null")].x == s.parse_color("#ebdbe8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__pitch__null")].x == s.parse_color("#aae0cbff"));
        CHECK(style_buffer[feature_to_style.at("fill__landuse__residential__null")].x == s.parse_color("#e0dfdfff"));
        CHECK(style_buffer[feature_to_style.at("fill__transportation__bridge__null__bridge__1")].x == s.parse_color("#b8b8b8ff"));
        CHECK(style_buffer[feature_to_style.at("fill__transportation__pier__null")].x == s.parse_color("#f6f1e5ff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__lake__null__0")].x == s.parse_color("#aad3dfff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__pond__null__0")].x == s.parse_color("#aad3dfff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__pond__null__1")].x == s.parse_color("#acdafbd8"));
        CHECK(style_buffer[feature_to_style.at("fill__water__river__null__0")].x == s.parse_color("#aad3dfff"));
        CHECK(style_buffer[feature_to_style.at("fill__water__swimming_pool__null__0")].x == s.parse_color("#aad3dfff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__1__paved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__bridge__1__paved")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__no__paved__1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__no__yes__paved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__paved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__tunnel__-1")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__tunnel__-1__paved")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__unpaved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__bridge__yes__1__paved")].x == s.parse_color("#ffffffff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__ford__yes__paved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__paved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__minor__null__yes__yes__paved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__1__paved__1")].x == s.parse_color("#dc2a67ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__1__paved__1")].x == s.parse_color("#dc2a67ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__1__paved__1")].x == s.parse_color("#dc2a67ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__bridge__1__1__paved__1")].x == s.parse_color("#000000ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__motorway__null__tunnel__-1__1__paved__1")].x == s.parse_color("#c24e6bff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__1")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__1__unpaved")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__2")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__no__no__unpaved")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__unpaved")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__no__yes__no__2")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__unpaved")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__no__yes")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__unpaved")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__path__path__yes__no__unpaved")].x == s.parse_color("#fa8072ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null")].x == s.parse_color("#a07400ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__1__1__paved")].x == s.parse_color("#a06b00ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1")].x == s.parse_color("#000000ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__bridge__1__paved")].x == s.parse_color("#000000ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__primary__null__paved")].x == s.parse_color("#a07400ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__raceway__null__unpaved")].x == s.parse_color("#febec8ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__1__paved")].x == s.parse_color("#707d05ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__bridge__1__paved")].x == s.parse_color("#000000ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__secondary__null__paved")].x == s.parse_color("#707d05ff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__alley")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__1__unpaved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__alley")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__1")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__no__no__no__paved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__paved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__unpaved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__service__null__yes")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null")].x == s.parse_color("#8f8f8fff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__tertiary__null__paved")].x == s.parse_color("#8f8f8fff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__1__unpaved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__no__yes__no__unpaved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__unpaved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__transportation__track__null__yes__yes__no__unpaved")].x == s.parse_color("#bbbbbbff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__river__null__0")].x == s.parse_color("#a0c8f0ff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__stream__null__0")].x == s.parse_color("#a0c8f0ff"));
        CHECK(style_buffer[feature_to_style.at("line__waterway__stream__null__1")].x == s.parse_color("#a0c8f0ff"));

        // // DEBUG show all keys to styles
        // for (const auto& el : feature_to_style) {

        //     if (el.second == -1u)
        //         continue;

        //     // easy to renew check:
        //     std::cout << "CHECK(style_buffer[feature_to_style.at(\"" << el.first << "\")].x == s.parse_color(\"#" << std::hex << style_buffer[el.second].x << "\"));" << std::endl;

        //     // simple output to check:
        //     // qDebug() << el.first << " " << el.second;
        // }
        // std::cout << std::endl << std::endl;
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

        CHECK(style_buffer[feature_to_style.at("line__BEV_BEZIRK_L_BEZIRKSGRENZE__0")].x == s.parse_color("#EAE0EF"));
        CHECK(style_buffer[feature_to_style.at("line__BEV_GEMEINDE_L_GEMEINDEGRENZE__0")].x == s.parse_color("#EAE0EF"));
        CHECK(style_buffer[feature_to_style.at("fill__GEBAEUDE_F_AGG__0")].x == s.parse_color("#EDCACA"));
        CHECK(style_buffer[feature_to_style.at("fill__GEBAEUDE_F_AGG__1")].x == s.parse_color("#E6B8B8"));
        CHECK(style_buffer[feature_to_style.at("fill__GEWAESSER_F_GEWF__1")].x == s.parse_color("#B3D9FF"));
        CHECK(style_buffer[feature_to_style.at("fill__GEWAESSER_F_GEWF__3")].x == s.parse_color("#B3D9FF"));
        CHECK(style_buffer[feature_to_style.at("line__GEWAESSER_L_GEWL __4")].x == s.parse_color("#B3D9FF"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__0")].x == s.parse_color("#CD8966"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__1")].x == s.parse_color("#CD8966"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_BRÜCKE__3")].x == s.parse_color("#CDAA66"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_BAUWERK_L_TUNNEL_BRUNNENCLA__0")].x == s.parse_color("#CD8966"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__0")].x == s.parse_color("#CD8966"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__1")].x == s.parse_color("#CD8966"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__3")].x == s.parse_color("#CDAA66"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__4")].x == s.parse_color("#B2B2B2"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__5")].x == s.parse_color("#B2B2B2"));
        CHECK(style_buffer[feature_to_style.at("line__GIP_L_GIP_144__6")].x == s.parse_color("#B2B2B2"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__0")].x == s.parse_color("#EFEBE9"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__1")].x == s.parse_color("rgba(136,204,102,0.25)"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__2")].x == s.parse_color("rgba(235,255,170,0.25)"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__3")].x == s.parse_color("rgba(71,179,18,0.25)"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__5")].x == s.parse_color("rgba(255,255,255,0.25)"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__6")].x == s.parse_color("rgba(163,255,115,0.25)"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__7")].x == s.parse_color("rgba(153,125,77,0.25)"));
        CHECK(style_buffer[feature_to_style.at("fill__NUTZUNG_L15_12__8")].x == s.parse_color("rgba(102,153,77,0.25)"));

        // // DEBUG show all keys to styles
        // for (const auto& el : feature_to_style) {

        //     if (el.second == -1u)
        //         continue;

        //     // easy to renew check:
        //     std::cout << "CHECK(style_buffer[feature_to_style.at(\"" << el.first << "\")].x == s.parse_color(\"#" << std::hex << style_buffer[el.second].x << "\"));" << std::endl;

        //     // simple output to check:
        //     // qDebug() << el.first << " " << el.second;
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
}
