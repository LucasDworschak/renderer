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

#pragma once

#include <QHash>
#include <QObject>
#include <QString>
#include <QVariantHash>
#include <cstdint>
#include <glm/glm.hpp>
#include <nucleus/tile/types.h>
#include <radix/tile.h>

namespace nucleus::vector_tile {

struct PointOfInterest {
    Q_GADGET
public:
    enum class Type { Unknown = 0, Peak, Settlement, AlpineHut, Webcam, NumberOfElements };
    Q_ENUM(Type)
    uint64_t id = -1;
    Type type = Type::Unknown;
    QString name;
    glm::dvec3 lat_long_alt = glm::dvec3(0);
    glm::dvec3 world_space_pos = glm::dvec3(0);
    float importance = 0;
    QVariantHash attributes;
};

using PointOfInterestCollection = std::vector<PointOfInterest>;
using PointOfInterestCollectionPtr = std::shared_ptr<const PointOfInterestCollection>;

struct PoiTile {
    tile::Id id;
    vector_tile::PointOfInterestCollectionPtr data;
};
static_assert(tile::NamedTile<PoiTile>);
// using PointOfInterestTileCollection = std::unordered_map<tile::Id, PointOfInterestCollectionPtr, tile::Id::Hasher>;

struct Triangle {
    unsigned int top_index;
    unsigned int middle_index;
    unsigned int bottom_index;

    unsigned int style_index;

    bool operator==(const Triangle& other) const = default;
};

// helper for catch2
inline std::ostream& operator<<(std::ostream& os, const Triangle& t)
{
    return os << "{ top: " << std::to_string(t.top_index) << ", middle: " << std::to_string(t.middle_index) << ", bottom: " << std::to_string(t.bottom_index)
              << ", style: " << std::to_string(t.style_index) << " }";
}

struct VectorLayer {
    Q_GADGET
public:
    std::vector<glm::vec2> vertices;
    std::vector<Triangle> triangles;
    std::vector<glm::vec2> vertex_normals;

    std::vector<std::unordered_set<uint8_t>> cell_to_triangle;
};

using VectorLayerCollection = std::unordered_map<tile::Id, VectorLayer, tile::Id::Hasher>;
using VectorLayerPtr = std::shared_ptr<const VectorLayer>;

struct VectorLayerTile {
    tile::Id id;
    vector_tile::VectorLayerPtr data;
};

static_assert(tile::NamedTile<VectorLayerTile>);

using VectorLayerTileCollection = std::unordered_map<tile::Id, std::shared_ptr<VectorLayerTile>, tile::Id::Hasher>;

} // namespace nucleus::vector_tile
