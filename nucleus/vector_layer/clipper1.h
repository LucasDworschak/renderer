#include <vector>
#include <glm/glm.hpp>

namespace nucleus::vector_layer {

using glm::vec2;
using Polygon = std::vector<vec2>;
using namespace std;

// Check if point is inside clip edge
bool inside(const vec2& p, int edge, float value) {
    switch (edge) {
        case 0: return p.x >= value; // Left
        case 1: return p.x <= value; // Right
        case 2: return p.y >= value; // Bottom
        case 3: return p.y <= value; // Top
        default: return false;
    }
}

// Compute intersection point of edge with clip boundary
vec2 intersect(const vec2& a, const vec2& b, int edge, float value) {
    float t;
    switch (edge) {
        case 0: // Left
            t = (value - a.x) / (b.x - a.x);
            return vec2(value, a.y + t * (b.y - a.y));
        case 1: // Right
            t = (value - a.x) / (b.x - a.x);
            return vec2(value, a.y + t * (b.y - a.y));
        case 2: // Bottom
            t = (value - a.y) / (b.y - a.y);
            return vec2(a.x + t * (b.x - a.x), value);
        case 3: // Top
            t = (value - a.y) / (b.y - a.y);
            return vec2(a.x + t * (b.x - a.x), value);
        default:
            return vec2(0);
    }
}

// Sutherland-Hodgman polygon clipping with a rectangular AABB
Polygon clip_polygon_with_cell(const Polygon& polygon, const glm::uvec4& cell, const float scale)
{
    Polygon output;
    for (const vec2& v : polygon) {
        output.push_back(v * scale);
    }

    auto cell_min = glm::vec2(cell.x, cell.y);
    auto cell_max = glm::vec2(cell.z, cell.w);

    for (int edge = 0; edge < 4; ++edge) {
        Polygon input = output;
        output.clear();
        if (input.empty()) break;

        vec2 S = input.back();
        for (const vec2& E : input) {
            if (inside(E, edge, (edge % 2 == 0) ? cell_min[edge / 2] : cell_max[edge / 2])) {
                if (!inside(S, edge, (edge % 2 == 0) ? cell_min[edge / 2] : cell_max[edge / 2])) {
                    output.push_back(intersect(S, E, edge, (edge % 2 == 0) ? cell_min[edge / 2] : cell_max[edge / 2]));
                }
                output.push_back(E);
            } else if (inside(S, edge, (edge % 2 == 0) ? cell_min[edge / 2] : cell_max[edge / 2])) {
                output.push_back(intersect(S, E, edge, (edge % 2 == 0) ? cell_min[edge / 2] : cell_max[edge / 2]));
            }
            S = E;
        }
    }
    return output;
}

} // namespace nucleus::vector_layer
