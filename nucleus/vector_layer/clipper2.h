#include <vector>
#include <glm/glm.hpp>
#include <functional>

namespace nucleus::vector_layer {

// o1 model
//----------------------------------------------------
// Clip a polygon against a single infinite line.
// This helper is used in the Sutherland–Hodgman step.
//   - edgeFunc returns a positive value if a point
//     is inside the clipping boundary, and negative
//     (or zero) if it's outside.
//   - intersectFunc finds the intersection of a line
//     segment with the boundary line.
//----------------------------------------------------
static std::vector<glm::vec2> ClipAgainstEdge(
    const std::vector<glm::vec2> &inputPolygon,
    std::function<bool(const glm::vec2 &)> insideFunc,
    std::function<glm::vec2(const glm::vec2 &, const glm::vec2 &)> intersectFunc)
{
    std::vector<glm::vec2> outputPolygon;
    if (inputPolygon.empty())
        return outputPolygon;

    // Current point and the previous one
    glm::vec2 prev = inputPolygon.back();
    bool prevInside = insideFunc(prev);

    for (const auto &curr : inputPolygon)
    {
        bool currInside = insideFunc(curr);
        if (prevInside && currInside)
        {
            // Both endpoints are inside -> current vertex is added
            outputPolygon.push_back(curr);
        }
        else if (prevInside && !currInside)
        {
            // Previous is inside, current is outside -> clip at boundary
            outputPolygon.push_back(intersectFunc(prev, curr));
        }
        else if (!prevInside && currInside)
        {
            // Previous is outside, current is inside -> add intersection & current
            outputPolygon.push_back(intersectFunc(prev, curr));
            outputPolygon.push_back(curr);
        }
        // else if both outside, do nothing

        prev = curr;
        prevInside = currInside;
    }

    return outputPolygon;
}

//----------------------------------------------------
// Clipping a polygon with an axis-aligned rectangle
// using Sutherland–Hodgman. The rectangle is defined
// by [minX, maxX] x [minY, maxY].
//----------------------------------------------------
std::vector<glm::vec2> clip_polygon_with_cell2(const std::vector<glm::vec2>& polygon, const glm::uvec4& cell, const float scale)
{
    std::vector<glm::vec2> scaled;
    for (const glm::vec2& v : polygon) {
        scaled.push_back(v * scale);
    }

    // 1. Clip against left boundary: x = minX
    auto clippedLeft = ClipAgainstEdge(
        scaled,
        [cell](const glm::vec2& p) { return p.x >= cell.x; },
        [cell](const glm::vec2& p1, const glm::vec2& p2) {
            // Intersection with vertical line x = cell.x
            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            if (std::abs(dx) < 1e-8f)
            {
                // Edge case: vertical line - just return p1
                return p1;
            }
            float t = (cell.x - p1.x) / dx;
            return glm::vec2(cell.x, p1.y + t * dy);
        });

    // 2. Clip against right boundary: x = maxX
    auto clippedRight = ClipAgainstEdge(
        clippedLeft,
        [cell](const glm::vec2& p) { return p.x <= cell.z; },
        [cell](const glm::vec2& p1, const glm::vec2& p2) {
            // Intersection with vertical line x = cell.z
            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            if (std::abs(dx) < 1e-8f)
            {
                return p1;
            }
            float t = (cell.z - p1.x) / dx;
            return glm::vec2(cell.z, p1.y + t * dy);
        });

    // 3. Clip against bottom boundary: y = minY
    auto clippedBottom = ClipAgainstEdge(
        clippedRight,
        [cell](const glm::vec2& p) { return p.y >= cell.y; },
        [cell](const glm::vec2& p1, const glm::vec2& p2) {
            // Intersection with horizontal line y = cell.y
            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            if (std::abs(dy) < 1e-8f)
            {
                return p1;
            }
            float t = (cell.y - p1.y) / dy;
            return glm::vec2(p1.x + t * dx, cell.y);
        });

    // 4. Clip against top boundary: y = maxY
    auto clippedTop = ClipAgainstEdge(
        clippedBottom,
        [cell](const glm::vec2& p) { return p.y <= cell.w; },
        [cell](const glm::vec2& p1, const glm::vec2& p2) {
            // Intersection with horizontal line y = cell.w
            float dx = p2.x - p1.x;
            float dy = p2.y - p1.y;
            if (std::abs(dy) < 1e-8f)
            {
                return p1;
            }
            float t = (cell.w - p1.y) / dy;
            return glm::vec2(p1.x + t * dx, cell.w);
        });

    return clippedTop;
}

// int main()
// {
//     // Example polygons (using glm::vec2)
//     // Let's say we have two polygons:
//     std::vector<glm::vec2> polygon1 = {
//         glm::vec2(0.0f, 0.0f),
//         glm::vec2(5.0f, 0.0f),
//         glm::vec2(5.0f, 4.0f),
//         glm::vec2(3.0f, 5.0f),
//         glm::vec2(0.0f, 3.0f)
//     };

//     std::vector<glm::vec2> polygon2 = {
//         glm::vec2(2.0f, 1.0f),
//         glm::vec2(6.0f, 1.0f),
//         glm::vec2(6.0f, 5.0f),
//         glm::vec2(2.0f, 5.0f)
//     };

//     // Store them in a container
//     std::vector<std::vector<glm::vec2>> polygons = { polygon1, polygon2 };

//     // Define the grid
//     // For simplicity, let's say we have a 2x2 grid from (0,0) to (6,6).
//     // Each cell is 3x3 in size:
//     //  Cells:
//     //    (0,0)-(3,3), (3,0)-(6,3)
//     //    (0,3)-(3,6), (3,3)-(6,6)
//     float gridcell.min.x = 0.0f;
//     float gridMinY = 0.0f;
//     float gridMaxX = 6.0f;
//     float gridMaxY = 6.0f;
//     float cellSize = 3.0f;
//     int cellsX = 2; // how many cells horizontally
//     int cellsY = 2; // how many cells vertically

//     // For each cell, clip each polygon
//     for (int gy = 0; gy < cellsY; ++gy)
//     {
//         for (int gx = 0; gx < cellsX; ++gx)
//         {
//             float cellMinX = gridMinX + gx * cellSize;
//             float cellMaxX = cellMinX + cellSize;
//             float cellMinY = gridMinY + gy * cellSize;
//             float cellMaxY = cellMinY + cellSize;

//             std::cout << "Clipping against cell ["
//                       << cellMinX << ", " << cellMinY << "] - ["
//                       << cellMaxX << ", " << cellMaxY << "]\n";

//             for (size_t i = 0; i < polygons.size(); ++i)
//             {
//                 auto &poly = polygons[i];
//                 auto clipped = clip_polygon_with_cell(poly, cellMinX, cellMinY, cellMaxX, cellMaxY);

//                 std::cout << "  Polygon " << i << " clipped vertices:\n";
//                 if (clipped.empty())
//                 {
//                     std::cout << "    - No intersection with this cell.\n";
//                 }
//                 else
//                 {
//                     for (auto &v : clipped)
//                     {
//                         std::cout << "    (" << v.x << ", " << v.y << ")\n";
//                     }
//                 }
//             }
//         }
//     }

// return 0;
// }
} // namespace nucleus::vector_layer
