/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2023 Adam Celarek
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

#include <array>
#include <cmath>
#include <cstdint>

#include <glm/glm.hpp>

namespace nucleus::utils::bit_coding {
inline glm::vec2 to_f16f16(const glm::u8vec4& v)
{
    const auto dec = [](uint8_t v1, uint8_t v2) {
        return float((v1 << 8) | v2) / 65535.f;
    };
    return { dec(v[0], v[1]), dec(v[2], v[3]) };
};

inline const glm::vec4 u32_to_f8_4(uint32_t number)
{
    return { float((number >> 24) & 255) / 255.0f, float((number >> 16) & 255) / 255.0f, float((number >> 8) & 255) / 255.0f, float(number & 255) / 255.0f };
}

inline uint32_t f8_4_to_u32(const glm::vec4& v)
{
    // rounding before casting to int to prevent rounding down 0.999 values that might happen due to float errors
    return (int(round(v[0] * 255.0f)) << 24) | (int(round(v[1] * 255.0f)) << 16) | (int(round(v[2] * 255.0f)) << 8) | int(round(v[3] * 255.0f));
}

inline glm::vec<2, uint32_t> u16_u24_u24_to_u32_2(const uint16_t u16Bit, const uint32_t u24Bit1, const uint32_t u24Bit2)
{
    uint32_t a = u16Bit << 16;

    a = a | (u24Bit1 >> 8);
    uint32_t b = u24Bit1 << 24;

    b = b | (u24Bit2 & ((1u << 24) - 1));

    return { a, b };
}
inline glm::vec<3, uint32_t> u32_2_to_u16_u24_u24(glm::vec<2, uint32_t> packed)
{
    glm::vec<3, uint32_t> out;

    out.x = packed.x >> 16;

    out.y = (packed.x & ((1u << 16) - 1)) << 8;
    out.y |= packed.y >> 24;

    out.z = packed.y & ((1u << 24) - 1);

    return out;
}

inline uint32_t u24_u8_to_u32(uint32_t a, uint8_t b)
{
    uint32_t out = a << 8;

    out |= b;

    return out;
}

inline glm::vec<2, uint32_t> u32_to_u24_u8(uint32_t packed)
{
    glm::vec<2, uint32_t> out;

    out.x = packed >> 8;
    out.y = packed & 255u;

    return out;
}

} // namespace nucleus::utils::bit_coding
