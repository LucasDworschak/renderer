/*****************************************************************************
 * AlpineMaps.org
 * Copyright (C) 2024 Adam Celarek
 * Copyright (C) 2024 Gerald Kimmersdorfer
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

#include "ColourTexture.h"

#include <QtGlobal>
#include <array>
#include <cstdint>
#include <stdexcept>

#define GOOFYTC_IMPLEMENTATION
#include <GoofyTC/goofy_tc.h>


namespace {

std::vector<uint8_t> to_dxt1(const nucleus::Raster<glm::u8vec4>& image)
{
    assert(image.width() == image.height());
    assert(image.width() % 16 == 0);
    assert(image.size_per_line() * image.height() == image.width() * image.height() * 4);
    assert(image.size_in_bytes() == image.width() * image.height() * 4);

    struct alignas(16) AlignedBlock
    {
        std::array<uint8_t, 16> data;
    };
    static_assert(sizeof(AlignedBlock) == 16);

    const auto n_bytes_in = size_t(image.size_in_bytes());
    const auto n_bytes_out = image.width() * image.height() / 2;
    assert(n_bytes_in % sizeof(AlignedBlock) == 0);

    auto aligned_in = std::vector<AlignedBlock>(n_bytes_in / sizeof(AlignedBlock));
    auto data_ptr = reinterpret_cast<uint8_t*>(aligned_in.data());
    std::copy(image.bytes(), image.bytes() + n_bytes_in, data_ptr);

    std::vector<uint8_t> compressed(n_bytes_out);
    const auto result = goofy::compressDXT1(compressed.data(), data_ptr, (uint32_t)image.width(), (uint32_t)image.height(), (uint32_t)image.width() * 4);
    assert(result == 0);
    Q_UNUSED(result);

    return compressed;
}

std::vector<uint8_t> to_etc1(const nucleus::Raster<glm::u8vec4>& image)
{
    assert(image.width() == image.height());
    assert(image.width() % 16 == 0);
    assert(image.size_per_line() * image.height() == image.width() * image.height() * 4);
    assert(image.size_in_bytes() == image.width() * image.height() * 4);

    struct alignas(16) AlignedBlock
    {
        std::array<uint8_t, 16> data;
    };
    static_assert(sizeof(AlignedBlock) == 16);

    const auto n_bytes_in = size_t(image.size_in_bytes());
    const auto n_bytes_out = image.width() * image.height() / 2;
    assert(n_bytes_in % sizeof(AlignedBlock) == 0);

    auto aligned_in = std::vector<AlignedBlock>(n_bytes_in / sizeof(AlignedBlock));
    auto data_ptr = reinterpret_cast<uint8_t*>(aligned_in.data());
    std::copy(image.bytes(), image.bytes() + n_bytes_in, data_ptr);

    std::vector<uint8_t> compressed(n_bytes_out);
    const auto result = goofy::compressETC1(compressed.data(), data_ptr, (uint32_t)image.width(), (uint32_t)image.height(), (uint32_t)image.width() * 4);
    assert(result == 0);
    Q_UNUSED(result);

    return compressed;
}

std::vector<uint8_t> to_uncompressed_rgba(const nucleus::Raster<glm::u8vec4>& image)
{
    size_t size_in_bytes = image.size_in_bytes();
    assert(size_in_bytes == (size_t)image.width() * image.height() * 4);
    std::vector<uint8_t> data(size_in_bytes);
    // Note: AND Another copy... We copy the data two times (once in stb_image_loader.cpp)
    std::copy(image.bytes(), image.bytes() + image.size_in_bytes(), data.data());
    return data;
}

std::vector<uint8_t> to_compressed(const nucleus::Raster<glm::u8vec4>& image, nucleus::utils::ColourTexture::Format algorithm)
{
    using Algorithm = nucleus::utils::ColourTexture::Format;
    assert(image.width() == image.height());
    assert(image.width() % 4 == 0 || image.width() == 2 || image.width() == 1);

    switch (algorithm) {
    case Algorithm::Uncompressed_RGBA:
        return to_uncompressed_rgba(image);
    case nucleus::utils::ColourTexture::Format::DXT1: {
        if (image.width() >= 16)
            return to_dxt1(image);
        glm::u32vec4 avg = {};
        for (const auto& c : image.buffer()) {
            avg += glm::u32vec4(c);
        }
        avg /= image.buffer_length();
        auto v = to_dxt1(nucleus::resize(image, { 16, 16 }, glm::u8vec4(avg)));
        if (image.width() == 8) {
            for (auto i = 16u; i < 32; ++i)
                v[i] = v[i + 16u];
            v.resize(32);
        } else if (image.width() <= 4) {
            v.resize(8);
        }
        return v;
    }
    case nucleus::utils::ColourTexture::Format::ETC1: {
        if (image.width() >= 16)
            return to_etc1(image);
        glm::u32vec4 avg = {};
        for (const auto& c : image.buffer()) {
            avg += glm::u32vec4(c);
        }
        avg /= image.buffer_length();
        auto v = to_etc1(nucleus::resize(image, { 16, 16 }, glm::u8vec4(avg)));
        if (image.width() == 8) {
            for (auto i = 16u; i < 32; ++i)
                v[i] = v[i + 16u];
            v.resize(32);
        } else if (image.width() <= 4) {
            v.resize(8);
        }
        return v;
    }
    }
    throw std::runtime_error("Unsupported algorithm for nucleus::Raster<glm::u8vec4>");
}
} // namespace

nucleus::utils::ColourTexture::ColourTexture(const nucleus::Raster<glm::u8vec4>& image, Format format)
    : m_data(to_compressed(image, format))
    , m_width(unsigned(image.width()))
    , m_height(unsigned(image.height()))
    , m_format(format)
{
}

nucleus::utils::MipmappedColourTexture nucleus::utils::generate_mipmapped_colour_texture(
    const nucleus::Raster<glm::u8vec4>& texture, ColourTexture::Format format)
{
    auto mip_levels = nucleus::generate_mipmap(texture);
    nucleus::utils::MipmappedColourTexture colour_texture = {};
    for (const auto& level : mip_levels) {
        colour_texture.emplace_back(level, format);
    }
    return colour_texture;
}
