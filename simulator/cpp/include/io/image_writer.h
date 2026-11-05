#pragma once
// =============================================================================
// io/image_writer.h — 图像写出（PNG，启用 zlib 时；否则 PPM 回退）
// =============================================================================
// 用途：CPU 软光栅化结果的落盘，供 ffmpeg 编码视频或直接查看。
//   - writePng: 8-bit truecolor PNG（IHDR/IDAT/IEND，zlib 压缩），需 FICTION_HAS_ZLIB
//   - writePpm: 二进制 P6，无任何依赖（ffmpeg 可直接读取）
//   - writeImage: 按扩展名分派；未启用 zlib 时 PNG 请求回退为 PPM
// =============================================================================

#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#ifdef FICTION_HAS_ZLIB
#include <zlib.h>
#endif

namespace fiction {

/// 二进制 PPM（P6，RGB8）
inline bool writePpm(const std::string& path, int w, int h,
                     const std::vector<unsigned char>& rgb) {
    if (w <= 0 || h <= 0 || rgb.size() != static_cast<size_t>(w) * h * 3) return false;
    std::ofstream out(path, std::ios::binary);
    if (!out) return false;
    out << "P6\n" << w << " " << h << "\n255\n";
    out.write(reinterpret_cast<const char*>(rgb.data()),
              static_cast<std::streamsize>(rgb.size()));
    out.close(); // 关闭时刷新缓冲区，捕获磁盘满等延迟写入错误。
    return static_cast<bool>(out);
}

#ifdef FICTION_HAS_ZLIB
/// 8-bit truecolor PNG（无隔行）
inline bool writePng(const std::string& path, int w, int h,
                     const std::vector<unsigned char>& rgb) {
    if (w <= 0 || h <= 0 || rgb.size() != static_cast<size_t>(w) * h * 3) return false;

    // 每行前置 filter byte 0
    const size_t rowBytes = static_cast<size_t>(w) * 3;
    std::vector<unsigned char> raw((rowBytes + 1) * static_cast<size_t>(h));
    for (int y = 0; y < h; ++y) {
        raw[(rowBytes + 1) * y] = 0;
        std::memcpy(&raw[(rowBytes + 1) * y + 1], &rgb[rowBytes * y], rowBytes);
    }

    uLongf compSize = compressBound(static_cast<uLong>(raw.size()));
    std::vector<unsigned char> comp(compSize);
    if (compress2(comp.data(), &compSize, raw.data(),
                  static_cast<uLong>(raw.size()), 6) != Z_OK)
        return false;
    comp.resize(compSize);

    std::ofstream out(path, std::ios::binary);
    if (!out) return false;

    auto writeU32 = [&out](std::uint32_t v) {
        unsigned char b[4] = {static_cast<unsigned char>(v >> 24),
                              static_cast<unsigned char>(v >> 16),
                              static_cast<unsigned char>(v >> 8),
                              static_cast<unsigned char>(v)};
        out.write(reinterpret_cast<const char*>(b), 4);
    };
    auto writeChunk = [&](const char type[4], const unsigned char* data, size_t len) {
        writeU32(static_cast<std::uint32_t>(len));
        out.write(type, 4);
        if (len) out.write(reinterpret_cast<const char*>(data),
                           static_cast<std::streamsize>(len));
        uLong crc = crc32(0L, Z_NULL, 0);
        crc = crc32(crc, reinterpret_cast<const Bytef*>(type), 4);
        if (len) crc = crc32(crc, data, static_cast<uInt>(len));
        writeU32(static_cast<std::uint32_t>(crc));
    };

    static const unsigned char kSignature[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    out.write(reinterpret_cast<const char*>(kSignature), 8);

    unsigned char ihdr[13];
    ihdr[0] = static_cast<unsigned char>(w >> 24);
    ihdr[1] = static_cast<unsigned char>(w >> 16);
    ihdr[2] = static_cast<unsigned char>(w >> 8);
    ihdr[3] = static_cast<unsigned char>(w);
    ihdr[4] = static_cast<unsigned char>(h >> 24);
    ihdr[5] = static_cast<unsigned char>(h >> 16);
    ihdr[6] = static_cast<unsigned char>(h >> 8);
    ihdr[7] = static_cast<unsigned char>(h);
    ihdr[8] = 8;   // bit depth
    ihdr[9] = 2;   // color type: truecolor RGB
    ihdr[10] = 0;  // compression
    ihdr[11] = 0;  // filter
    ihdr[12] = 0;  // interlace
    writeChunk("IHDR", ihdr, sizeof(ihdr));
    writeChunk("IDAT", comp.data(), comp.size());
    writeChunk("IEND", nullptr, 0);
    out.close();
    return static_cast<bool>(out);
}
#endif // FICTION_HAS_ZLIB

/// 按扩展名分派（.ppm 或 .png；未启用 zlib 时 png 回退为 ppm 字节流）
inline bool writeImage(const std::string& path, int w, int h,
                       const std::vector<unsigned char>& rgb) {
    if (path.size() >= 4 && path.compare(path.size() - 4, 4, ".ppm") == 0)
        return writePpm(path, w, h, rgb);
#ifdef FICTION_HAS_ZLIB
    return writePng(path, w, h, rgb);
#else
    return writePpm(path, w, h, rgb);
#endif
}

/// 当前构建下 writeImage 实际使用的格式扩展名（视频帧命名用）
inline const char* imageExtension() {
#ifdef FICTION_HAS_ZLIB
    return "png";
#else
    return "ppm";
#endif
}

} // namespace fiction
