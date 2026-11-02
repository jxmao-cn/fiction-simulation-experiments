#pragma once
// =============================================================================
// io/video_encoder.h — ffmpeg 视频编码封装
// =============================================================================
// 用途：把 <framesDir>/frame_%06d.<ext> 图像序列编码为 mp4。
//   - 优先 libx264（yuv420p + faststart，兼容常见播放器/浏览器）
//   - 失败回退 mpeg4
//   - 无 ffmpeg 时返回 false，调用方可保留帧目录并提示手工命令
// =============================================================================

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

namespace fiction {

/// PATH 中是否存在 ffmpeg
inline bool ffmpegAvailable() {
    return std::system("command -v ffmpeg > /dev/null 2>&1") == 0;
}

/// ffmpeg image2 输入模式（%06d 对应 frame_000000.ext）
inline std::string ffmpegFramePattern(const std::string& framesDir, const std::string& ext) {
    return framesDir + "/frame_%06d." + ext;
}

namespace detail {

inline std::string ffmpegCommand(const std::string& pattern, double fps,
                                 const std::string& codecArgs, const std::string& outPath) {
    std::ostringstream os;
    os << "ffmpeg -y -loglevel error -framerate " << fps
       << " -i \"" << pattern << "\" " << codecArgs << " \"" << outPath << "\"";
    return os.str();
}

} // namespace detail

/// 编码图像序列为 mp4；成功返回 true。logCommand=true 时打印实际命令
inline bool encodeVideo(const std::string& pattern, double fps,
                        const std::string& outPath, bool logCommand = false) {
    if (!ffmpegAvailable()) return false;
    const std::string x264 = detail::ffmpegCommand(
        pattern, fps, "-c:v libx264 -pix_fmt yuv420p -crf 20 -movflags +faststart", outPath);
    if (logCommand) std::cout << "  [ffmpeg] " << x264 << "\n";
    if (std::system(x264.c_str()) == 0) return true;

    const std::string mpeg4 =
        detail::ffmpegCommand(pattern, fps, "-c:v mpeg4 -q:v 3", outPath);
    if (logCommand) std::cout << "  [ffmpeg] " << mpeg4 << "\n";
    return std::system(mpeg4.c_str()) == 0;
}

} // namespace fiction
