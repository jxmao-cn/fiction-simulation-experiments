#pragma once
// =============================================================================
// io/obj_writer.h — 表面三角面 OBJ 帧导出
// =============================================================================
// 用途：把当前节点位置 + 边界三角面写成 Wavefront OBJ（v/f，1-based），
//       命名与内容对齐 IPC 的输出习惯，便于 MeshLab / Blender 对照查看。
// 说明：OBJ 只含几何，不写法向/材质；失败返回 false（由调用方决定是否报错）。
// =============================================================================

#include "types.h"

#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace fiction {

/// 写出 [N,3] 位置与 [F,3] 三角面的 OBJ 文件
inline bool writeSurfaceObj(const std::string& path,
                            const Positions& x,
                            const std::vector<std::array<int, 3>>& faces) {
    std::ofstream out(path);
    if (!out) return false;

    out << "# fiction mass-spring surface export\n";
    out << std::setprecision(9);
    for (int i = 0; i < x.rows(); ++i) {
        out << "v " << x(i, 0) << ' ' << x(i, 1) << ' ' << x(i, 2) << '\n';
    }
    for (const auto& f : faces) {
        out << "f " << (f[0] + 1) << ' ' << (f[1] + 1) << ' ' << (f[2] + 1) << '\n';
    }
    out.close(); // 将缓冲区写入错误也纳入返回值，而非等析构时静默丢失。
    return static_cast<bool>(out);
}

} // namespace fiction
