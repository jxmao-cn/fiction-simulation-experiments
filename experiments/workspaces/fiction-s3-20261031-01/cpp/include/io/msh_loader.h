#pragma once
// =============================================================================
// io/msh_loader.h — MshIO 4.1 ASCII 四面体网格读取
// =============================================================================
// 职责：
//   - 读取 IPC（ipc-sim/IPC）input/tetMeshes/*.msh 使用的 MshIO 4.1 ASCII 格式
//   - 解析 $Nodes / $Elements（tet4，elementType=4）/ $Surface（IPC 边界三角面段）
//   - 缺少 $Surface 时，从四面体面出现次数提取边界三角面（仅出现一次的面）
//   - 不做任何物理假设，只输出几何：节点、四面体、边界三角面
//
// 支持范围（已完成实测：IPC/input/tetMeshes 下 38 个文件全部为 4.1 ASCII）：
//   - $MeshFormat: "4.1 0 8"（二进制、非 4.1 版本直接报错）
//   - $Nodes: 标准 entity block（parametric 节点不支持）
//   - $Elements: 每个 block 一行一个单元，elementType=4 收集为四面体、
//                elementType=2 收集为表面三角面、其余类型跳过
//   - $Surface: IPC 自定义段（count + 每行 3 个节点 tag）；与显式三角面单元二选一
//   - 未知 $Section 自动跳过直到对应的 $EndSection
//
// 用法：
//   fiction::MshMesh mesh = fiction::loadMsh("cube.msh");
//   mesh.x          // [N,3] 参考位置
//   mesh.tets       // [T,4] 0-based 节点索引
//   mesh.surface    // [S,3] 0-based 边界三角面（可为空，见 fallback）
// =============================================================================

#include "types.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fiction {

// -----------------------------------------------------------------------------
// 输出数据结构
// -----------------------------------------------------------------------------
struct MshMesh {
    Positions x;                              ///< [N,3] 节点坐标（按 tag 出现顺序重排）
    std::vector<std::array<int, 4>> tets;     ///< [T,4] 四面体单元（0-based 索引）
    std::vector<std::array<int, 3>> surface;  ///< [S,3] 边界三角面（0-based 索引）
    int numTets() const { return static_cast<int>(tets.size()); }
    int numSurface() const { return static_cast<int>(surface.size()); }
    bool empty() const { return x.rows() == 0; }
};

// -----------------------------------------------------------------------------
// 内部：行读取器（统一处理 \r、行号、报错信息）
// -----------------------------------------------------------------------------
namespace detail {

class MshLineReader {
public:
    MshLineReader(std::istream& in, std::string source)
        : in_(in), source_(std::move(source)) {}

    /// 读取下一个非空行（去除行尾 \r 与首尾空白）；EOF 返回 false
    bool next(std::string& line) {
        while (std::getline(in_, line)) {
            ++lineNo_;
            while (!line.empty() && (line.back() == '\r' || line.back() == '\n'))
                line.pop_back();
            if (!line.empty()) return true;
        }
        return false;
    }

    /// 读取一个非空行，失败即报错
    std::string require(const std::string& what) {
        std::string line;
        if (!next(line)) throwError("文件提前结束（期望 " + what + "）");
        return line;
    }

    [[noreturn]] void throwError(const std::string& msg) const {
        throw std::runtime_error(source_ + ":" + std::to_string(lineNo_) + ": " + msg);
    }

    int lineNo() const { return lineNo_; }
    const std::string& source() const { return source_; }

private:
    std::istream& in_;
    std::string source_;
    int lineNo_ = 0;
};

/// 按空白切分一行
inline std::vector<std::string> tokenize(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> out;
    std::string tok;
    while (iss >> tok) out.push_back(tok);
    return out;
}

inline double toDouble(const std::string& tok, const MshLineReader& r) {
    try {
        size_t pos = 0;
        const double v = std::stod(tok, &pos);
        if (pos != tok.size()) r.throwError("无法解析浮点数: " + tok);
        return v;
    } catch (const std::exception&) {
        r.throwError("无法解析浮点数: " + tok);
    }
}

inline long toLong(const std::string& tok, const MshLineReader& r) {
    try {
        size_t pos = 0;
        const long v = std::stol(tok, &pos);
        if (pos != tok.size()) r.throwError("无法解析整数: " + tok);
        return v;
    } catch (const std::exception&) {
        r.throwError("无法解析整数: " + tok);
    }
}

/// 解析 $MeshFormat 段（调用时已读掉段名行）
inline void parseMeshFormat(MshLineReader& r) {
    const std::string line = r.require("MeshFormat 内容");
    const auto tok = tokenize(line);
    if (tok.size() < 3) r.throwError("MeshFormat 行应为 'version fileType dataSize'");
    if (tok[0] != "4.1")
        r.throwError("仅支持 MshIO 4.1，当前版本: " + tok[0]);
    if (toLong(tok[1], r) != 0)
        r.throwError("仅支持 ASCII 格式（fileType=0），当前为二进制");
    const std::string end = r.require("$EndMeshFormat");
    if (end != "$EndMeshFormat") r.throwError("期望 $EndMeshFormat，实际: " + end);
}

/// 解析 $Nodes 段
inline void parseNodes(MshLineReader& r, MshMesh& mesh,
                       std::unordered_map<long, int>& tagToIndex) {
    const auto head = tokenize(r.require("Nodes 头行"));
    if (head.size() < 4) r.throwError("Nodes 头行应为 'numEntityBlocks numNodes minTag maxTag'");
    const long numNodes = toLong(head[1], r);
    if (numNodes <= 0) r.throwError("节点数必须为正");

    std::vector<Vec3> coords(static_cast<size_t>(numNodes));
    int nextIndex = 0;

    const long numBlocks = toLong(head[0], r);
    for (long b = 0; b < numBlocks; ++b) {
        const auto bh = tokenize(r.require("Nodes block 头行"));
        if (bh.size() < 4)
            r.throwError("Nodes block 头行应为 'entityDim entityTag parametric numNodesInBlock'");
        const long parametric = toLong(bh[2], r);
        const long blockNodes = toLong(bh[3], r);
        if (parametric != 0)
            r.throwError("不支持 parametric 节点块");
        if (blockNodes <= 0) r.throwError("block 内节点数必须为正");

        // 先列出 block 内全部 tag，再列出对应坐标
        std::vector<long> tags(static_cast<size_t>(blockNodes));
        for (long i = 0; i < blockNodes; ++i) {
            const auto t = tokenize(r.require("节点 tag"));
            if (t.empty()) r.throwError("空节点 tag");
            tags[static_cast<size_t>(i)] = toLong(t[0], r);
        }
        for (long i = 0; i < blockNodes; ++i) {
            const auto t = tokenize(r.require("节点坐标"));
            if (t.size() < 3) r.throwError("节点坐标应为 3 个浮点数");
            if (nextIndex >= numNodes) r.throwError("节点数量超过头行声明");
            const long tag = tags[static_cast<size_t>(i)];
            if (tagToIndex.count(tag)) r.throwError("节点 tag 重复: " + std::to_string(tag));
            tagToIndex[tag] = nextIndex;
            coords[static_cast<size_t>(nextIndex)] =
                Vec3(toDouble(t[0], r), toDouble(t[1], r), toDouble(t[2], r));
            ++nextIndex;
        }
    }
    if (nextIndex != numNodes)
        r.throwError("实际节点数与头行声明不一致（声明 " + std::to_string(numNodes)
                     + "，读取 " + std::to_string(nextIndex) + "）");

    const std::string end = r.require("$EndNodes");
    if (end != "$EndNodes") r.throwError("期望 $EndNodes，实际: " + end);

    mesh.x.resize(numNodes, kDim);
    for (long i = 0; i < numNodes; ++i)
        mesh.x.row(static_cast<int>(i)) = coords[static_cast<size_t>(i)].transpose();
}

/// 解析 $Elements 段；type=4 收四面体，type=2 收表面三角面，其余跳过
inline void parseElements(MshLineReader& r,
                          std::vector<std::array<int, 3>>& triElements,
                          MshMesh& mesh,
                          const std::unordered_map<long, int>& tagToIndex) {
    const auto head = tokenize(r.require("Elements 头行"));
    if (head.size() < 4) r.throwError("Elements 头行应为 'numEntityBlocks numElements minTag maxTag'");

    auto resolve = [&](long tag) -> int {
        const auto it = tagToIndex.find(tag);
        if (it == tagToIndex.end())
            r.throwError("单元引用了不存在的节点 tag: " + std::to_string(tag));
        return it->second;
    };

    const long numBlocks = toLong(head[0], r);
    for (long b = 0; b < numBlocks; ++b) {
        const auto bh = tokenize(r.require("Elements block 头行"));
        if (bh.size() < 4)
            r.throwError("Elements block 头行应为 'entityDim entityTag elementType numElementsInBlock'");
        const long elemType = toLong(bh[2], r);
        const long blockElems = toLong(bh[3], r);
        if (blockElems <= 0) r.throwError("block 内单元数必须为正");

        for (long e = 0; e < blockElems; ++e) {
            const auto t = tokenize(r.require("单元行"));
            if (t.size() < 2) r.throwError("单元行缺少节点 tag");
            if (elemType == 4) {
                if (t.size() < 5) r.throwError("四面体行应为 'tag n0 n1 n2 n3'");
                mesh.tets.push_back({resolve(toLong(t[1], r)), resolve(toLong(t[2], r)),
                                     resolve(toLong(t[3], r)), resolve(toLong(t[4], r))});
            } else if (elemType == 2) {
                if (t.size() < 4) r.throwError("三角面行应为 'tag n0 n1 n2'");
                triElements.push_back({resolve(toLong(t[1], r)), resolve(toLong(t[2], r)),
                                       resolve(toLong(t[3], r))});
            }
            // 其余 elementType（1/3/5/8/...）本工具不需要，跳过该行
        }
    }

    const std::string end = r.require("$EndElements");
    if (end != "$EndElements") r.throwError("期望 $EndElements，实际: " + end);
}

/// 解析 $Surface 段（IPC 自定义：count + 每行 3 个节点 tag）
inline void parseSurface(MshLineReader& r, std::vector<std::array<int, 3>>& surface,
                         const std::unordered_map<long, int>& tagToIndex) {
    const auto head = tokenize(r.require("Surface 头行"));
    if (head.empty()) r.throwError("Surface 头行缺少三角面数量");
    const long count = toLong(head[0], r);
    if (count < 0) r.throwError("Surface 三角面数量不能为负");

    auto resolve = [&](long tag) -> int {
        const auto it = tagToIndex.find(tag);
        if (it == tagToIndex.end())
            r.throwError("Surface 引用了不存在的节点 tag: " + std::to_string(tag));
        return it->second;
    };

    for (long i = 0; i < count; ++i) {
        const auto t = tokenize(r.require("Surface 三角面行"));
        if (t.size() < 3) r.throwError("Surface 三角面行应为 3 个节点 tag");
        surface.push_back({resolve(toLong(t[0], r)), resolve(toLong(t[1], r)),
                           resolve(toLong(t[2], r))});
    }

    const std::string end = r.require("$EndSurface");
    if (end != "$EndSurface") r.throwError("期望 $EndSurface，实际: " + end);
}

/// 跳过未知段：一直读到匹配的 $End<name>
inline void skipSection(MshLineReader& r, const std::string& section) {
    const std::string endTag = "$End" + section.substr(1);
    std::string line;
    while (r.next(line)) {
        if (line == endTag) return;
    }
    r.throwError("未知段 " + section + " 缺少 " + endTag);
}

/// 从四面体提取边界三角面（面在两个四面体间共享，只出现一次的面即为边界）
inline std::vector<std::array<int, 3>> extractBoundaryFaces(
    const std::vector<std::array<int, 4>>& tets) {
    struct FaceInfo {
        int count = 0;
        std::array<int, 3> oriented{};
    };
    std::unordered_map<std::uint64_t, FaceInfo> faces;
    faces.reserve(tets.size() * 4);

    auto keyOf = [](int a, int b, int c) {
        std::array<int, 3> s{a, b, c};
        std::sort(s.begin(), s.end());
        return (static_cast<std::uint64_t>(s[0]) << 40) |
               (static_cast<std::uint64_t>(s[1]) << 20) |
               static_cast<std::uint64_t>(s[2]);
    };

    // 每个四面体的 4 个面（保持单元内定向）
    static constexpr int kFaceVerts[4][3] = {
        {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
    for (const auto& tet : tets) {
        for (const auto& f : kFaceVerts) {
            const int a = tet[f[0]], b = tet[f[1]], c = tet[f[2]];
            FaceInfo& info = faces[keyOf(a, b, c)];
            ++info.count;
            if (info.count == 1) info.oriented = {a, b, c};
        }
    }

    std::vector<std::array<int, 3>> boundary;
    boundary.reserve(faces.size() / 2 + 8);
    for (const auto& [key, info] : faces) {
        (void)key;
        if (info.count == 1) boundary.push_back(info.oriented);
    }
    // 排序保证不同运行/编译器下表面面片顺序稳定（OBJ 导出与测试可复现）
    std::sort(boundary.begin(), boundary.end());
    return boundary;
}

/// 真正的解析入口（流版本）
inline MshMesh parseMsh(std::istream& in, const std::string& source) {
    MshLineReader r(in, source);
    MshMesh mesh;
    std::unordered_map<long, int> tagToIndex;
    std::vector<std::array<int, 3>> triElements;
    std::vector<std::array<int, 3>> surfaceSection;
    bool sawNodes = false;
    bool sawElements = false;

    std::string line;
    while (r.next(line)) {
        if (line == "$MeshFormat") {
            parseMeshFormat(r);
        } else if (line == "$Nodes") {
            parseNodes(r, mesh, tagToIndex);
            sawNodes = true;
        } else if (line == "$Elements") {
            parseElements(r, triElements, mesh, tagToIndex);
            sawElements = true;
        } else if (line == "$Surface") {
            parseSurface(r, surfaceSection, tagToIndex);
        } else if (!line.empty() && line[0] == '$') {
            skipSection(r, line);
        } else {
            r.throwError("段外出现无法识别的内容: " + line);
        }
    }

    if (!sawNodes) throw std::runtime_error(source + ": 缺少 $Nodes 段");
    if (!sawElements) throw std::runtime_error(source + ": 缺少 $Elements 段");
    if (mesh.tets.empty())
        throw std::runtime_error(source + ": 未找到 elementType=4（tet4）四面体单元");

    // 表面来源优先级：显式三角面单元 > $Surface 段 > 边界提取
    if (!triElements.empty()) {
        mesh.surface = std::move(triElements);
    } else if (!surfaceSection.empty()) {
        mesh.surface = std::move(surfaceSection);
    } else {
        mesh.surface = extractBoundaryFaces(mesh.tets);
    }
    return mesh;
}

} // namespace detail

// -----------------------------------------------------------------------------
// 公共接口
// -----------------------------------------------------------------------------
/// 从文件读取 .msh（失败抛 std::runtime_error，信息含文件名与行号）
inline MshMesh loadMsh(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("无法打开 msh 文件: " + path);
    return detail::parseMsh(in, path);
}

/// 从流读取 .msh（便于测试；source 仅用于报错信息）
inline MshMesh loadMsh(std::istream& in, const std::string& source = "<stream>") {
    return detail::parseMsh(in, source);
}

} // namespace fiction
