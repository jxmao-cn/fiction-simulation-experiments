#pragma once
// =============================================================================
// scene/multibody_scene.h — S6a 无接触多体全局场景装配
// =============================================================================
// 职责：
//   - 把同一份四面体资产按“纯平移”复制成 nb 个体，拼成一套全局状态：
//     全局节点、全局 tet、全局集中质量、全局 DBC 索引均按 body_id 连续排列
//   - 为每个体记录 body_id / 平移 / 节点区间 / tet 区间 / 表面区间 / DBC 节点数 / 质量
//   - 提供只读资产体检（bbox、体积、质量、正定向、表面闭合）
//
// 约束（S6a 预登记）：
//   - 首轮只平移：rotation 固定单位阵、scale 固定 1；不实现接触/堆叠
//   - 各体之间只共享“同一个时间步入口”，不共享 tet / 节点 / 约束；
//     因此全局 Hessian 应当块对角（按 body_id 检查，期望跨体条目为 0）
// =============================================================================

#include "io/msh_loader.h"
#include "scene/mass_spring_builder.h"
#include "types.h"

#include <array>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace fiction {

// -----------------------------------------------------------------------------
// 只读资产体检
// -----------------------------------------------------------------------------
struct MeshAssetAudit {
    int numNodes = 0;
    int numTets = 0;
    int numSurfaceTriangles = 0;
    std::array<double, 6> bbox{};
    double volume = 0.0;          ///< sum |det Dm| / 6
    double mass = 0.0;            ///< density * volume
    double minTetVolume = 0.0;
    double maxTetVolume = 0.0;
    double minSignedDet = 0.0;    ///< det([x1-x0, x2-x0, x3-x0]) 的最小/最大值
    double maxSignedDet = 0.0;
    int nonPositiveTets = 0;      ///< signed det <= 0 的 tet 数
    int degenerateTets = 0;       ///< |det| < 1e-18 的 tet 数
    bool positiveOrientation = false;
    int uniqueSurfaceNodes = 0;
    int uniqueSurfaceEdges = 0;
    int surfaceEulerCharacteristic = 0;   ///< V - E + F（闭合亏格 0 应为 2）
    bool surfaceClosedManifold = false;   ///< 每条表面边恰好被 2 个三角面共享
    int boundaryFacesFromTets = 0;        ///< 只出现一次的 tet 面数
    bool surfaceMatchesTetBoundary = false;
    double minSurfaceArea = 0.0;
    double maxSurfaceArea = 0.0;
};

/// 只读体检：不改动 mesh，不写文件
inline MeshAssetAudit auditMesh(const MshMesh& mesh, double density) {
    if (mesh.empty() || mesh.tets.empty())
        throw std::runtime_error("auditMesh: 空网格或没有四面体");
    MeshAssetAudit a;
    a.numNodes = static_cast<int>(mesh.x.rows());
    a.numTets = mesh.numTets();
    a.numSurfaceTriangles = mesh.numSurface();
    a.bbox = boundingBox(mesh.x);
    a.minTetVolume = 1e300;
    a.maxTetVolume = -1e300;
    a.minSignedDet = 1e300;
    a.maxSignedDet = -1e300;

    for (const auto& tet : mesh.tets) {
        const Vec3 p0 = mesh.x.row(tet[0]).transpose();
        Eigen::Matrix3d Dm;
        Dm.col(0) = mesh.x.row(tet[1]).transpose() - p0;
        Dm.col(1) = mesh.x.row(tet[2]).transpose() - p0;
        Dm.col(2) = mesh.x.row(tet[3]).transpose() - p0;
        const double det = Dm.determinant();
        const double vol = std::abs(det) / 6.0;
        a.volume += vol;
        a.minTetVolume = std::min(a.minTetVolume, vol);
        a.maxTetVolume = std::max(a.maxTetVolume, vol);
        a.minSignedDet = std::min(a.minSignedDet, det);
        a.maxSignedDet = std::max(a.maxSignedDet, det);
        if (!(det > 0.0)) ++a.nonPositiveTets;
        if (std::abs(det) < 1e-18) ++a.degenerateTets;
    }
    a.mass = density * a.volume;
    a.positiveOrientation = (a.nonPositiveTets == 0);

    // --- 表面闭合：每条无向边恰好被 2 个三角面共享 ---
    auto edgeKey = [](int u, int v) -> std::uint64_t {
        if (u > v) std::swap(u, v);
        return (static_cast<std::uint64_t>(u) << 32) | static_cast<std::uint32_t>(v);
    };
    auto triKey = [](std::array<int, 3> t) -> std::uint64_t {
        std::sort(t.begin(), t.end());
        return (static_cast<std::uint64_t>(t[0]) << 40) |
               (static_cast<std::uint64_t>(t[1]) << 20) |
               static_cast<std::uint64_t>(t[2]);
    };

    std::unordered_map<std::uint64_t, int> edgeUse;
    std::unordered_set<int> surfaceNodes;
    a.minSurfaceArea = 1e300;
    a.maxSurfaceArea = -1e300;
    for (const auto& f : mesh.surface) {
        for (int k = 0; k < 3; ++k) {
            const int u = f[k], v = f[(k + 1) % 3];
            ++edgeUse[edgeKey(u, v)];
            surfaceNodes.insert(u);
        }
        const Vec3 e1 = (mesh.x.row(f[1]) - mesh.x.row(f[0])).transpose();
        const Vec3 e2 = (mesh.x.row(f[2]) - mesh.x.row(f[0])).transpose();
        const double area = 0.5 * e1.cross(e2).norm();
        a.minSurfaceArea = std::min(a.minSurfaceArea, area);
        a.maxSurfaceArea = std::max(a.maxSurfaceArea, area);
    }
    if (!mesh.surface.empty()) {
        bool closed = true;
        for (const auto& [k, c] : edgeUse) {
            (void)k;
            if (c != 2) closed = false;
        }
        a.surfaceClosedManifold = closed;
        a.uniqueSurfaceNodes = static_cast<int>(surfaceNodes.size());
        a.uniqueSurfaceEdges = static_cast<int>(edgeUse.size());
        a.surfaceEulerCharacteristic =
            a.uniqueSurfaceNodes - a.uniqueSurfaceEdges + a.numSurfaceTriangles;

        // tet 面出现次数：1 次即边界（与文件里的 $Surface 比对）
        std::unordered_map<std::uint64_t, int> faceCount;
        static constexpr int kFaceVerts[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
        for (const auto& tet : mesh.tets) {
            for (const auto& fv : kFaceVerts)
                ++faceCount[triKey({tet[fv[0]], tet[fv[1]], tet[fv[2]]})];
        }
        std::unordered_set<std::uint64_t> boundary;
        for (const auto& [k, c] : faceCount)
            if (c == 1) boundary.insert(k);
        a.boundaryFacesFromTets = static_cast<int>(boundary.size());
        bool match = (a.boundaryFacesFromTets == a.numSurfaceTriangles);
        if (match) {
            for (const auto& f : mesh.surface)
                if (!boundary.count(triKey(f))) { match = false; break; }
        }
        a.surfaceMatchesTetBoundary = match;
    }
    return a;
}

// -----------------------------------------------------------------------------
// 文件哈希（FNV-1a 64，仅用于把资产和结果绑定；非密码学）
// -----------------------------------------------------------------------------
inline std::string fnv1aFileHash(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("无法打开文件做哈希: " + path);
    std::uint64_t h = 1469598103934665603ULL;
    char buf[1 << 16];
    while (in.read(buf, sizeof(buf)) || in.gcount() > 0) {
        const std::streamsize n = in.gcount();
        for (std::streamsize i = 0; i < n; ++i) {
            h ^= static_cast<unsigned char>(buf[i]);
            h *= 1099511628211ULL;
        }
    }
    char hex[32];
    std::snprintf(hex, sizeof(hex), "%016llx", static_cast<unsigned long long>(h));
    return std::string(hex);
}

// -----------------------------------------------------------------------------
// 多体场景
// -----------------------------------------------------------------------------
struct BodyTransform {
    Vec3 translation = Vec3::Zero();  ///< 纯平移（S6a 首轮不做旋转/缩放）
};

struct BodyRecord {
    int bodyId = 0;
    std::string assetHash;
    Vec3 translation = Vec3::Zero();
    double scale = 1.0;
    bool rotated = false;                 ///< S6a 恒为 false
    int nodeBegin = 0, nodeEnd = 0;       ///< 全局节点区间 [begin, end)
    int tetBegin = 0, tetEnd = 0;
    int surfaceBegin = 0, surfaceEnd = 0;
    int dbcCount = 0;
    double volume = 0.0;
    double mass = 0.0;
    std::array<double, 6> bbox{};         ///< 该体参考位形的全局 bbox

    int numNodes() const { return nodeEnd - nodeBegin; }
    int numTets() const { return tetEnd - tetBegin; }
};

struct MultiBodyScene {
    Positions x0;                              ///< 全局参考位置 [Ng,3]
    std::vector<std::array<int, 4>> tets;      ///< 全局 tet（0-based 全局节点）
    std::vector<std::array<int, 3>> surface;   ///< 全局表面三角面
    Eigen::VectorXd masses;                    ///< 全局集中质量 [Ng]
    std::vector<int> nodeBody;                 ///< [Ng] 节点 -> body_id
    std::vector<int> tetBody;                  ///< [Nt] tet -> body_id
    std::vector<int> dbcIndices;               ///< 全局 DBC 节点（按 body 拼接）
    std::vector<Vec3> dbcTargets;              ///< 与 dbcIndices 对应的绝对目标（= 参考位置）
    std::vector<BodyRecord> bodies;
    double totalVolume = 0.0;
    double totalMass = 0.0;

    int numNodes() const { return static_cast<int>(x0.rows()); }
    int numTets() const { return static_cast<int>(tets.size()); }
    int numBodies() const { return static_cast<int>(bodies.size()); }

    /// 组装全局 DBCInfo（每体顶部 cap 固定在参考位置）
    DBCInfo makeDBC(double stiffness) const {
        DBCInfo dbc;
        dbc.stiffness = stiffness;
        dbc.indices = dbcIndices;
        dbc.targetPositions = dbcTargets;
        dbc.targetVelocities.assign(dbcIndices.size(), Vec3::Zero());
        dbc.limits.assign(dbcIndices.size(), Vec3::Zero());
        dbc.isDBC.assign(static_cast<size_t>(numNodes()), 0);
        for (int i : dbcIndices) dbc.isDBC[i] = 1;
        return dbc;
    }
};

/// 装配全局多体场景。
/// @param mesh        单一资产（所有体共用同一拓扑）
/// @param assetHash   资产哈希（写入体记录）
/// @param transforms  每体的纯平移；空则退化为 1 体（平移 0）
/// @param density     密度 kg/m^3
/// @param dbcCapZMin  局部坐标 z >= 该值的节点被固定（顶部 cap）
inline MultiBodyScene buildMultiBodyScene(const MshMesh& mesh,
                                          const std::string& assetHash,
                                          const std::vector<BodyTransform>& transforms,
                                          double density,
                                          double dbcCapZMin) {
    if (mesh.empty() || mesh.tets.empty())
        throw std::runtime_error("buildMultiBodyScene: 空网格或没有四面体");
    if (mesh.surface.empty())
        throw std::runtime_error("buildMultiBodyScene: 缺少表面三角面");

    const MassSpringModel base = buildMassSpringModel(mesh, density);
    const int nLocal = base.numNodes();
    const int tLocal = mesh.numTets();
    const int sLocal = mesh.numSurface();

    std::vector<BodyTransform> xf = transforms;
    if (xf.empty()) xf.push_back(BodyTransform{});

    MultiBodyScene scene;
    const int nb = static_cast<int>(xf.size());
    scene.bodies.reserve(nb);
    scene.x0.resize(nb * nLocal, kDim);
    scene.masses.resize(nb * nLocal);
    scene.nodeBody.assign(static_cast<size_t>(nb * nLocal), 0);
    scene.tets.reserve(static_cast<size_t>(nb) * tLocal);
    scene.tetBody.reserve(static_cast<size_t>(nb) * tLocal);
    scene.surface.reserve(static_cast<size_t>(nb) * sLocal);

    for (int b = 0; b < nb; ++b) {
        const Vec3 t = xf[b].translation;
        BodyRecord rec;
        rec.bodyId = b;
        rec.assetHash = assetHash;
        rec.translation = t;
        rec.scale = 1.0;
        rec.rotated = false;
        rec.nodeBegin = b * nLocal;
        rec.nodeEnd = rec.nodeBegin + nLocal;
        rec.tetBegin = b * tLocal;
        rec.tetEnd = rec.tetBegin + tLocal;
        rec.surfaceBegin = b * sLocal;
        rec.surfaceEnd = rec.surfaceBegin + sLocal;
        rec.volume = base.volume;
        rec.mass = base.totalMass;

        // 节点：纯平移；质量：同一资产按体复制
        for (int i = 0; i < nLocal; ++i) {
            const int g = rec.nodeBegin + i;
            scene.x0.row(g) = (base.x0.row(i).transpose() + t).transpose();
            scene.masses[g] = base.masses[i];
            scene.nodeBody[g] = b;
        }
        // tet / 表面：索引整体平移，保证每个 tet 的 4 个节点同 body
        for (int e = 0; e < tLocal; ++e) {
            const auto& lt = mesh.tets[e];
            scene.tets.push_back({lt[0] + rec.nodeBegin, lt[1] + rec.nodeBegin,
                                  lt[2] + rec.nodeBegin, lt[3] + rec.nodeBegin});
            scene.tetBody.push_back(b);
        }
        for (int f = 0; f < sLocal; ++f) {
            const auto& lf = mesh.surface[f];
            scene.surface.push_back({lf[0] + rec.nodeBegin, lf[1] + rec.nodeBegin,
                                     lf[2] + rec.nodeBegin});
        }
        // DBC：局部 z >= cap 的节点，固定在该体的参考位置
        for (int i = 0; i < nLocal; ++i) {
            if (mesh.x(i, 2) >= dbcCapZMin) {
                const int g = rec.nodeBegin + i;
                scene.dbcIndices.push_back(g);
                scene.dbcTargets.push_back(scene.x0.row(g).transpose());
                ++rec.dbcCount;
            }
        }
        rec.bbox = boundingBox(Eigen::Map<const Positions>(
            scene.x0.data() + static_cast<Eigen::Index>(rec.nodeBegin) * kDim,
            nLocal, kDim));
        scene.totalVolume += rec.volume;
        scene.totalMass += rec.mass;
        scene.bodies.push_back(rec);
    }
    return scene;
}

} // namespace fiction
