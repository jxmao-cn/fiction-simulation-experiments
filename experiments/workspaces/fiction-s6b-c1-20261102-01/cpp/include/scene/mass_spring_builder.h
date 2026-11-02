#pragma once
// =============================================================================
// scene/mass_spring_builder.h — 从 MshMesh 装配 mass-spring 模型
// =============================================================================
// 职责：
//   - 提取四面体全部唯一边（每个 tet 6 条边）作为弹簧拓扑
//   - 集中质量：m_i = rho * sum_{tet ∋ i} V_tet / 4
//   - 静止长度平方 l2、每边刚度 k（显式或 auto = E * 平均边长）
//   - 坐标盒选点（供 pin-box DBC / force-box 区域外力使用）
//
// 单位约定：长度 m，密度 kg/m^3，刚度 k 为每边线刚度（N/m），
//           auto 模式 k ≈ E * l_avg 使量级与 IPC 的 Young's modulus 对齐。
// =============================================================================

#include "io/msh_loader.h"
#include "types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace fiction {

// -----------------------------------------------------------------------------
// 几何辅助
// -----------------------------------------------------------------------------
/// 单个四面体体积（绝对值，对节点定向不敏感）
inline double tetVolume(const Positions& x, const std::array<int, 4>& tet) {
    const Vec3 a = x.row(tet[0]).transpose();
    const Vec3 b = x.row(tet[1]).transpose();
    const Vec3 c = x.row(tet[2]).transpose();
    const Vec3 d = x.row(tet[3]).transpose();
    return std::abs((b - a).dot((c - a).cross(d - a))) / 6.0;
}

/// 坐标包围盒 [xmin, ymin, zmin, xmax, ymax, zmax]
inline std::array<double, 6> boundingBox(const Positions& x) {
    std::array<double, 6> box = {1e300, 1e300, 1e300, -1e300, -1e300, -1e300};
    for (int i = 0; i < x.rows(); ++i) {
        for (int d = 0; d < kDim; ++d) {
            box[d] = std::min(box[d], x(i, d));
            box[3 + d] = std::max(box[3 + d], x(i, d));
        }
    }
    return box;
}

/// 选出位于坐标盒 [xmin, ymin, zmin, xmax, ymax, zmax] 内的节点索引（闭区间）
inline std::vector<int> nodesInBox(const Positions& x, const std::array<double, 6>& box) {
    std::vector<int> nodes;
    for (int i = 0; i < x.rows(); ++i) {
        bool inside = true;
        for (int d = 0; d < kDim; ++d) {
            if (x(i, d) < box[d] || x(i, d) > box[3 + d]) {
                inside = false;
                break;
            }
        }
        if (inside) nodes.push_back(i);
    }
    return nodes;
}

/// 模型包围盒中心（用于 viewer 相机对焦）
inline Vec3 boxCenter(const std::array<double, 6>& box) {
    return Vec3(0.5 * (box[0] + box[3]),
                0.5 * (box[1] + box[4]),
                0.5 * (box[2] + box[5]));
}

/// 模型包围盒对角线长度（用于 viewer 相机距离）
inline double boxDiagonal(const std::array<double, 6>& box) {
    return std::sqrt((box[3] - box[0]) * (box[3] - box[0]) +
                     (box[4] - box[1]) * (box[4] - box[1]) +
                     (box[5] - box[2]) * (box[5] - box[2]));
}

// -----------------------------------------------------------------------------
// Mass-spring 模型
// -----------------------------------------------------------------------------
struct MassSpringModel {
    Positions x0;                            ///< [N,3] 参考位置
    Eigen::VectorXd masses;                  ///< [N] 集中质量
    std::vector<std::array<int, 2>> edges;   ///< [E,2] 唯一边（弹簧）
    Eigen::VectorXd l2;                      ///< [E] 静止长度平方
    Eigen::VectorXd k;                       ///< [E] 每边刚度
    double volume = 0.0;                     ///< 四面体总体积
    double totalMass = 0.0;                  ///< 总质量
    double avgEdgeLength = 0.0;              ///< 平均静止边长（auto 刚度 / 报告用）
    std::array<double, 6> bbox{};            ///< 包围盒

    int numNodes() const { return static_cast<int>(x0.rows()); }
    int numEdges() const { return static_cast<int>(edges.size()); }
};

/// 装配 mass-spring 模型
/// @param mesh      已加载的四面体网格
/// @param density   密度 kg/m^3（IPC 默认 1000）
/// @param stiffness 每边刚度；<= 0 时自动取 youngs * 平均边长
/// @param youngs    auto 刚度使用的 Young's modulus（IPC 默认 1e5 Pa）
inline MassSpringModel buildMassSpringModel(const MshMesh& mesh,
                                            double density = 1000.0,
                                            double stiffness = 0.0,
                                            double youngs = 1e5) {
    if (mesh.tets.empty()) throw std::runtime_error("buildMassSpringModel: mesh 无四面体单元");
    if (density <= 0.0) throw std::runtime_error("buildMassSpringModel: density 必须为正");

    MassSpringModel model;
    model.x0 = mesh.x;
    const int N = static_cast<int>(mesh.x.rows());
    const int T = mesh.numTets();

    // --- 集中质量：rho * V_tet / 4 分摊到 4 个节点 ---
    model.masses = Eigen::VectorXd::Zero(N);
    for (const auto& tet : mesh.tets) {
        const double V = tetVolume(model.x0, tet);
        model.volume += V;
        for (int i = 0; i < 4; ++i) model.masses[tet[i]] += 0.25 * density * V;
    }
    model.totalMass = model.masses.sum();
    const double massFloor = 1e-12 * model.totalMass / std::max(N, 1);
    for (int i = 0; i < N; ++i)
        if (model.masses[i] < massFloor) model.masses[i] = massFloor;

    // --- 唯一边提取（每个四面体 6 条边）---
    const std::array<std::array<int, 2>, 6> kTetEdges = {{
        {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
    std::unordered_set<std::uint64_t> seen;
    seen.reserve(static_cast<size_t>(T) * 6);
    auto keyOf = [](int a, int b) -> std::uint64_t {
        if (a > b) std::swap(a, b);
        return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
    };
    model.edges.reserve(static_cast<size_t>(T) * 6);
    for (const auto& tet : mesh.tets) {
        for (const auto& e : kTetEdges) {
            const int a = tet[e[0]], b = tet[e[1]];
            if (seen.insert(keyOf(a, b)).second) model.edges.push_back({a, b});
        }
    }

    // --- 静止长度、刚度 ---
    const int E = model.numEdges();
    model.l2.resize(E);
    model.k.resize(E);
    double sumLen = 0.0;
    for (int e = 0; e < E; ++e) {
        const Vec3 d = model.x0.row(model.edges[e][0]).transpose() -
                       model.x0.row(model.edges[e][1]).transpose();
        const double len2 = d.squaredNorm();
        model.l2[e] = len2 > 1e-24 ? len2 : 1e-24;
        sumLen += std::sqrt(model.l2[e]);
    }
    model.avgEdgeLength = E > 0 ? sumLen / E : 0.0;

    const double kValue = stiffness > 0.0 ? stiffness : youngs * model.avgEdgeLength;
    model.k.setConstant(E, kValue);

    model.bbox = boundingBox(model.x0);
    return model;
}

} // namespace fiction
