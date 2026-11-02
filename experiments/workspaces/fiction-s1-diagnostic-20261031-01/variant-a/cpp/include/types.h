#pragma once
// =============================================================================
// types.h — fiction 纯 3D 公共类型定义（2D 适配器锁 z 版）
// =============================================================================
// 变更说明（纯 3D 主线）：
//   - 废弃运行时 Dim 推断，固定 kDim=3，Positions=[N,3] RowMajor
//   - 2D 视为 z=0 退化，需经 upgrade2DTo3D 补零并配合 lockZAxis 锁住 z 自由度
//   - DBCInfo 不再含 dim 字段，target/limits 固定 Vec3
//   - flatten/unflatten 基于 RowMajor 零拷贝 Map，移除动态 Dim 重载
// =============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cassert>

namespace fiction {

// -----------------------------------------------------------------------------
// 维度常量
// -----------------------------------------------------------------------------
constexpr int kDim = 3;

// -----------------------------------------------------------------------------
// 基础向量 / 矩阵别名
// -----------------------------------------------------------------------------

/// 单节点三维向量
using Vec3 = Eigen::Vector3d;
/// 遗留二维别名，仅作 2D 适配器入参
using Vec2 = Eigen::Vector2d;
/// 动态向量（遗留兼容，新代码请用 Vec3）
using VecX = Eigen::VectorXd;

/// 所有节点的位置矩阵，固定 [N,3] 行主序，每行一个节点
/// 内存布局 RowMajor：[x0,y0,z0, x1,y1,z1, ...] 与 FlatVector 连续一致
using Positions = Eigen::Matrix<double, Eigen::Dynamic, kDim, Eigen::RowMajor>;

/// 所有节点的速度矩阵，固定 [N,3]
using Velocities = Positions;

/// 展平后的自由度向量，维度 [3*N]，排序 [x0,y0,z0, x1,y1,z1, ...]
using FlatVector = Eigen::VectorXd;

/// 稀疏矩阵类型
using SparseMatrix = Eigen::SparseMatrix<double>;

/// 三元组
using Triplet = Eigen::Triplet<double>;

/// 质量向量 [N]
using MassVector = Eigen::VectorXd;

// -----------------------------------------------------------------------------
// 维度辅助 — 纯 3D 版固定为 3
// -----------------------------------------------------------------------------
inline int spatialDim(const Positions& /*m*/) { return kDim; }
inline int spatialDim(const Eigen::MatrixXd& m) {
    // 兼容旧调用：若传入动态矩阵，仍返回 3 并断言列数为 3 或 2（2D 适配）
    if (m.cols() == 2 || m.cols() == 3) return kDim;
    return kDim;
}
inline int numNodes(const Positions& m) { return static_cast<int>(m.rows()); }
inline int numNodes(const Eigen::MatrixXd& m) { return static_cast<int>(m.rows()); }

// -----------------------------------------------------------------------------
// 2D ↔ 3D 适配器
// -----------------------------------------------------------------------------
/// [N,2] → [N,3] 补零（z=0），用于复用旧 2D 网格/算例
inline Positions upgrade2DTo3D(const Eigen::MatrixXd& x2) {
    assert(x2.cols() == 2 && "upgrade2DTo3D expects [N,2]");
    const int N = static_cast<int>(x2.rows());
    Positions x3(N, kDim);
    for (int i = 0; i < N; ++i) {
        x3(i, 0) = x2(i, 0);
        x3(i, 1) = x2(i, 1);
        x3(i, 2) = 0.0;
    }
    return x3;
}

/// [N,3] → [N,2] 丢弃 z，用于离线导出/对比
inline Eigen::MatrixXd downgrade3DTo2D(const Positions& x3) {
    const int N = static_cast<int>(x3.rows());
    Eigen::MatrixXd x2(N, 2);
    for (int i = 0; i < N; ++i) {
        x2(i, 0) = x3(i, 0);
        x2(i, 1) = x3(i, 1);
    }
    return x2;
}

/// Vec2 → Vec3 补零
inline Vec3 upgradeVec2To3(const Vec2& v2) { return Vec3(v2.x(), v2.y(), 0.0); }
inline Vec3 upgradeVec2To3(const Eigen::VectorXd& v) {
    if (v.size() >= 3) return Vec3(v[0], v[1], v[2]);
    if (v.size() == 2) return Vec3(v[0], v[1], 0.0);
    if (v.size() == 1) return Vec3(v[0], 0.0, 0.0);
    return Vec3::Zero();
}

// -----------------------------------------------------------------------------
// 狄利克雷边界条件 — 纯 3D 版
// -----------------------------------------------------------------------------
struct DBCInfo {
    std::vector<int> indices;                 ///< 被约束的节点索引
    std::vector<Vec3> targetVelocities;       ///< 每个约束节点的目标速度 [3]
    std::vector<Vec3> limits;                 ///< 每个约束节点的限位 [3]
    std::vector<char> isDBC;                  ///< 节点级标记 size N
    double stiffness = 1e6;                   ///< 罚刚度
    /// 可选：绝对目标位置轨迹（如旋转驱动 DBC）。非空且 size==indices.size() 时，
    /// 每步 DBC 目标直接取该值（优先于 velocity/limit 机制）；可逐帧在调用方更新。
    std::vector<Vec3> targetPositions;

    // --- 便捷添加 ---
    void addDBC(int nodeIdx, const Vec3& vel, const Vec3& lim) {
        indices.push_back(nodeIdx);
        targetVelocities.push_back(vel);
        limits.push_back(lim);
    }
    void addDBC(int nodeIdx, const Vec2& vel2, const Vec2& lim2) {
        addDBC(nodeIdx, upgradeVec2To3(vel2), upgradeVec2To3(lim2));
    }
    void addDBC(int nodeIdx, const VecX& velX, const VecX& limX) {
        addDBC(nodeIdx, upgradeVec2To3(velX), upgradeVec2To3(limX));
    }
    /// 兼容旧 setFromVec 接口（现统一转 3D）
    void setFromVec2(const std::vector<int>& idx,
                     const std::vector<Vec2>& vels,
                     const std::vector<Vec2>& lims) {
        indices = idx;
        targetVelocities.resize(vels.size());
        for (size_t i = 0; i < vels.size(); ++i) targetVelocities[i] = upgradeVec2To3(vels[i]);
        limits.resize(lims.size());
        for (size_t i = 0; i < lims.size(); ++i) limits[i] = upgradeVec2To3(lims[i]);
    }
    void setFromVec3(const std::vector<int>& idx,
                     const std::vector<Vec3>& vels,
                     const std::vector<Vec3>& lims) {
        indices = idx;
        targetVelocities = vels;
        limits = lims;
    }
    void ensureIsDBCSize(int N) {
        if ((int)isDBC.size() != N) isDBC.assign(N, 0);
        for (int id : indices) if (id >= 0 && id < N) isDBC[id] = 1;
    }
};

// 2D 适配器：z 轴锁定请使用 energies/zlock_energy.h 的 ZLockEnergy（仅罚 z 分量，
// 不触发 DBC 自由度消去，x/y 运动不受影响）。整节点 DBC 会把 x/y 一并锁死，勿用于 z 锁。

// -----------------------------------------------------------------------------
// 仿真参数 — 纯 3D 固定 dim=3
// -----------------------------------------------------------------------------
struct SimParams {
    double h = 0.01;
    double tol = 1e-2;
    int maxNewtonIter = 100;
    static constexpr int dim = kDim;
};

// -----------------------------------------------------------------------------
// 布局转换 — 纯 3D RowMajor 零拷贝 Map
// -----------------------------------------------------------------------------
/// [N,3] → [3*N] 行主序展平
inline FlatVector flatten(const Positions& m) {
    const int N = static_cast<int>(m.rows());
    FlatVector v(static_cast<Eigen::Index>(N * kDim));
    // RowMajor 的 data 已是 [x0,y0,z0, x1,y1,z1...] 顺序，可直接 Map 拷贝
    Eigen::Map<const FlatVector> mp(m.data(), m.size());
    v = mp;
    return v;
}
/// 重载：兼容旧 Eigen::MatrixXd 入参（用于 2D 适配器内部）
inline FlatVector flatten(const Eigen::MatrixXd& m) {
    // 动态矩阵可能是 ColMajor，需按行主序手动展开以保持 [x0,y0,z0...]
    const int N = static_cast<int>(m.rows());
    const int D = static_cast<int>(m.cols());
    if (D == kDim) {
        // 若是 3 列，按 3*N 行主序展开
        FlatVector v(static_cast<Eigen::Index>(N * kDim));
        for (int i = 0; i < N; ++i)
            for (int d = 0; d < kDim; ++d)
                v(kDim * i + d) = m(i, d);
        return v;
    } else if (D == 2) {
        // 2D 适配：补零展平为 3*N
        FlatVector v(static_cast<Eigen::Index>(N * kDim));
        for (int i = 0; i < N; ++i) {
            v(kDim * i + 0) = m(i, 0);
            v(kDim * i + 1) = m(i, 1);
            v(kDim * i + 2) = 0.0;
        }
        return v;
    }
    FlatVector v(static_cast<Eigen::Index>(N * D));
    for (int i = 0; i < N; ++i)
        for (int d = 0; d < D; ++d)
            v(D * i + d) = m(i, d);
    return v;
}

/// [3*N] → [N,3]
inline Positions unflatten(const FlatVector& v) {
    const int total = static_cast<int>(v.size());
    const int N = total / kDim;
    Positions m(N, kDim);
    Eigen::Map<FlatVector> mp(m.data(), m.size());
    const int copy = std::min<int>(total, static_cast<int>(m.size()));
    for (int i = 0; i < copy; ++i) mp[i] = v[i];
    return m;
}
/// 兼容旧接口 unflatten(v, Dim) — 忽略 Dim，固定 3
inline Positions unflatten(const FlatVector& v, int /*Dim*/) {
    return unflatten(v);
}
/// 兼容旧接口 reshape — 固定 3
inline Positions reshape(const FlatVector& v, int N, int D) {
    (void)D; // 强制 3
    Positions m(N, kDim);
    const int total = static_cast<int>(v.size());
    const int copy = std::min<int>(total, N * kDim);
    Eigen::Map<FlatVector> mp(m.data(), m.size());
    for (int i = 0; i < copy; ++i) mp[i] = v[i];
    return m;
}
inline Positions reshape(const FlatVector& v, int D) {
    (void)D;
    return unflatten(v);
}

}  // namespace fiction
