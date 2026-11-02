#pragma once
// =============================================================================
// energies/RegionForce.h — 坐标盒区域恒定外力能量（对标 IPC Neumann BC）
// =============================================================================
// 物理含义：对指定节点集合施加恒定加速度 a（如拉力/侧推），
// 等价于 IPC Optimizer.cpp 中的 NBC 线性项：
//   E = -sum_{i in region} m_i * a . x_i
//   grad_i = -m_i * a（常数）
//   Hess_i = 0（以 hessianWeight = h^2 进入 IP，h^2 * grad 对应 -h^2 * m_i * a）
// 与 Gravity 同构：重力即"全节点 + g"的特例。
//
// 注：Hessian 为零是线性势能的正确形式；牛顿求解的 PSD 性由 Inertia 与
//     弹簧 Hessian 提供，外力只改变梯度方向，不会破坏收敛框架。
// =============================================================================

#include "energy.h"

#include <utility>
#include <vector>

namespace fiction {

class RegionForceEnergy final : public IEnergy {
public:
    /// @param masses 质量向量 [N]（与 TimeIntegrator 使用同一份）
    /// @param nodes  受力节点索引集合（如坐标盒筛选结果）
    /// @param accel  恒定加速度向量（力 = m_i * accel）
    RegionForceEnergy(const Eigen::VectorXd& masses,
                      std::vector<int> nodes,
                      const Vec3& accel)
        : masses_(masses), nodes_(std::move(nodes)), accel_(accel) {}

    const char* name() const override { return "RegionForce"; }
    double hessianWeight(double h) const override { return h * h; }

    /// 受力节点数量（报告用）
    int numNodes() const { return static_cast<int>(nodes_.size()); }
    const Vec3& accel() const { return accel_; }

    double val(const Positions& x) const override {
        double E = 0.0;
        for (int i : nodes_) {
            if (i < 0 || i >= x.rows() || i >= masses_.size()) continue;
            const Vec3 xi = x.row(i).transpose();
            E += -masses_[i] * xi.dot(accel_);
        }
        return E;
    }

    void addGrad(const Positions& /*x*/, Positions& grad) const override {
        for (int i : nodes_) {
            if (i < 0 || i >= grad.rows() || i >= masses_.size()) continue;
            grad.row(i) += (-masses_[i] * accel_).transpose();
        }
    }

    void addHess(const Positions& /*x*/, std::vector<Triplet>& /*triplets*/) const override {
        // 线性势能 Hessian 为零
    }

private:
    Eigen::VectorXd masses_;
    std::vector<int> nodes_;
    Vec3 accel_;
};

/// 便捷工厂
inline EnergyPtr makeRegionForceEnergy(const Eigen::VectorXd& masses,
                                       std::vector<int> nodes,
                                       const Vec3& accel) {
    return std::make_shared<RegionForceEnergy>(masses, std::move(nodes), accel);
}

} // namespace fiction
