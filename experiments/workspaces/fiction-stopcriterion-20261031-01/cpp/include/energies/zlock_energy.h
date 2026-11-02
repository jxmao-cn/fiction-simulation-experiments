#pragma once
// =============================================================================
// zlock_energy.h — z 轴锁定能量（纯 3D 2D 退化适配器）
// =============================================================================
// 用途：
//   2D 算例经 upgrade2DTo3D 补零后，为防止 z 自由度被数值误差激发，
//   对每个节点加一个仅作用于 z 分量的罚项：
//     E = 0.5 * k * sum_i (x(i,2))^2
//   grad 仅 z 分量，hess 仅 (3*i+2) 对角，权重固定 1.0（不经 h² 缩放）
// 注意：这是能量而非 DBC，不会触发自由度消去，x/y 运动完全不受影响。
// =============================================================================

#include "energy.h"

namespace fiction {

class ZLockEnergy final : public IEnergy {
public:
    explicit ZLockEnergy(const Eigen::VectorXd& masses, double stiffness = 1e8)
        : masses_(masses), k_(stiffness) {}

    const char* name() const override { return "ZLock"; }
    double hessianWeight(double /*h*/) const override { return 1.0; }

    double val(const Positions& x) const override {
        const int N = static_cast<int>(masses_.size());
        double E = 0.0;
        for (int i = 0; i < N && i < (int)x.rows(); ++i) {
            E += 0.5 * k_ * x(i,2) * x(i,2);
        }
        return E;
    }

    void addGrad(const Positions& x, Positions& grad) const override {
        const int N = static_cast<int>(masses_.size());
        for (int i = 0; i < N && i < (int)grad.rows(); ++i) {
            grad(i,2) += k_ * x(i,2);
        }
        (void)x;
    }

    void addHess(const Positions& x, std::vector<Triplet>& triplets) const override {
        const int N = static_cast<int>(masses_.size());
        for (int i = 0; i < N && i < (int)x.rows(); ++i) {
            triplets.emplace_back(3*i+2, 3*i+2, k_);
        }
    }

private:
    Eigen::VectorXd masses_;   // 保留以对齐能量接口风格（质量加权可后续引入）
    double k_;
};

/// 便捷工厂：为 N 个节点创建 z 锁能量
inline EnergyPtr makeZLockEnergy(const Eigen::VectorXd& masses, double stiffness = 1e8) {
    return std::make_shared<ZLockEnergy>(masses, stiffness);
}

} // namespace fiction
