#pragma once
// =============================================================================
// Gravity.h — 重力能量（纯 3D 版）
// =============================================================================
// E = sum -m_i * x_i.dot(g)，梯度常数，无 Hessian
// =============================================================================

#include "energy.h"
#include "utils.h"

namespace fiction {

// -----------------------------------------------------------------------------
// ExampleGravityEnergy — 示例能量（重力占位实现，纯 3D）
//   E = sum -m_i * x_i.dot(g)，梯度常数，无 Hessian
// -----------------------------------------------------------------------------
class ExampleGravityEnergy final : public IEnergy {
public:
    explicit ExampleGravityEnergy(const Eigen::VectorXd& masses, const Vec3& gravity)
        : masses_(masses), gravity_(gravity) {}
    // 默认重力：z 向下（也可传入 x/y 方向）
    explicit ExampleGravityEnergy(const Eigen::VectorXd& masses)
        : masses_(masses), gravity_(0.0, 0.0, -9.81) {}

    const char* name() const override { return "ExampleGravity"; }
    double hessianWeight(double h) const override { return h * h; }

    double val(const Positions& x) const override {
        const int N = static_cast<int>(masses_.size());
        double E = 0.0;
        for (int i = 0; i < N && i < (int)x.rows(); ++i) {
            Vec3 xi(x(i,0), x(i,1), x(i,2));
            E += -masses_[i] * xi.dot(gravity_);
        }
        return E;
    }

    void addGrad(const Positions& x, Positions& grad) const override {
        const int N = static_cast<int>(masses_.size());
        for (int i = 0; i < N && i < (int)grad.rows(); ++i) {
            grad.row(i) += (-masses_[i] * gravity_).transpose();
        }
        (void)x;
    }

    void addHess(const Positions& x, std::vector<Triplet>& triplets) const override {
        (void)x; (void)triplets;
        // 重力 Hessian 为 0
    }

    double initStepSize(const Positions& /*x*/, const Positions& /*p*/) const override {
        return 1.0;
    }

private:
    Eigen::VectorXd masses_;
    Vec3 gravity_;
};

} // namespace fiction