#pragma once
// =============================================================================
// GravitationEnergy — 万有引力（纯 3D）
// =============================================================================
// 设计（对标 Python/UniversalGravitationEnergy.py）：
//   - 一组"吸引子"（如太阳）固定于 p_j，质量 M_j，位置不参与优化
//   - 其余节点（如地球）质量为 m_i，位置 x_i 参与优化
//   - 势能 E  = sum_i sum_j  -G * M_j * m_i / ||x_i - p_j||_2
//   - 梯度    dE/dx_i = sum_j  mu * d / r^3，  d = x_i - p_j， r = ||d||， mu = G*M_j*m_i
//   - Hessian d2E/dx_i^2 = sum_j [ mu/r^3 * I - 3*mu/r^5 * d d^T ]
//       精确 Hessian 非半正定（径向特征值 -2*mu/r^3，切向 +mu/r^3），
//       按框架惯例先 makePSD 再入 Newton 系统（与 Python project_hess_to_PSD=True 一致）
//   - 力 F = -grad 指向吸引子，大小 G*M_j*m_i / r^2
// =============================================================================
#include "energy.h"
#include "utils.h"

#include <cmath>
#include <utility>
#include <vector>

namespace fiction {

class GravitationEnergy final : public IEnergy {
public:
    struct Attractor {
        double mass;   ///< 吸引子质量 M_j（如太阳）
        Vec3 pos;      ///< 吸引子固定位置 p_j
    };

    /// @param masses     每个节点的质量 m_i，size N（与 TimeIntegrator 的 masses 一致）
    /// @param attractors 固定吸引子列表（点质量）
    GravitationEnergy(const Eigen::VectorXd& masses, std::vector<Attractor> attractors)
        : masses_(masses), attractors_(std::move(attractors)) {}

    const char* name() const override { return "GravitationEnergy"; }
    double hessianWeight(double h) const override { return h * h; }

    double val(const Positions& x) const override {
        const int N = static_cast<int>(masses_.size());
        double E = 0.0;
        for (int i = 0; i < N && i < (int)x.rows(); ++i) {
            const Vec3 xi = x.row(i).transpose();
            for (const Attractor& a : attractors_) {
                const Vec3 d = xi - a.pos;
                const double r = d.norm();
                if (r < kEps) continue;          // 防 r=0 除零
                E += -kG * a.mass * masses_[i] / r;
            }
        }
        return E;
    }

    void addGrad(const Positions& x, Positions& grad) const override {
        const int N = static_cast<int>(masses_.size());
        for (int i = 0; i < N && i < (int)grad.rows(); ++i) {
            const Vec3 xi = x.row(i).transpose();
            Vec3 g = Vec3::Zero();
            for (const Attractor& a : attractors_) {
                const Vec3 d = xi - a.pos;
                const double r = d.norm();
                if (r < kEps) continue;
                const double mu = kG * a.mass * masses_[i];
                g += (mu / (r * r * r)) * d;     // ∇(-mu/r) = +mu/r^3 * d，符号为正
            }
            grad.row(i) += g.transpose();
        }
    }

    void addHess(const Positions& x, std::vector<Triplet>& triplets) const override {
        const int N = static_cast<int>(masses_.size());
        for (int i = 0; i < N && i < (int)x.rows(); ++i) {
            const Vec3 xi = x.row(i).transpose();
            Eigen::Matrix3d Hi = Eigen::Matrix3d::Zero();
            for (const Attractor& a : attractors_) {
                const Vec3 d = xi - a.pos;
                const double r = d.norm();
                if (r < kEps) continue;
                const double mu = kG * a.mass * masses_[i];
                const double r3 = r * r * r;
                const double r5 = r3 * r * r;
                Hi += (mu / r3) * Eigen::Matrix3d::Identity()
                    - (3.0 * mu / r5) * (d * d.transpose());
            }
            Hi = makePSD(Hi);                    // 非 SPD -> PSD 投影（框架惯例）
            for (int rI = 0; rI < kDim; ++rI)
                for (int cI = 0; cI < kDim; ++cI)
                    triplets.emplace_back(kDim * i + rI, kDim * i + cI, Hi(rI, cI));
        }
    }

    double initStepSize(const Positions& /*x*/, const Positions& /*p*/) const override {
        return 1.0;
    }

private:
    static constexpr double kG   = 6.674e-11;    ///< 万有引力常量（与 Python 一致）
    static constexpr double kEps = 1e-12;        ///< 距离下限，防除零

    Eigen::VectorXd masses_;                     ///< 节点质量 m_i，size N
    std::vector<Attractor> attractors_;          ///< 固定吸引子（点质量）
};

} // namespace fiction
