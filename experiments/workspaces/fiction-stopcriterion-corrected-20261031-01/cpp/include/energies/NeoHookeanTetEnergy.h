#pragma once
// =============================================================================
// energies/NeoHookeanTetEnergy.h — 四面体 P1 有限元 Neo-Hookean 能量（纯 3D 版）
// =============================================================================
// 背景：
//   复现 Elastogen 论文 cantilever 实验（docs/elastogen-rods-experiments.md）所需的
//   体弹性能量；公式与 IPC 的 NeoHookeanEnergy（SIGGRAPH 2020, Minchen Li）一致：
//     E = V0 * [ μ/2 (F:F − 3) − μ·log J + λ/2 (log J)² ]
//     P = ∂E/∂F = μ (F − F⁻ᵀ) + λ (log J) F⁻ᵀ
//   P1 常应变四面体：F = Ds · G，Ds = [x1−x0, x2−x0, x3−x0]（3×3 列块），
//   G = Dm⁻¹，Dm = [x̂1−x̂0, x̂2−x0, x̂3−x̂0]，V0 = |det Dm| / 6
//
// 梯度（节点 k = 1,2,3）：
//   dE/dx_k = V0 · (P · u_k)，u_k = G 的第 k 行（列向量）
//   dE/dx_0 = −Σ_{k=1..3} dE/dx_k
// Hessian（12×12，F 空间二阶契约到 x 空间）：
//   H_F[(i,j),(m,n)] = μ δ_{im} δ_{jn}
//                    + (μ − λ·logJ) Φ_{in} Φ_{mj}      （Φ = F⁻ᵀ）
//                    + λ Φ_{ij} Φ_{mn}
//   B: 9×12 算子（每个节点自由度 → dF 列展开），H_x = V0 · Bᵀ H_F B，
//   再做 makePSD（12×12）投影后进入全局装配。
//
// 与 mass-spring 的关系：
//   该能量带泊松比（λ 由 ν 决定），是有连续先验的体模型；
//   mass-spring（边弹簧）无泊松比，只能定性对照。
// =============================================================================

#include "energy.h"
#include "utils.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace fiction {

class NeoHookeanTetEnergy final : public IEnergy {
public:
    /// @param tets     四面体连接（0-based 顶点索引）
    /// @param x0       参考位置 [N,3]
    /// @param youngs   Young's modulus E (Pa)
    /// @param poisson  泊松比 ν ∈ (0, 0.5)
    NeoHookeanTetEnergy(const std::vector<std::array<int, 4>>& tets,
                        const Positions& x0,
                        double youngs,
                        double poisson)
        : tets_(tets) {
        // Lame 参数（3D 各向同性）：μ = E/(2(1+ν)), λ = Eν/((1+ν)(1−2ν))
        mu_ = youngs / (2.0 * (1.0 + poisson));
        lambda_ = youngs * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));

        const int T = static_cast<int>(tets_.size());
        v0_.assign(T, 0.0);
        gj_.assign(T, Eigen::Matrix3d::Identity());
        for (int t = 0; t < T; ++t) {
            const auto& tet = tets_[t];
            const Eigen::Vector3d p0 = x0.row(tet[0]).transpose();
            Eigen::Matrix3d Dm;
            Dm.col(0) = x0.row(tet[1]).transpose() - p0;
            Dm.col(1) = x0.row(tet[2]).transpose() - p0;
            Dm.col(2) = x0.row(tet[3]).transpose() - p0;
            const double detDm = Dm.determinant();
            const double vol = std::abs(detDm) / 6.0;
            if (vol > 1e-14) {
                v0_[t] = vol;
                gj_[t] = Dm.inverse();
            } else {
                v0_[t] = 0.0; // 退化四面体：能量恒为 0
            }
        }
    }

    const char* name() const override { return "NeoHookeanTet"; }
    double hessianWeight(double h) const override { return h * h; }

    double val(const Positions& x) const override {
        double E_sum = 0.0;
        for (size_t t = 0; t < tets_.size(); ++t) {
            if (v0_[t] <= 0.0) continue;
            const Eigen::Matrix3d F = computeF(x, t);
            const double J = F.determinant();
            const double Jc = std::max(J, kJClamp);
            const double lc = std::log(Jc);
            E_sum += v0_[t] * (0.5 * mu_ * (F.squaredNorm() - 3.0)
                               - mu_ * lc + 0.5 * lambda_ * lc * lc);
        }
        return E_sum;
    }

    void addGrad(const Positions& x, Positions& grad) const override {
        for (size_t t = 0; t < tets_.size(); ++t) {
            if (v0_[t] <= 0.0) continue;
            const auto& tet = tets_[t];
            const Eigen::Matrix3d F = computeF(x, t);
            const double J = F.determinant();
            const double lc = std::log(std::max(J, kJClamp));
            const Eigen::Matrix3d Phi = F.inverse().transpose(); // F⁻ᵀ
            const Eigen::Matrix3d P = mu_ * (F - Phi) + lambda_ * lc * Phi;

            // dE/dx_k (k=1..3) = V0 * P * u_k，u_k = G 的第 k 行
            Eigen::Vector3d g1 = v0_[t] * (P * gj_[t].row(0).transpose());
            Eigen::Vector3d g2 = v0_[t] * (P * gj_[t].row(1).transpose());
            Eigen::Vector3d g3 = v0_[t] * (P * gj_[t].row(2).transpose());
            const Eigen::Vector3d g0 = -(g1 + g2 + g3);

            for (int d = 0; d < 3; ++d) {
                grad(tet[0], d) += g0(d);
                grad(tet[1], d) += g1(d);
                grad(tet[2], d) += g2(d);
                grad(tet[3], d) += g3(d);
            }
        }
    }

    void addHess(const Positions& x, std::vector<Triplet>& trips) const override {
        for (size_t t = 0; t < tets_.size(); ++t) {
            if (v0_[t] <= 0.0) continue;
            Eigen::Matrix<double, 12, 12> H = makePSD(computeHBlock(x, t));
            scatterBlock(tets_[t], H, trips);
        }
    }

    /// 精确 Hessian（未做 makePSD 投影）——供有限差分校验；对压缩/剪切位形非 PSD 属正常
    void addHessRaw(const Positions& x, std::vector<Triplet>& trips) const {
        for (size_t t = 0; t < tets_.size(); ++t) {
            if (v0_[t] <= 0.0) continue;
            scatterBlock(tets_[t], computeHBlock(x, t), trips);
        }
    }

    double initStepSize(const Positions& /*x*/, const Positions& /*p*/) const override {
        return 1.0;
    }

private:
    static constexpr double kJClamp = 1e-12; // log(J) 数值下限

    Eigen::Matrix3d computeF(const Positions& x, size_t t) const {
        const auto& tet = tets_[t];
        const Eigen::Vector3d p0 = x.row(tet[0]).transpose();
        Eigen::Matrix3d Ds;
        Ds.col(0) = x.row(tet[1]).transpose() - p0;
        Ds.col(1) = x.row(tet[2]).transpose() - p0;
        Ds.col(2) = x.row(tet[3]).transpose() - p0;
        return Ds * gj_[t];
    }

    /// 单四面体 x 空间 12×12 精确 Hessian 块（未投影）
    Eigen::Matrix<double, 12, 12> computeHBlock(const Positions& x, size_t t) const {
        const Eigen::Matrix3d F = computeF(x, t);
        const double J = F.determinant();
        const double lc = std::log(std::max(J, kJClamp));
        const Eigen::Matrix3d Phi = F.inverse().transpose();

        // F 空间 9×9 Hessian（dF 用列展开 F_{mn} → m*3+n）
        Eigen::Matrix<double, 9, 9> Hf;
        Hf.setZero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int m = 0; m < 3; ++m)
                    for (int n = 0; n < 3; ++n) {
                        Hf(i * 3 + j, m * 3 + n) =
                            mu_ * double(i == m && j == n)
                            + (mu_ - lambda_ * lc) * Phi(i, n) * Phi(m, j)
                            + lambda_ * Phi(i, j) * Phi(m, n);
                    }

        // dF/dx 契约：dF(m,n) = Σ_k dDs(m,k) · G(k,n)
        //   a==0: dDs(m,k) = −δ_{m,r}（三列同为 −e_r）
        //   a≥1 : dDs(m,k) =  δ_{m,r}·[k == a−1]
        Eigen::Matrix<double, 9, 12> B;
        B.setZero();
        for (int a = 0; a < 4; ++a) {
            for (int r = 0; r < 3; ++r) {
                for (int m = 0; m < 3; ++m) {
                    if (m != r) continue;
                    for (int k = 0; k < 3; ++k) {
                        double dDs_mk;
                        if (a == 0) dDs_mk = -1.0;
                        else dDs_mk = (k == a - 1) ? 1.0 : 0.0;
                        for (int n = 0; n < 3; ++n) {
                            B(m * 3 + n, 3 * a + r) += dDs_mk * gj_[t](k, n);
                        }
                    }
                }
            }
        }
        const double v0 = v0_[t];
        return v0 * (B.transpose() * Hf * B);
    }

    static void scatterBlock(const std::array<int, 4>& tet,
                             const Eigen::Matrix<double, 12, 12>& H,
                             std::vector<Triplet>& trips) {
        for (int a = 0; a < 4; ++a) {
            for (int b = 0; b < 4; ++b) {
                for (int r = 0; r < 3; ++r) {
                    for (int c = 0; c < 3; ++c) {
                        trips.emplace_back(3 * tet[a] + r, 3 * tet[b] + c,
                                           H(3 * a + r, 3 * b + c));
                    }
                }
            }
        }
    }

    std::vector<std::array<int, 4>> tets_;
    std::vector<double> v0_;            // 每四面体体积（退化为 0）
    std::vector<Eigen::Matrix3d> gj_;   // Dm⁻¹
    double mu_ = 0.0;
    double lambda_ = 0.0;
};

} // namespace fiction
