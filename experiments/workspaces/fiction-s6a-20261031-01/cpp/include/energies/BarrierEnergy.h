#pragma once
// =============================================================================
// BarrierEnergy.h — Barrier 接触能量（纯 3D 版），对标 Python/BarrierEnergy.py
// =============================================================================
// P_b(x) = Σ_i A_i * dhat * b(d_i)，其中
//   b(d)     = κ/2 * (s-1) * ln(s)，s = d/dhat        （对标 py:14）
//   b'(d)    = κ/2 * (ln(s)/dhat + (s-1)/d)           （对标 py:32 括号项）
//   b''(d)   = κ * (d+dhat) / (2*d²*dhat)             （对标 py:51 系数）
// 仅当 0 < d < dhat 时激活；d ≥ dhat 时能量/梯度/Hessian 全为 0。
//
// 两种平面（对标 Python 的 floor / ceiling 两段）：
//   - 固定平面（floor）  ：d_i = n·(x_i - o)
//   - 节点锚定平面（ceiling）：d_i = n·(x_i - x_driver)，driver 为锚定节点
//     （Python 中即 x[-1] 驱动的天花板，n=(0,-1)）
//
// 平面有向距离 d 是 x 的线性函数 ⇒ grad d = n（常单位向量），hess d = 0，故：
//   grad：g_i = A_i*dhat*b'(d)*n；ceiling 时 driver 节点取等值反向（对标 py:40-41）
//   hess：floor 为单节点 3×3 块 A_i*dhat*b''(d)*n nᵀ（对标 py:52-56）；
//         ceiling 为 (i, driver) 的 6×6 块 [[H,-H],[-H,H]]（对标 py:62-70）
//
// initStepSize：数值 CCD（对标 py:74-88 init_step_size）——
//   对所有靠近平面且朝平面运动的节点，限制 alpha 使 d + alpha*(n·p) ≥ 0.1*d，
//   安全系数 0.9，保证线搜索内 log 定义域（d>0）不被破坏。
//
// 与 Python 的差异（妥善完善项）：
//   - grad 严格累加（+=），符合 IEnergy 多能量叠加约定（Python floor 段为赋值，
//     因其 g 初始为零且每个节点至多命中一个平面）
//   - 激活条件加 d > 0 下界保护，避免 log(≤0) 产生 NaN（Python 无保护）
// =============================================================================

#include "energy.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace fiction {

class BarrierEnergy final : public IEnergy {
public:
    // -----------------------------------------------------------------------------
    // 平面定义：固定平面（driverIdx<0）用 origin；节点锚定平面（driverIdx>=0）
    // 的原点取该节点当前位置，随仿真移动
    // -----------------------------------------------------------------------------
    struct Plane {
        Vec3 n;              // 单位法向量，d = n·(x - o)，安全侧为正
        Vec3 origin;         // 平面上一点（driverIdx >= 0 时忽略）
        int  driverIdx = -1; // >= 0：平面锚定在该节点（ceiling）；-1：固定平面（floor）
    };

    // -----------------------------------------------------------------------------
    // @param contact_area 每节点接触面积权重 A_i，size N
    // @param planes       平面列表（floor/ceiling 混合）
    // @param dhat         激活距离（对标 Python dhat = 0.01）
    // @param kappa        罚刚度（对标 Python kappa = 1e5）
    // -----------------------------------------------------------------------------
    BarrierEnergy(const Eigen::VectorXd& contact_area,
                  const std::vector<Plane>& planes,
                  double dhat = 0.01,
                  double kappa = 1e5)
        : contact_area_(contact_area), planes_(planes), dhat_(dhat), kappa_(kappa) {}

    const char* name() const override { return "BarrierEnergy"; }
    double hessianWeight(double h) const override { return h * h; }

    // -----------------------------------------------------------------------------
    // val — 对标 py:7-23 val（floor + ceiling 求和）
    // -----------------------------------------------------------------------------
    double val(const Positions& x) const override {
        const int N = nodeCount(x);
        double E = 0.0;
        for (const Plane& pl : planes_) {
            for (int i = 0; i < N; ++i) {
                if (i == pl.driverIdx) continue; // ceiling 锚点自身不参与（对标 py:17）
                const double d = signedDist(pl, x, i);
                if (d > 0.0 && d < dhat_) {
                    E += contact_area_[i] * dhat_ * b(d);
                }
            }
        }
        return E;
    }

    // -----------------------------------------------------------------------------
    // addGrad — 对标 py:25-43 grad；grad d = n，故 g_i = A_i*dhat*b'(d)*n
    // -----------------------------------------------------------------------------
    void addGrad(const Positions& x, Positions& grad) const override {
        const int N = nodeCount(x);
        for (const Plane& pl : planes_) {
            for (int i = 0; i < N; ++i) {
                if (i == pl.driverIdx) continue;
                const double d = signedDist(pl, x, i);
                if (d <= 0.0 || d >= dhat_) continue;
                const Vec3 g = (contact_area_[i] * dhat_ * bPrime(d)) * pl.n;
                grad.row(i) += g.transpose();
                if (pl.driverIdx >= 0) {
                    grad.row(pl.driverIdx) -= g.transpose(); // 反作用（对标 py:41）
                }
            }
        }
    }

    // -----------------------------------------------------------------------------
    // addHess — 对标 py:45-72 hess（IJV scatter → 此处 Triplet scatter）
    // -----------------------------------------------------------------------------
    void addHess(const Positions& x, std::vector<Triplet>& triplets) const override {
        const int N = nodeCount(x);
        for (const Plane& pl : planes_) {
            for (int i = 0; i < N; ++i) {
                if (i == pl.driverIdx) continue;
                const double d = signedDist(pl, x, i);
                if (d <= 0.0 || d >= dhat_) continue;
                // d 线性 ⇒ hess d = 0，局部块 = A_i*dhat*b''(d)*n nᵀ
                const Eigen::Matrix3d H3 = (contact_area_[i] * dhat_ * bDoublePrime(d))
                                         * (pl.n * pl.n.transpose());
                if (pl.driverIdx < 0) {
                    // floor：单节点 3×3 对角块
                    const Eigen::Matrix3d H = makePSD(H3);
                    for (int r = 0; r < 3; ++r)
                        for (int c = 0; c < 3; ++c)
                            triplets.emplace_back(3 * i + r, 3 * i + c, H(r, c));
                } else {
                    // ceiling：(i, driver) 的 6×6 边块 [[H,-H],[-H,H]]
                    Eigen::Matrix<double, 6, 6> Hl;
                    Hl <<  H3, -H3,
                          -H3,  H3;
                    Hl = makePSD(Hl);
                    const int idx[2] = { i, pl.driverIdx };
                    for (int nI = 0; nI < 2; ++nI)
                        for (int nJ = 0; nJ < 2; ++nJ)
                            for (int r = 0; r < 3; ++r)
                                for (int c = 0; c < 3; ++c)
                                    triplets.emplace_back(3 * idx[nI] + r, 3 * idx[nJ] + c,
                                                          Hl(nI * 3 + r, nJ * 3 + c));
                }
            }
        }
    }

    // -----------------------------------------------------------------------------
    // initStepSize — 数值 CCD，对标 py:74-88 init_step_size
    //   pn = n·p（ceiling 用相对速度 n·(p_i - p_driver)）< 0（朝平面运动）时，
    //   alpha ≤ 0.9 * d / (-pn)，保证 x + alpha*p 不穿透平面
    // -----------------------------------------------------------------------------
    double initStepSize(const Positions& x, const Positions& p) const override {
        const int N = nodeCount(x);
        double alpha = 1.0;
        for (const Plane& pl : planes_) {
            for (int i = 0; i < N; ++i) {
                if (i == pl.driverIdx) continue;
                Vec3 relP = p.row(i).transpose();
                if (pl.driverIdx >= 0) relP -= p.row(pl.driverIdx).transpose();
                const double pn = pl.n.dot(relP);
                if (pn < 0.0) {
                    const double d = signedDist(pl, x, i);
                    if (d > 0.0) alpha = std::min(alpha, 0.9 * d / -pn);
                }
            }
        }
        return alpha;
    }

private:
    int nodeCount(const Positions& x) const {
        return std::min<int>(static_cast<int>(x.rows()),
                             static_cast<int>(contact_area_.size()));
    }

    Vec3 planeOrigin(const Plane& pl, const Positions& x) const {
        if (pl.driverIdx >= 0) return x.row(pl.driverIdx).transpose();
        return pl.origin;
    }

    // 有向距离 d = n·(x_i - o)，对标 py:11 / py:18
    double signedDist(const Plane& pl, const Positions& x, int i) const {
        return pl.n.dot(x.row(i).transpose() - planeOrigin(pl, x));
    }

    // b(d) = κ/2 * (s-1) * ln(s)，s = d/dhat（仅 0 < d < dhat 调用）
    double b(double d) const {
        const double s = d / dhat_;
        return 0.5 * kappa_ * (s - 1.0) * std::log(s);
    }

    // b'(d) = κ/2 * (ln(s)/dhat + (s-1)/d)（乘积法则；1/(s*dhat) = 1/d）
    double bPrime(double d) const {
        const double s = d / dhat_;
        return 0.5 * kappa_ * (std::log(s) / dhat_ + (s - 1.0) / d);
    }

    // b''(d) = κ * (d+dhat) / (2*d²*dhat)，与 py:51 系数一致
    double bDoublePrime(double d) const {
        return kappa_ * (d + dhat_) / (2.0 * d * d * dhat_);
    }

    Eigen::VectorXd contact_area_;   ///< A_i，size N
    std::vector<Plane> planes_;      ///< floor / ceiling 平面集合
    double dhat_;                    ///< 激活距离
    double kappa_;                   ///< 罚刚度
};

} // namespace fiction
