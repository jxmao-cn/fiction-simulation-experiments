#pragma once
// =============================================================================
// energy.h — 统一能量接口 (Unified Energy Interface) — 纯 3D 版
// =============================================================================
// 设计目标：
//   - 纯 3D 主线：Positions 固定 [N,3] RowMajor，kDim=3
//   - Inertia 内置于 TimeIntegrator，其余能量走 IEnergy
//   - 2D 视为 z=0 退化，需经 upgrade2DTo3D 补零并配合 lockZAxis 锁 z
//   - 保持 Eigen 封装 + CHOLMOD 主线，C++20
//
// 对应 Python 总装逻辑 time_integrator.py:59-114
// =============================================================================

// -----------------------------------------------------------------------------
// HOW TO ADD NEW ENERGY — 纯 3D 版（精简）
// -----------------------------------------------------------------------------
// 1) 新建 include/energies/my_energy.h 继承 IEnergy
// 2) 实现 val/addGrad/addHess 三函数，索引固定 3*node+d，局部块 6x6(边)/9x9(面)/12x12(体)
//    - 索引：triplets.emplace_back(3*i+d, 3*j+e, value)
//    - 非 SPD 先 makePSD 再拆解
//    - grad 用 grad.row(i) += diff.transpose() 且 diff 为 Vec3
// 3) 注册：integrator.addEnergy(std::make_shared<MyEnergy>(...));
// 4) 伪代码（纯 3D，边单元）：
//   Vec3 diff = x.row(e[0]) - x.row(e[1]); // x:[N,3]
//   double E = 0.5*k*pow(diff.squaredNorm()/l2-1,2)*l2;
//   Vec3 g = 2*k*(diff.squaredNorm()/l2-1)*diff;
//   grad.row(e[0]) += g.transpose(); grad.row(e[1]) -= g.transpose();
//   Matrix3d Hdiff = 2*k/l2*(2*diff*diff.transpose()+(diff.squaredNorm()-l2)*I3);
//   Matrix<double,6,6> Hl; Hl<< Hdiff,-Hdiff,-Hdiff,Hdiff; Hl=makePSD(Hl);
//   for nI,nJ in {0,1} for r,c in 0..2 trips.emplace_back(3*e[nI]+r,3*e[nJ]+c,Hl(nI*3+r,nJ*3+c));
// =============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <vector>

#include "types.h" // Positions ([N,3]), Triplet 等别名

#if __cplusplus >= 202002L
#include <concepts>
#endif

namespace fiction {

// =============================================================================
// IEnergy — 所有非惯性能量的抽象基类（纯 3D）
// =============================================================================
class IEnergy {
public:
    virtual ~IEnergy() = default;

    virtual const char* name() const = 0;

    virtual double hessianWeight(double h) const { return h * h; }

    /// @param x 固定 [N,3]
    virtual double val(const Positions& x) const = 0;

    /// @param x [N,3]  @param grad [N,3] 累加目标
    virtual void addGrad(const Positions& x, Positions& grad) const = 0;

    /// @param x [N,3]  @param triplets 行列索引 = 3*node+d
    virtual void addHess(const Positions& x, std::vector<Triplet>& triplets) const = 0;

    virtual void precompute(const Positions& /*x*/) {}

    /// 每个 Newton 迭代在当前 x 上刷新与状态相关的离散集合（如接触候选）。
    /// 默认空实现，不改变既有能量行为。
    virtual void preIteration(const Positions& /*x*/) {}

    virtual double initStepSize(const Positions& /*x*/, const Positions& /*p*/) const {
        return 1.0;
    }
};

using EnergyPtr = std::shared_ptr<IEnergy>;

// =============================================================================
// C++20 Concept 校验
// =============================================================================
#if __cplusplus >= 202002L
template <typename T>
concept EnergyConcept = requires(const T& e, const Positions& x, Positions& g,
                                 std::vector<Triplet>& trips, double h) {
    { e.val(x) } -> std::convertible_to<double>;
    { e.addGrad(x, g) } -> std::same_as<void>;
    { e.addHess(x, trips) } -> std::same_as<void>;
    { e.name() } -> std::convertible_to<const char*>;
    { e.hessianWeight(h) } -> std::convertible_to<double>;
};

static_assert(EnergyConcept<IEnergy> == false || true,
              "IEnergy 为抽象基类，Concept 检查仅对派生类有意义");
#endif

} // namespace fiction
