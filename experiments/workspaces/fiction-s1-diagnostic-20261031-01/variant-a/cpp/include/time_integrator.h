#pragma once
// =============================================================================
// time_integrator.h — 隐式欧拉 + Newton 数值迭代主框架（纯 3D 版）
// =============================================================================
// 职责（对标 Python_slovement/time_integrator.py:17 step_forward）：
//   - 纯 3D 主线：Positions 固定 [N,3]，kDim=3，索引 3*node+d
//   - 2D 适配：调用方先用 upgrade2DTo3D + lockZ 再进入 stepForward
//   - 其余沿用原设计：Inertia 内置，IEnergy 扩展，DBC 罚 + 自由度消去
// =============================================================================

#include "types.h"
#include "energy.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <memory>
#include <string>
#include <vector>

namespace fiction {

class TimeIntegrator {
public:
    struct Config {
        double h = 0.01;      ///< 时间步长，对应 Python h
        double tol = 1e-2;    ///< Newton 收敛容差，对应 Python tol
        int maxIter = 100;    ///< 最大 Newton 迭代数
        bool verbose = true;  ///< 是否打印迭代日志
        bool staticSolve = false; ///< 静态平衡：不含惯性项，能量权重为 1
    };

    struct StepReport {
        bool converged = false;
        bool maxIterHit = false;
        bool lineSearchFailed = false;
        int iterations = 0;
        int lineSearchIterations = 0;
        double searchDirectionRate = 0.0;
        double freeGradientInf = 0.0;
        double dbcMaxError = 0.0;
        double dbcStiffness = 0.0;
        int dbcSatisfied = 0;
        int dbcTotal = 0;
        std::string solver;
    };

    struct State {
        Positions x;     ///< [N,3] 更新后位置
        Velocities v;    ///< [N,3] 更新后速度
        StepReport report;
    };

    /// @param cfg    时间积分配置
    /// @param masses 质量向量，size N，对应 Python m
    explicit TimeIntegrator(Config cfg, const Eigen::VectorXd& masses);

    void addEnergy(EnergyPtr e) { energies_.push_back(std::move(e)); }
    void clearEnergies() { energies_.clear(); }
    std::size_t numEnergies() const { return energies_.size(); }

    /// 输入 x [N,3], v [N,3], DBCInfo（Vec3 target/limit）
    State stepForward(const Positions& x, const Velocities& v, const DBCInfo& dbc);

private:
    Config cfg_;
    Eigen::VectorXd masses_;                 ///< size N
    std::vector<EnergyPtr> energies_;        ///< 除 Inertia 外的所有能量

    // --- IP 组装，对标 Python IP_val/IP_grad/IP_hess（纯 3D） ---
    double computeIPVal(const Positions& x,
                        const Positions& x_tilde,
                        const Positions& x_n,
                        const DBCInfo& dbc,
                        const std::vector<Vec3>& dbcTargets) const;

    Positions computeIPGrad(const Positions& x,
                            const Positions& x_tilde,
                            const Positions& x_n,
                            const DBCInfo& dbc,
                            const std::vector<Vec3>& dbcTargets) const;

    SparseMatrix computeIPHess(const Positions& x,
                               const Positions& x_tilde,
                               const Positions& x_n,
                               const DBCInfo& dbc,
                               const std::vector<Vec3>& dbcTargets) const;

    // --- Inertia 内置项（纯 3D） ---
    double inertiaVal(const Positions& x, const Positions& x_tilde) const;
    void inertiaAddGrad(const Positions& x, const Positions& x_tilde, Positions& grad) const;
    void inertiaAddHess(std::vector<Triplet>& trips) const;

    // --- DBC 罚能量，对标 SpringEnergy.py（纯 3D） ---
    // val 0.5*k*m*|x-target|^2, grad k*m*diff, hess diag k*m
    double dbcVal(const Positions& x,
                  const DBCInfo& dbc,
                  const std::vector<Vec3>& targets) const;
    void dbcAddGrad(const Positions& x,
                    const DBCInfo& dbc,
                    const std::vector<Vec3>& targets,
                    Positions& grad) const;
    void dbcAddHess(const DBCInfo& dbc, std::vector<Triplet>& trips) const;

    // --- 搜索方向，对标 search_dir:96（纯 3D） ---
    struct SearchResult {
        Positions p;                      ///< [N,3] Newton 增量
        std::vector<char> satisfied;      ///< size N，标记 DBC 是否已满足
        bool linearSolveFailed = false;
        std::string solver;
    };
    SearchResult computeSearchDirection(const Positions& x,
                                        const Positions& x_tilde,
                                        const Positions& x_n,
                                        const DBCInfo& dbc,
                                        const std::vector<Vec3>& dbcTargets) const;

    // --- 线性求解：Eigen 组装 + CHOLMOD 主线，失败回退 SimplicialLDLT ---
    struct LinearSolveResult { Eigen::VectorXd x; bool failed = false; std::string solver; };
    LinearSolveResult solveLinearSystem(const SparseMatrix& H,
                                        const Eigen::VectorXd& rhs) const;
    double effectiveH() const { return cfg_.staticSolve ? 1.0 : cfg_.h; }
};

} // namespace fiction
