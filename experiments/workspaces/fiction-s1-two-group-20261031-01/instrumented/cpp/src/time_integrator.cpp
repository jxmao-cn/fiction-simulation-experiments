// =============================================================================
// time_integrator.cpp — 隐式欧拉 Newton 迭代实现（纯 3D 版）
// =============================================================================
// 设计约束（纯 3D 主线）：
//   - Positions 固定 [N,3] RowMajor，kDim=3，所有自由度索引 = 3*node+d
//   - 2D 适配：upgrade2DTo3D + lockZAxis 在调用方完成，本文件不再做 Dim 推断
//   - Inertia 内置，其余能量统一走 IEnergy
//   - 保持 Eigen 封装 + CHOLMOD 主线，C++20
// =============================================================================

#include "time_integrator.h"

#include <Eigen/SparseCholesky> // SimplicialLDLT
// Eigen CholmodSupport expects <cholmod.h> at include root; on macOS brew
// the header lives at <suitesparse/cholmod.h>. Support both layouts.
#if defined(FICTION_USE_CHOLMOD) && __has_include(<Eigen/CholmodSupport>)
#  if __has_include(<cholmod.h>)
#    include <Eigen/CholmodSupport>
#    define FICTION_HAS_CHOLMOD 1
#  elif __has_include(<suitesparse/cholmod.h>)
#    include <Eigen/CholmodSupport>
#    define FICTION_HAS_CHOLMOD 1
#  endif
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace fiction {

// -----------------------------------------------------------------------------
// 工具：DBC 相关辅助（纯 3D）
// -----------------------------------------------------------------------------
namespace {

inline bool isDBCNode(const DBCInfo& dbc, int node, int N) {
    if (node < 0 || node >= N) return false;
    if (!dbc.isDBC.empty()) {
        if ((int)dbc.isDBC.size() == N) return dbc.isDBC[node] != 0;
        if ((int)dbc.isDBC.size() > node) return dbc.isDBC[node] != 0;
    }
    for (int idx : dbc.indices) if (idx == node) return true;
    return false;
}

// fetchLimit/fetchVel 在纯 3D 下：limits/targetVelocities 与 indices 一一对应
inline Vec3 fetchLimit(const DBCInfo& dbc, size_t k, int idx, int N) {
    if (!dbc.limits.empty()) {
        if (dbc.limits.size() == dbc.indices.size() && k < dbc.limits.size())
            return dbc.limits[k];
        if ((int)dbc.limits.size() == N && idx >= 0 && idx < (int)dbc.limits.size())
            return dbc.limits[idx];
        if (!dbc.limits.empty()) return dbc.limits[0];
    }
    return Vec3::Zero();
}

inline Vec3 fetchVel(const DBCInfo& dbc, size_t k, int idx, int N) {
    if (!dbc.targetVelocities.empty()) {
        if (dbc.targetVelocities.size() == dbc.indices.size() && k < dbc.targetVelocities.size())
            return dbc.targetVelocities[k];
        if ((int)dbc.targetVelocities.size() == N && idx >= 0 && idx < (int)dbc.targetVelocities.size())
            return dbc.targetVelocities[idx];
        if (!dbc.targetVelocities.empty()) return dbc.targetVelocities[0];
    }
    return Vec3::Zero();
}

} // anonymous namespace

// =============================================================================
// 构造
// =============================================================================
TimeIntegrator::TimeIntegrator(Config cfg, const Eigen::VectorXd& masses)
    : cfg_(cfg), masses_(masses) {
    if (!(cfg_.h > 0.0) || !std::isfinite(cfg_.h) || !(cfg_.tol > 0.0) ||
        !std::isfinite(cfg_.tol) || cfg_.maxIter < 0)
        throw std::invalid_argument("TimeIntegrator: invalid h/tol/maxIter");
    if (masses_.size() == 0 || !masses_.allFinite() || (masses_.array() <= 0.0).any())
        throw std::invalid_argument("TimeIntegrator: masses must be finite and positive");
}

// =============================================================================
// Inertia 内置项（纯 3D）
// =============================================================================
double TimeIntegrator::inertiaVal(const Positions& x, const Positions& x_tilde) const {
    const int N = static_cast<int>(masses_.size());
    double E = 0.0;
    for (int i = 0; i < N; ++i) {
        Vec3 diff = x.row(i).transpose() - x_tilde.row(i).transpose();
        E += 0.5 * masses_[i] * diff.squaredNorm();
    }
    return E;
}

void TimeIntegrator::inertiaAddGrad(const Positions& x, const Positions& x_tilde, Positions& grad) const {
    const int N = static_cast<int>(masses_.size());
    for (int i = 0; i < N; ++i) {
        grad.row(i) += masses_[i] * (x.row(i) - x_tilde.row(i));
    }
}

void TimeIntegrator::inertiaAddHess(std::vector<Triplet>& trips) const {
    const int N = static_cast<int>(masses_.size());
    trips.reserve(trips.size() + N * kDim);
    for (int i = 0; i < N; ++i) {
        const double m = masses_[i];
        for (int d = 0; d < kDim; ++d) {
            const int idx = kDim * i + d;
            trips.emplace_back(idx, idx, m);
        }
    }
}

// =============================================================================
// DBC 罚能量（SpringEnergy.py 对标，纯 3D）
//   val  0.5*k*m*|x-target|^2
//   grad k*m*(x-target)
//   hess diag k*m
// =============================================================================
double TimeIntegrator::dbcVal(const Positions& x,
                              const DBCInfo& dbc,
                              const std::vector<Vec3>& targets) const {
    double E = 0.0;
    for (size_t k = 0; k < dbc.indices.size(); ++k) {
        const int idx = dbc.indices[k];
        if (idx < 0 || idx >= (int)masses_.size()) continue;
        if (k >= targets.size()) continue;
        Vec3 diff = x.row(idx).transpose() - targets[k];
        E += 0.5 * dbc.stiffness * masses_[idx] * diff.squaredNorm();
    }
    return E;
}

void TimeIntegrator::dbcAddGrad(const Positions& x,
                                const DBCInfo& dbc,
                                const std::vector<Vec3>& targets,
                                Positions& grad) const {
    for (size_t k = 0; k < dbc.indices.size(); ++k) {
        const int idx = dbc.indices[k];
        if (idx < 0 || idx >= (int)masses_.size()) continue;
        if (k >= targets.size()) continue;
        Vec3 diff = x.row(idx).transpose() - targets[k];
        Vec3 g = dbc.stiffness * masses_[idx] * diff;
        grad.row(idx) += g.transpose();
    }
}

void TimeIntegrator::dbcAddHess(const DBCInfo& dbc, std::vector<Triplet>& trips) const {
    trips.reserve(trips.size() + (int)dbc.indices.size() * kDim);
    for (size_t k = 0; k < dbc.indices.size(); ++k) {
        const int idx = dbc.indices[k];
        if (idx < 0 || idx >= (int)masses_.size()) continue;
        const double km = dbc.stiffness * masses_[idx];
        for (int d = 0; d < kDim; ++d) {
            const int j = kDim * idx + d;
            trips.emplace_back(j, j, km);
        }
    }
}

// =============================================================================
// IP 组装：对标 IP_val / IP_grad / IP_hess（纯 3D）
// =============================================================================
double TimeIntegrator::computeIPVal(const Positions& x,
                                    const Positions& x_tilde,
                                    const Positions& x_n,
                                    const DBCInfo& dbc,
                                    const std::vector<Vec3>& dbcTargets) const {
    (void)x_n;
    double E = cfg_.staticSolve ? 0.0 : inertiaVal(x, x_tilde);
    for (const auto& e : energies_) {
        const double w = e->hessianWeight(effectiveH());
        E += w * e->val(x);
    }
    E += dbcVal(x, dbc, dbcTargets);
    return E;
}

Positions TimeIntegrator::computeIPGrad(const Positions& x,
                                        const Positions& x_tilde,
                                        const Positions& x_n,
                                        const DBCInfo& dbc,
                                        const std::vector<Vec3>& dbcTargets) const {
    (void)x_n;
    const int N = static_cast<int>(masses_.size());
    Positions grad = Positions::Zero(N, kDim);
    if (!cfg_.staticSolve) inertiaAddGrad(x, x_tilde, grad);
    for (const auto& e : energies_) {
        Positions contrib = Positions::Zero(N, kDim);
        e->addGrad(x, contrib);
        if (contrib.rows() != N || contrib.cols() != kDim) {
            Positions tmp = Positions::Zero(N, kDim);
            const int rn = std::min<int>(N, (int)contrib.rows());
            const int rd = std::min<int>(kDim, (int)contrib.cols());
            tmp.block(0,0,rn,rd) = contrib.block(0,0,rn,rd);
            contrib = tmp;
        }
        const double w = e->hessianWeight(effectiveH());
        if (w != 1.0) contrib *= w;
        grad += contrib;
    }
    dbcAddGrad(x, dbc, dbcTargets, grad);
    return grad;
}

SparseMatrix TimeIntegrator::computeIPHess(const Positions& x,
                                           const Positions& x_tilde,
                                           const Positions& x_n,
                                           const DBCInfo& dbc,
                                           const std::vector<Vec3>& dbcTargets) const {
    (void)x_tilde; (void)x_n; (void)dbcTargets;
    const int N = static_cast<int>(masses_.size());
    const int dof = N * kDim;
    std::vector<Triplet> trips;
    trips.reserve(dof * 4 + 4 * (int)energies_.size() * kDim * 2);

    if (!cfg_.staticSolve) inertiaAddHess(trips);

    for (const auto& e : energies_) {
        const size_t before = trips.size();
        e->addHess(x, trips);
        const double w = e->hessianWeight(effectiveH());
        if (w != 1.0) {
            for (size_t i = before; i < trips.size(); ++i) {
                const auto& t = trips[i];
                trips[i] = Triplet(t.row(), t.col(), t.value() * w);
            }
        }
    }

    dbcAddHess(dbc, trips);

    SparseMatrix H(dof, dof);
    H.setFromTriplets(trips.begin(), trips.end());
    return H;
}

// =============================================================================
// 线性求解：Eigen 封装 + CHOLMOD 主线，回退 SimplicialLDLT
// =============================================================================
TimeIntegrator::LinearSolveResult TimeIntegrator::solveLinearSystem(
    const SparseMatrix& H, const Eigen::VectorXd& rhs) const {
    LinearSolveResult result;
#ifdef FICTION_USE_CHOLMOD
    Eigen::CholmodSupernodalLLT<SparseMatrix> chol;
    chol.compute(H);
    if (chol.info() == Eigen::Success) {
        result.x = chol.solve(rhs);
        if (chol.info() == Eigen::Success && result.x.allFinite()) {
            result.solver = "CHOLMOD";
            return result;
        }
    }
#endif
    // FICTION_USE_CHOLMOD=OFF 时此处确实只走 Eigen LDLT。
    Eigen::SimplicialLDLT<SparseMatrix> ldlt;
    ldlt.compute(H);
    if (ldlt.info() == Eigen::Success) {
        result.x = ldlt.solve(rhs);
        if (ldlt.info() == Eigen::Success && result.x.allFinite()) {
            result.solver =
#ifdef FICTION_USE_CHOLMOD
                "SimplicialLDLT(fallback)";
#else
                "SimplicialLDLT";
#endif
            return result;
        }
    }
    // 保留旧实现针对数值奇异系统的微小对角正则化兜底，并显式报告实际路径。
    std::vector<Triplet> regularizedTrips;
    regularizedTrips.reserve(H.nonZeros() + H.rows());
    for (int col = 0; col < H.outerSize(); ++col)
        for (SparseMatrix::InnerIterator it(H, col); it; ++it)
            regularizedTrips.emplace_back(it.row(), it.col(), it.value());
    for (int i = 0; i < H.rows(); ++i) regularizedTrips.emplace_back(i, i, 1e-10);
    SparseMatrix Hreg(H.rows(), H.cols());
    Hreg.setFromTriplets(regularizedTrips.begin(), regularizedTrips.end());
    Eigen::SimplicialLDLT<SparseMatrix> regSolver;
    regSolver.compute(Hreg);
    if (regSolver.info() == Eigen::Success) {
        result.x = regSolver.solve(rhs);
        if (regSolver.info() == Eigen::Success && result.x.allFinite()) {
            result.solver = "SimplicialLDLT(regularized)";
            return result;
        }
    }
    result.failed = true;
    result.solver = "failed";
    result.x = Eigen::VectorXd::Constant(rhs.size(), std::numeric_limits<double>::quiet_NaN());
    return result;
}

// =============================================================================
// 搜索方向：对标 search_dir:96-114（纯 3D）
// =============================================================================
TimeIntegrator::SearchResult TimeIntegrator::computeSearchDirection(
    const Positions& x,
    const Positions& x_tilde,
    const Positions& x_n,
    const DBCInfo& dbc,
    const std::vector<Vec3>& dbcTargets) const {

    const int N = static_cast<int>(masses_.size());
    const int dof = N * kDim;

    SparseMatrix H = computeIPHess(x, x_tilde, x_n, dbc, dbcTargets);
    Positions gradMat = computeIPGrad(x, x_tilde, x_n, dbc, dbcTargets);
    FlatVector g = flatten(gradMat); // [3*N]

    std::vector<char> satisfied(N, 0);
    for (size_t k = 0; k < dbc.indices.size(); ++k) {
        const int idx = dbc.indices[k];
        if (idx < 0 || idx >= N) continue;
        if (k >= dbcTargets.size()) continue;
        Vec3 diff = x.row(idx).transpose() - dbcTargets[k];
        if (diff.norm() / effectiveH() < cfg_.tol) satisfied[idx] = 1;
    }

    FlatVector gMod = g;
    for (int i = 0; i < N; ++i) {
        if (isDBCNode(dbc, i, N) && satisfied[i]) {
            for (int d = 0; d < kDim; ++d) gMod[kDim * i + d] = 0.0;
        }
    }

    std::vector<Triplet> newTrips;
    newTrips.reserve(H.nonZeros() + N * kDim);
    std::vector<char> diagAdded(dof, 0);

    for (int col = 0; col < H.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(H, col); it; ++it) {
            const int row = it.row();
            const int c = it.col();
            const bool rowSat = isDBCNode(dbc, row / kDim, N) && satisfied[row / kDim];
            const bool colSat = isDBCNode(dbc, c / kDim, N) && satisfied[c / kDim];
            if (rowSat || colSat) {
                if (row == c) {
                    newTrips.emplace_back(row, c, 1.0);
                    diagAdded[row] = 1;
                }
            } else {
                newTrips.emplace_back(row, c, it.value());
            }
        }
    }
    for (int i = 0; i < N; ++i) {
        if (isDBCNode(dbc, i, N) && satisfied[i]) {
            for (int d = 0; d < kDim; ++d) {
                const int idx = kDim * i + d;
                if (!diagAdded[idx]) newTrips.emplace_back(idx, idx, 1.0);
            }
        }
    }

    SparseMatrix Hmod(dof, dof);
    Hmod.setFromTriplets(newTrips.begin(), newTrips.end());

    FlatVector rhs = -gMod;
    LinearSolveResult solve = solveLinearSystem(Hmod, rhs);
    Positions p = solve.failed ? Positions::Constant(N, kDim,
        std::numeric_limits<double>::quiet_NaN()) : unflatten(solve.x);
    return {p, satisfied, solve.failed, solve.solver};
}

// =============================================================================
// 主循环：对标 Python step_forward:17-57（纯 3D）
// =============================================================================
TimeIntegrator::State TimeIntegrator::stepForward(const Positions& x_in,
                                                  const Velocities& v_in,
                                                  const DBCInfo& dbcIn) {
    const int N = static_cast<int>(masses_.size());
    if (x_in.rows() != N || v_in.rows() != N || x_in.cols() != kDim || v_in.cols() != kDim)
        throw std::invalid_argument("TimeIntegrator::stepForward dimension mismatch");
    if (!x_in.allFinite() || !v_in.allFinite())
        throw std::invalid_argument("TimeIntegrator::stepForward input contains nonfinite values");
    DBCInfo dbc = dbcIn;
    if (!std::isfinite(dbc.stiffness) || dbc.stiffness < 0.0)
        throw std::invalid_argument("TimeIntegrator::stepForward invalid DBC stiffness");

    Positions x_n = x_in;
    Positions x = x_in;
    Positions x_tilde = cfg_.staticSolve ? x_in : x_in + v_in * cfg_.h;

    for (auto& e : energies_) e->precompute(x_n);

    // DBC_target 计算（纯 3D）
    std::vector<Vec3> dbcTargets;
    dbcTargets.reserve(dbc.indices.size());
    const bool hasAbsTargets = static_cast<int>(dbc.targetPositions.size()) ==
                               static_cast<int>(dbc.indices.size());
    for (size_t k = 0; k < dbc.indices.size(); ++k) {
        const int idx = dbc.indices[k];
        if (idx < 0 || idx >= N) {
            dbcTargets.emplace_back(Vec3::Zero());
            continue;
        }
        if (hasAbsTargets) {
            // 绝对目标位置优先（如 --twist 旋转轨迹驱动）：调用方逐帧更新
            dbcTargets.push_back(dbc.targetPositions[k]);
            continue;
        }
        Vec3 xn = x_n.row(idx).transpose();
        Vec3 lim = fetchLimit(dbc, k, idx, N);
        Vec3 vel = fetchVel(dbc, k, idx, N);
        Vec3 tgt;
        if ((lim - xn).dot(vel) > 0) tgt = xn + effectiveH() * vel;
        else tgt = xn;
        dbcTargets.push_back(tgt);
    }

    double E_last = computeIPVal(x, x_tilde, x_n, dbc, dbcTargets);
    SearchResult sr = computeSearchDirection(x, x_tilde, x_n, dbc, dbcTargets);
    Positions p = sr.p;
    std::vector<char> satisfied = sr.satisfied;
    StepReport report;
    report.dbcTotal = static_cast<int>(dbc.indices.size());
    report.dbcStiffness = dbc.stiffness;
    report.solver = sr.solver;
    if (!std::isfinite(E_last) || !p.allFinite()) report.lineSearchFailed = !sr.linearSolveFailed;

    auto infNormPos = [](const Positions& pp) -> double {
        if (pp.size() == 0) return 0.0;
        FlatVector f = flatten(pp);
        return f.lpNorm<Eigen::Infinity>();
    };

    // 可选纯观测轨迹：诊断量从当前状态重新计算，不参与迭代决策。
    std::vector<std::string> traceRows;
    bool stiffnessChanged = false;
    auto captureTrace = [&](int stateIter, const char* stage, double usedAlpha,
                            int halvings, double dxInf, bool trialEqualsX) {
        if (cfg_.traceCsvPath.empty()) return;
        Positions physicalGrad = computeIPGrad(x, x_tilde, x_n, dbc, dbcTargets);
        Positions penaltyGrad = Positions::Zero(N, kDim);
        dbcAddGrad(x, dbc, dbcTargets, penaltyGrad);
        physicalGrad -= penaltyGrad;
        double matInf = 0.0, penInf = 0.0, totalInf = 0.0;
        double freeInf = 0.0, dbcError = 0.0;
        double balNum2 = 0.0, matNorm2 = 0.0, penNorm2 = 0.0;
        int unsatCount = 0;
        for (int i = 0; i < N; ++i) {
            if (!isDBCNode(dbc, i, N)) {
                freeInf = std::max(freeInf, physicalGrad.row(i).cwiseAbs().maxCoeff());
                continue;
            }
            if (i < (int)satisfied.size() && satisfied[i]) continue;
            ++unsatCount;
            const Vec3 gm = physicalGrad.row(i).transpose();
            const Vec3 gp = penaltyGrad.row(i).transpose();
            const Vec3 gt = gm + gp;
            matInf = std::max(matInf, gm.cwiseAbs().maxCoeff());
            penInf = std::max(penInf, gp.cwiseAbs().maxCoeff());
            totalInf = std::max(totalInf, gt.cwiseAbs().maxCoeff());
            balNum2 += gt.squaredNorm(); matNorm2 += gm.squaredNorm(); penNorm2 += gp.squaredNorm();
        }
        for (size_t k = 0; k < dbc.indices.size(); ++k) {
            const int idx = dbc.indices[k];
            if (idx >= 0 && idx < N && k < dbcTargets.size())
                dbcError = std::max(dbcError, (x.row(idx).transpose() - dbcTargets[k]).norm());
        }
        const double balance = unsatCount == 0 ? 0.0 : std::sqrt(balNum2) /
            std::max({std::sqrt(matNorm2), std::sqrt(penNorm2), 1e-30});
        double materialEnergy = cfg_.staticSolve ? 0.0 : inertiaVal(x, x_tilde);
        for (const auto& e : energies_) materialEnergy += e->hessianWeight(effectiveH()) * e->val(x);
        const double penaltyEnergy = dbcVal(x, dbc, dbcTargets);
        const Positions objectiveGrad = physicalGrad + penaltyGrad;
        const double gradDotP = objectiveGrad.size() == 0 ? 0.0 :
            flatten(objectiveGrad).dot(flatten(p));
        const int satCount = (int)std::count(satisfied.begin(), satisfied.end(), 1);
        std::ostringstream row;
        row << std::setprecision(17) << cfg_.traceLabel << ',' << stage << ',' << stateIter << ','
            << dbc.stiffness << ',' << (p.allFinite() ? infNormPos(p) / effectiveH() : std::numeric_limits<double>::max()) << ','
            << satCount << ',' << report.dbcTotal << ',' << dbcError << ','
            << freeInf << ',' << matInf << ',' << penInf << ',' << totalInf << ','
            << balance << ',' << unsatCount << ',' << usedAlpha << ',' << halvings << ',' << dxInf << ','
            << (trialEqualsX ? 1 : 0) << ',' << (E_last + 0.0) << ',' << materialEnergy << ','
            << penaltyEnergy << ',' << report.dbcStiffness << ',' << (stiffnessChanged ? 1 : 0) << ',' << gradDotP;
        traceRows.push_back(row.str());
    };
    captureTrace(0, "initial", std::numeric_limits<double>::quiet_NaN(), 0, 0.0, false);

    int iter = 0;
    while (!sr.linearSolveFailed && p.allFinite() && std::isfinite(E_last) &&
           (infNormPos(p) / effectiveH() > cfg_.tol ||
            std::count(satisfied.begin(), satisfied.end(), 1) != report.dbcTotal)) {
        if (iter >= cfg_.maxIter) { report.maxIterHit = true; break; }

        const double res = infNormPos(p) / effectiveH();
        if (cfg_.verbose) {
            std::cout << "Iteration " << iter << ":\n";
            std::cout << "residual = " << res << "\n";
        }

        if (res <= cfg_.tol &&
            std::count(satisfied.begin(), satisfied.end(), 1) != (int)dbc.indices.size()) {
            dbc.stiffness *= 2.0;
            stiffnessChanged = true;
            if (cfg_.verbose) std::cout << "DBC not satisfied, doubling stiffness to " << dbc.stiffness << "\n";
            E_last = computeIPVal(x, x_tilde, x_n, dbc, dbcTargets);
            sr = computeSearchDirection(x, x_tilde, x_n, dbc, dbcTargets);
            p = sr.p;
            satisfied = sr.satisfied;
            report.dbcStiffness = dbc.stiffness;
            report.solver = sr.solver;
            if (sr.linearSolveFailed || !p.allFinite() || !std::isfinite(E_last)) {
                report.lineSearchFailed = !sr.linearSolveFailed;
                break;
            }
        }

        double alpha = 1.0;
        for (const auto& e : energies_) {
            const double a = e->initStepSize(x, p);
            if (a < alpha) alpha = a;
        }
        if (alpha > 1.0) alpha = 1.0;
        if (alpha <= 0.0) alpha = 1e-12;

        int lsIter = 0;
        const int maxLS = 60;
        bool accepted = false;
        while (lsIter < maxLS) {
            Positions xTrial = x + alpha * p;
            const double E_trial = computeIPVal(xTrial, x_tilde, x_n, dbc, dbcTargets);
            if (std::isfinite(E_trial) && E_trial <= E_last) { accepted = true; break; }
            alpha *= 0.5;
            ++lsIter;
            if (alpha < 1e-14) break;
        }
        report.lineSearchIterations += lsIter;
        if (!accepted) { report.lineSearchFailed = true; break; }

        if (cfg_.verbose) {
            std::cout << "step size = " << alpha << "\n";
            if (lsIter >= maxLS) std::cout << "[Warning] line search reached maxLS\n";
            if (alpha < 1e-12) std::cout << "[Warning] alpha very small: " << alpha << "\n";
        }

        const Positions oldX = x;
        x += alpha * p;
        const double actualDxInf = infNormPos(x - oldX);
        const bool trialEqualsX = (x.array() == oldX.array()).all();
        E_last = computeIPVal(x, x_tilde, x_n, dbc, dbcTargets);

        sr = computeSearchDirection(x, x_tilde, x_n, dbc, dbcTargets);
        p = sr.p;
        satisfied = sr.satisfied;

        ++iter;
        report.iterations = iter;
        report.dbcStiffness = dbc.stiffness;
        report.solver = sr.solver;
        captureTrace(iter, "accepted", alpha, lsIter, actualDxInf, trialEqualsX);
        stiffnessChanged = false;
    }

    report.iterations = iter;
    report.searchDirectionRate = p.allFinite() ? infNormPos(p) / effectiveH()
        : std::numeric_limits<double>::max();
    report.dbcSatisfied = static_cast<int>(std::count(satisfied.begin(), satisfied.end(), 1));
    report.dbcMaxError = 0.0;
    for (size_t k = 0; k < dbc.indices.size(); ++k) {
        const int idx = dbc.indices[k];
        if (idx >= 0 && idx < N && k < dbcTargets.size())
            report.dbcMaxError = std::max(report.dbcMaxError,
                std::isfinite((x.row(idx).transpose() - dbcTargets[k]).norm())
                    ? (x.row(idx).transpose() - dbcTargets[k]).norm()
                    : std::numeric_limits<double>::max());
    }
    Positions finalGrad = computeIPGrad(x, x_tilde, x_n, dbc, dbcTargets);
    // 自由梯度诊断排除 DBC 罚项，只统计非 DBC 节点。
    Positions dbcGrad = Positions::Zero(N, kDim);
    dbcAddGrad(x, dbc, dbcTargets, dbcGrad);
    finalGrad -= dbcGrad;
    const bool finalGradientFinite = finalGrad.allFinite();
    if (!finalGradientFinite) report.lineSearchFailed = true;
    report.freeGradientInf = 0.0;
    for (int i = 0; i < N; ++i) if (!isDBCNode(dbc, i, N))
        report.freeGradientInf = std::max(report.freeGradientInf,
            finalGrad.row(i).allFinite() ? finalGrad.row(i).cwiseAbs().maxCoeff()
                                        : std::numeric_limits<double>::max());
    if (finalGrad.rows() > 0 && !finalGradientFinite)
        report.freeGradientInf = std::numeric_limits<double>::max();
    Velocities v_new(N, kDim);
    if (cfg_.staticSolve) v_new.setZero();
    else v_new = (x - x_n) / cfg_.h;
    report.converged = !sr.linearSolveFailed && !report.lineSearchFailed &&
        !report.maxIterHit && p.allFinite() && std::isfinite(report.freeGradientInf) &&
        std::isfinite(E_last) && x.allFinite() && v_new.allFinite() && finalGradientFinite &&
        report.searchDirectionRate <= cfg_.tol && report.dbcSatisfied == report.dbcTotal;

    if (!cfg_.traceCsvPath.empty()) {
        if (!report.converged && !traceRows.empty()) {
            const auto firstComma = traceRows.back().find(',');
            const auto secondComma = traceRows.back().find(',', firstComma + 1);
            if (firstComma != std::string::npos && secondComma != std::string::npos)
                traceRows.back().replace(firstComma + 1, secondComma - firstComma - 1, "terminal");
        }
        std::ofstream trace(cfg_.traceCsvPath);
        if (!trace) throw std::runtime_error("cannot open trace CSV: " + cfg_.traceCsvPath);
        trace << "label,stage,iter,k,search_rate,dbc_satisfied,dbc_total,dbc_error_max,free_gradient_inf,"
                 "material_force_inf_unsatisfied,penalty_force_inf_unsatisfied,total_force_inf_unsatisfied,"
                 "force_balance_relative_unsatisfied,unsatisfied_nodes,alpha,halvings,actual_dx_inf,"
                 "trial_equals_x,energy,material_energy,penalty_energy,dbc_stiffness,stiffness_changed,"
                 "state_gradient_dot_current_direction,final_converged,final_max_iter_hit,final_line_search_failed\n";
        for (const auto& row : traceRows)
            trace << row << ',' << (report.converged ? 1 : 0) << ',' << (report.maxIterHit ? 1 : 0)
                  << ',' << (report.lineSearchFailed ? 1 : 0) << '\n';
        if (!trace) throw std::runtime_error("failed writing trace CSV: " + cfg_.traceCsvPath);

        const std::string nodePath = cfg_.traceCsvPath + ".nodes.csv";
        std::ofstream nodes(nodePath);
        if (!nodes) throw std::runtime_error("cannot open node trace CSV: " + nodePath);
        nodes << "node,is_dbc,satisfied,xerror_x,xerror_y,xerror_z,material_force_x,material_force_y,material_force_z,"
                 "penalty_force_x,penalty_force_y,penalty_force_z,total_force_x,total_force_y,total_force_z\n";
        Positions nodePhysical = computeIPGrad(x, x_tilde, x_n, dbc, dbcTargets);
        Positions nodePenalty = Positions::Zero(N, kDim);
        dbcAddGrad(x, dbc, dbcTargets, nodePenalty);
        nodePhysical -= nodePenalty;
        for (int i = 0; i < N; ++i) {
            Vec3 error = Vec3::Zero();
            for (size_t k = 0; k < dbc.indices.size(); ++k)
                if (dbc.indices[k] == i && k < dbcTargets.size()) { error = x.row(i).transpose() - dbcTargets[k]; break; }
            const Vec3 fm = -nodePhysical.row(i).transpose();
            const Vec3 fp = -nodePenalty.row(i).transpose();
            const Vec3 ft = fm + fp;
            nodes << i << ',' << (isDBCNode(dbc, i, N) ? 1 : 0) << ','
                  << (isDBCNode(dbc, i, N) && i < (int)satisfied.size() && satisfied[i] ? 1 : 0)
                  << ',' << error.x() << ',' << error.y() << ',' << error.z()
                  << ',' << fm.x() << ',' << fm.y() << ',' << fm.z()
                  << ',' << fp.x() << ',' << fp.y() << ',' << fp.z()
                  << ',' << ft.x() << ',' << ft.y() << ',' << ft.z() << '\n';
        }
        if (!nodes) throw std::runtime_error("failed writing node trace CSV: " + nodePath);
    }

    if (cfg_.verbose) {
        const double finalRes = infNormPos(p) / effectiveH();
        const int satCnt = std::count(satisfied.begin(), satisfied.end(), 1);
        std::cout << "Newton finished in " << iter << " iterations, final residual=" << finalRes
                  << " satisfied " << satCnt << "/" << dbc.indices.size()
                  << " stiffness=" << dbc.stiffness << " Dim=" << kDim << "\n";
    }

    return {x, v_new, report};
}

} // namespace fiction
