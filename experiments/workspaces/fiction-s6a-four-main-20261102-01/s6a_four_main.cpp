// =============================================================================
// examples/s6a_multibody_check.cpp — S6a 无接触多体全局场景检查（1 体 / 4 体）
// =============================================================================
// 预登记见 docs/s3-s6a-execution-design-20261031.md 第 5 节：
//   - 同一份 sphere1K 资产按纯平移复制成 1 / 4 个体，共用一套全局状态与时间步入口
//   - NH E=1e5 Pa、nu=0.3、rho=1000、重力 (0,0,-9.81)、dt=0.005、10 步
//   - tol=1e-8、maxIter=100、DBC 初始 k=1e8、每体固定局部 z>=0.45 的顶部 cap
//   - 检查：全局 Hessian 按 body_id 无跨体条目；4 体轨迹 = 4 次独立单体求解（<1e-8 m）
// 本程序比较无接触多体与独立单体轨迹。
// =============================================================================
#include "energies/Gravity.h"
#include "energies/NeoHookeanTetEnergy.h"
#include "io/msh_loader.h"
#include "scene/mass_spring_builder.h"
#include "scene/multibody_scene.h"
#include "time_integrator.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace fiction;

constexpr double kYoungs = 1.0e5;
constexpr double kPoisson = 0.3;
constexpr double kDensity = 1000.0;
constexpr double kGravityZ = -9.81;
constexpr double kDt = 0.005;
constexpr int kSteps = 10;
constexpr double kTol = 1e-8;
constexpr int kMaxIter = 100;
constexpr double kDbcK = 1e8;
constexpr double kCapZMin = 0.45;      ///< 局部 z >= 0.45 的顶部 cap 被固定
constexpr double kSpacing = 1.5;       ///< > 1.2 × 直径（直径 = 1.0）
constexpr double kCompareTol = 1e-8;   ///< 4 体 vs 独立单体 的最大位置差门槛

/// 诊断可覆盖的运行规格（默认值即预登记规格；诊断运行只在 diagnose/ 下另行保存）
struct RunSpec {
    int steps = kSteps;
    int maxIter = kMaxIter;
    double tol = kTol;
    // --- 停止判据诊断（默认 step = 原判据，行为不变） ---
    std::string stopMode = "step";
    double gradTol = 1e-7;
    double gradScale = std::numeric_limits<double>::quiet_NaN(); ///< NaN -> 用该场景总重
    double actualStepTol = 1e-12;
    int actualStepCount = 2;
};

// -----------------------------------------------------------------------------
// 每个 tet 的 J = det(F)，用于逐步 min/max J 诊断
// -----------------------------------------------------------------------------
class TetJacobian {
public:
    explicit TetJacobian(const MultiBodyScene& scene) : tets_(&scene.tets) {
        const int T = scene.numTets();
        G_.resize(static_cast<size_t>(T));
        for (int t = 0; t < T; ++t) {
            const auto& tet = (*tets_)[t];
            const Vec3 p0 = scene.x0.row(tet[0]).transpose();
            Eigen::Matrix3d Dm;
            Dm.col(0) = scene.x0.row(tet[1]).transpose() - p0;
            Dm.col(1) = scene.x0.row(tet[2]).transpose() - p0;
            Dm.col(2) = scene.x0.row(tet[3]).transpose() - p0;
            if (std::abs(Dm.determinant()) > 1e-18) G_[t] = Dm.inverse();
            else G_[t] = Eigen::Matrix3d::Identity();
        }
    }

    void eval(const Positions& x, double& minJ, double& maxJ) const {
        minJ = std::numeric_limits<double>::infinity();
        maxJ = -minJ;
        for (size_t t = 0; t < tets_->size(); ++t) {
            const auto& tet = (*tets_)[t];
            const Vec3 p0 = x.row(tet[0]).transpose();
            Eigen::Matrix3d Ds;
            Ds.col(0) = x.row(tet[1]).transpose() - p0;
            Ds.col(1) = x.row(tet[2]).transpose() - p0;
            Ds.col(2) = x.row(tet[3]).transpose() - p0;
            const double J = (Ds * G_[t]).determinant();
            minJ = std::min(minJ, J);
            maxJ = std::max(maxJ, J);
        }
    }

private:
    const std::vector<std::array<int, 4>>* tets_ = nullptr;
    std::vector<Eigen::Matrix3d> G_;
};

// -----------------------------------------------------------------------------
// 逐步指标
// -----------------------------------------------------------------------------
struct StepRow {
    int step = 0;
    bool converged = false;
    int iterations = 0;
    int lineSearchIterations = 0;
    double minJ = 0.0, maxJ = 0.0;
    double maxDisp = 0.0, maxSpeed = 0.0, dbcError = 0.0, wallSeconds = 0.0;
    bool maxIterHit = false, lineSearchFailed = false;
    double searchDirectionRate = 0.0;   ///< 停止判据 \(||p||_inf / h\)
    double freeGradientInf = 0.0;       ///< 自由节点梯度无穷范数（不含罚项）
    double forceResidualInf = 0.0, relativeForce = 0.0;
    bool gatePass = false;
    double stateChangeInf = 0.0;        ///< 本步实际最大节点位移（诊断）
    int dbcSatisfied = 0, dbcTotal = 0;
    std::string solver;
    std::string stopReason;             ///< step|grad|actual|combined|maxIter|lineSearch|linearSolve
};

struct SceneRun {
    std::vector<StepRow> rows;
    std::vector<Positions> traj;   ///< traj[0] = 初态，traj[s] = 第 s 步后
    Positions initialX;
    double totalWall = 0.0;
};

struct BodyBox {
    std::array<double, 6> box{};
};

inline BodyBox bodyBox(const Positions& x, const BodyRecord& rec) {
    BodyBox b;
    b.box = {1e300, 1e300, 1e300, -1e300, -1e300, -1e300};
    for (int i = rec.nodeBegin; i < rec.nodeEnd; ++i)
        for (int d = 0; d < kDim; ++d) {
            b.box[d] = std::min(b.box[d], x(i, d));
            b.box[3 + d] = std::max(b.box[3 + d], x(i, d));
        }
    return b;
}

/// 两个轴对齐 bbox 的分离下界：max_axis gap（>0 表示该轴完全分开，
/// 且欧氏距离 >= 该 gap，可用作“无接触”的严格下界证据）
inline double boxGapLowerBound(const BodyBox& a, const BodyBox& b) {
    double best = -std::numeric_limits<double>::infinity();
    for (int d = 0; d < kDim; ++d) {
        const double gap = std::max(b.box[d] - a.box[3 + d], a.box[d] - b.box[3 + d]);
        best = std::max(best, gap);
    }
    return best;
}

/// 表面节点两两最小距离（经验值；不是严格下界，用于对照 bbox 下界）
inline double minSurfaceNodeDistance(const Positions& x,
                                     const std::vector<std::array<int, 3>>& surface,
                                     const BodyRecord& a, const BodyRecord& b) {
    std::vector<int> na, nb;
    for (const auto& f : surface) {
        for (int k = 0; k < 3; ++k) {
            const int g = f[k];
            if (g >= a.nodeBegin && g < a.nodeEnd) na.push_back(g);
            if (g >= b.nodeBegin && g < b.nodeEnd) nb.push_back(g);
        }
    }
    std::sort(na.begin(), na.end());
    na.erase(std::unique(na.begin(), na.end()), na.end());
    std::sort(nb.begin(), nb.end());
    nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    double best = std::numeric_limits<double>::infinity();
    for (int i : na)
        for (int j : nb) {
            const double d = (x.row(i) - x.row(j)).norm();
            if (d < best) best = d;
        }
    return best;
}

// -----------------------------------------------------------------------------
// 运行一个全局场景：steps 步，每步记录指标；返回逐步轨迹
// -----------------------------------------------------------------------------
SceneRun runScene(const MultiBodyScene& scene, const std::string& traceCsvPath,
                  const RunSpec& spec) {
    SceneRun out;
    TimeIntegrator::Config cfg;
    cfg.h = kDt;
    cfg.tol = spec.tol;
    cfg.maxIter = spec.maxIter;
    cfg.verbose = false;
    cfg.traceCsvPath = traceCsvPath;
    cfg.traceLabel = "S6a";
    // 修正停止判据：gradScale 缺省取单体重（1 体 = 518.4769746175648 * 9.81 N）
    cfg.stopMode = spec.stopMode;
    cfg.gradTolRel = spec.gradTol;
    cfg.gradScale = std::isfinite(spec.gradScale) ? spec.gradScale
                                                  : scene.totalMass / scene.bodies.size() * std::abs(kGravityZ);
    cfg.actualStepTol = spec.actualStepTol;
    cfg.actualStepCount = spec.actualStepCount;

    auto elastic = std::make_shared<NeoHookeanTetEnergy>(scene.tets, scene.x0, kYoungs, kPoisson);
    auto gravity = std::make_shared<ExampleGravityEnergy>(scene.masses, Vec3(0.0, 0.0, kGravityZ));
    TimeIntegrator solver(cfg, scene.masses);
    solver.addEnergy(elastic);
    solver.addEnergy(gravity);

    const DBCInfo dbc = scene.makeDBC(kDbcK);
    TetJacobian jac(scene);

    Positions x = scene.x0;
    Velocities v = Velocities::Zero(scene.numNodes(), kDim);
    out.initialX = x;
    out.traj.push_back(x);

    const auto batchBegin = std::chrono::steady_clock::now();
    for (int step = 1; step <= spec.steps; ++step) {
        const auto t0 = std::chrono::steady_clock::now();
        const Positions xPrev = x;
        const auto state = solver.stepForward(x, v, dbc);
        if (!traceCsvPath.empty()) {
            const std::string prefix=traceCsvPath+".step-"+std::to_string(step);
            std::filesystem::copy_file(traceCsvPath,prefix+".csv");
            std::filesystem::copy_file(traceCsvPath+".nodes.csv",prefix+".nodes.csv");
        }
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        x = state.x;
        v = state.v;

        StepRow row;
        row.step = step;
        row.converged = state.report.converged;
        row.iterations = state.report.iterations;
        row.lineSearchIterations = state.report.lineSearchIterations;
        row.dbcError = state.report.dbcMaxError;
        row.maxIterHit = state.report.maxIterHit;
        row.lineSearchFailed = state.report.lineSearchFailed;
        row.solver = state.report.solver;
        row.stopReason = state.report.stopReason;
        row.wallSeconds = wall;
        row.searchDirectionRate = state.report.searchDirectionRate;
        row.freeGradientInf = state.report.freeGradientInf;
        row.forceResidualInf = state.report.forceResidualInf;
        row.relativeForce = row.forceResidualInf / cfg.gradScale;
        row.dbcSatisfied = state.report.dbcSatisfied;
        row.dbcTotal = state.report.dbcTotal;
        jac.eval(x, row.minJ, row.maxJ);
        row.maxDisp = 0.0;
        row.maxSpeed = 0.0;
        row.stateChangeInf = 0.0;
        for (int i = 0; i < scene.numNodes(); ++i) {
            row.maxDisp = std::max(row.maxDisp, (x.row(i) - out.initialX.row(i)).norm());
            row.maxSpeed = std::max(row.maxSpeed, v.row(i).norm());
            row.stateChangeInf = std::max(row.stateChangeInf, (x.row(i) - xPrev.row(i)).norm());
        }
        row.gatePass = row.converged && !row.maxIterHit && !row.lineSearchFailed && x.allFinite() && v.allFinite() && std::isfinite(row.minJ) && std::isfinite(row.maxJ) && row.minJ > .9 && row.dbcError < 1e-6 && row.dbcSatisfied == row.dbcTotal && std::isfinite(row.relativeForce) && row.relativeForce < 1e-6;
        out.rows.push_back(row);
        out.traj.push_back(x);
        if (!row.gatePass) break;
    }
    out.totalWall = std::chrono::duration<double>(std::chrono::steady_clock::now() - batchBegin).count();
    return out;
}

// -----------------------------------------------------------------------------
// 全局 Hessian 块对角检查
//   直接调用能量 addHess 拿 triplet，另按 integrator 的公式补惯性/罚对角项，
//   然后按节点 body_id 检查 (i,j) 是否跨体；同时检查每个 tet 的 4 节点同体。
// -----------------------------------------------------------------------------
struct BlockCheck {
    int numBodies = 0;
    int numNodes = 0;
    int numTets = 0;
    int totalTriplets = 0;
    int crossBodyTriplets = 0;
    double crossBodyMaxAbs = 0.0;
    int assembledNonzeros = 0;
    int assembledCrossEntries = 0;
    double assembledCrossMaxAbs = 0.0;
    int tetsSpanningBodies = 0;
    bool nodeBodyContiguous = true;
    bool pass = false;
};

BlockCheck checkBlockDiagonal(const MultiBodyScene& scene, const Positions& x) {
    BlockCheck c;
    c.numBodies = scene.numBodies();
    c.numNodes = scene.numNodes();
    c.numTets = scene.numTets();

    // 与 TimeIntegrator 相同的装配项：惯性对角 + NH Hessian(h^2) + DBC 罚对角(k*m)
    const std::vector<std::array<int, 4>>& tets = scene.tets;
    const Positions& x0 = scene.x0;
    NeoHookeanTetEnergy elastic(tets, x0, kYoungs, kPoisson);
    std::vector<Triplet> trips;
    trips.reserve(static_cast<size_t>(c.numNodes) * kDim * 4 +
                  static_cast<size_t>(c.numTets) * 144);

    for (int i = 0; i < c.numNodes; ++i)
        for (int d = 0; d < kDim; ++d)
            trips.emplace_back(kDim * i + d, kDim * i + d, scene.masses[i]);

    const size_t before = trips.size();
    elastic.addHess(x, trips);
    const double w = kDt * kDt;   // hessianWeight
    for (size_t i = before; i < trips.size(); ++i) {
        const auto& t = trips[i];
        trips[i] = Triplet(t.row(), t.col(), t.value() * w);
    }
    for (int idx : scene.dbcIndices)
        for (int d = 0; d < kDim; ++d) {
            const int j = kDim * idx + d;
            trips.emplace_back(j, j, kDbcK * scene.masses[idx]);
        }

    c.totalTriplets = static_cast<int>(trips.size());
    for (const auto& t : trips) {
        const int bi = scene.nodeBody[t.row() / kDim];
        const int bj = scene.nodeBody[t.col() / kDim];
        if (bi != bj) {
            ++c.crossBodyTriplets;
            c.crossBodyMaxAbs = std::max(c.crossBodyMaxAbs, std::abs(t.value()));
        }
    }

    const int dof = c.numNodes * kDim;
    SparseMatrix H(dof, dof);
    H.setFromTriplets(trips.begin(), trips.end());
    H.makeCompressed();
    c.assembledNonzeros = static_cast<int>(H.nonZeros());
    for (int col = 0; col < H.outerSize(); ++col)
        for (SparseMatrix::InnerIterator it(H, col); it; ++it) {
            const int bi = scene.nodeBody[it.row() / kDim];
            const int bj = scene.nodeBody[it.col() / kDim];
            if (bi != bj) {
                ++c.assembledCrossEntries;
                c.assembledCrossMaxAbs = std::max(c.assembledCrossMaxAbs, std::abs(it.value()));
            }
        }

    // 能量层面替代检查：每个 tet 的 4 个节点必须同 body
    for (const auto& tet : tets) {
        const int b = scene.nodeBody[tet[0]];
        for (int k = 1; k < 4; ++k)
            if (scene.nodeBody[tet[k]] != b) { ++c.tetsSpanningBodies; break; }
    }

    // 全局节点区间必须按 body 连续（拼接正确性）
    for (int i = 1; i < c.numNodes; ++i)
        if (scene.nodeBody[i] < scene.nodeBody[i - 1]) c.nodeBodyContiguous = false;

    c.pass = (c.crossBodyTriplets == 0) && (c.assembledCrossEntries == 0) &&
             (c.tetsSpanningBodies == 0) && c.nodeBodyContiguous;
    return c;
}

// -----------------------------------------------------------------------------
// 4 次独立单体求解：相同局部初态/载荷/边界
// -----------------------------------------------------------------------------
struct IndependentRun {
    std::vector<Positions> traj;   ///< 局部坐标；traj[0] 初态
    std::vector<StepRow> rows;
};

IndependentRun runIndependent(const MshMesh& mesh, const MassSpringModel& base,
                              const RunSpec& spec) {
    IndependentRun out;
    TimeIntegrator::Config cfg;
    cfg.h = kDt;
    cfg.tol = spec.tol;
    cfg.maxIter = spec.maxIter;
    cfg.verbose = false;
    // 独立单体对照必须与全局场景使用同一停止判据，否则对照不成立。
    // 调用方传入的 gradScale 是单体力尺度，独立复制保持不变。
    cfg.stopMode = spec.stopMode;
    cfg.gradTolRel = spec.gradTol;
    cfg.gradScale = std::isfinite(spec.gradScale) ? spec.gradScale
                                                  : base.masses.sum() * std::abs(kGravityZ);
    cfg.actualStepTol = spec.actualStepTol;
    cfg.actualStepCount = spec.actualStepCount;

    auto elastic = std::make_shared<NeoHookeanTetEnergy>(mesh.tets, mesh.x, kYoungs, kPoisson);
    auto gravity = std::make_shared<ExampleGravityEnergy>(base.masses, Vec3(0.0, 0.0, kGravityZ));
    TimeIntegrator solver(cfg, base.masses);
    solver.addEnergy(elastic);
    solver.addEnergy(gravity);

    DBCInfo dbc;
    dbc.stiffness = kDbcK;
    dbc.isDBC.assign(static_cast<size_t>(mesh.x.rows()), 0);
    for (int i = 0; i < mesh.x.rows(); ++i) {
        if (mesh.x(i, 2) >= kCapZMin) {
            dbc.indices.push_back(i);
            dbc.targetPositions.push_back(mesh.x.row(i).transpose());
            dbc.targetVelocities.push_back(Vec3::Zero());
            dbc.limits.push_back(Vec3::Zero());
            dbc.isDBC[i] = 1;
        }
    }
    if (dbc.indices.empty()) throw std::runtime_error("独立单体：顶部 cap 为空");

    Positions x = mesh.x;
    Velocities v = Velocities::Zero(x.rows(), kDim);
    out.traj.push_back(x);
    for (int step = 1; step <= spec.steps; ++step) {
        const auto state = solver.stepForward(x, v, dbc);
        x = state.x;
        v = state.v;
        StepRow row;
        row.step = step;
        row.converged = state.report.converged;
        row.iterations = state.report.iterations;
        row.lineSearchIterations = state.report.lineSearchIterations;
        row.dbcError = state.report.dbcMaxError;
        row.maxIterHit = state.report.maxIterHit;
        row.lineSearchFailed = state.report.lineSearchFailed;
        row.solver = state.report.solver;
        row.stopReason = state.report.stopReason;
        row.searchDirectionRate = state.report.searchDirectionRate;
        row.freeGradientInf = state.report.freeGradientInf;
        row.dbcSatisfied = state.report.dbcSatisfied;
        row.dbcTotal = state.report.dbcTotal;
        out.rows.push_back(row);
        out.traj.push_back(x);
    }
    return out;
}

// -----------------------------------------------------------------------------
// 输出辅助
// -----------------------------------------------------------------------------
void writeStepsCsv(const std::string& path, const std::vector<StepRow>& rows) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法创建 " + path);
    out << std::setprecision(17);
    out << "step,converged,iterations,line_search_iterations,minJ,maxJ,max_disp,max_speed,"
           "dbc_error,wall_seconds,max_iter_hit,line_search_failed,solver,"
           "search_direction_rate,free_gradient_inf,state_change_inf,dbc_satisfied,dbc_total,"
           "stop_reason,force_residual_inf_N,relative_force_residual,gate_pass\n";
    for (const auto& r : rows)
        out << r.step << ',' << (r.converged ? 1 : 0) << ',' << r.iterations << ','
            << r.lineSearchIterations << ',' << r.minJ << ',' << r.maxJ << ',' << r.maxDisp
            << ',' << r.maxSpeed << ',' << r.dbcError << ',' << r.wallSeconds << ','
            << (r.maxIterHit ? 1 : 0) << ',' << (r.lineSearchFailed ? 1 : 0) << ','
            << r.solver << ',' << r.searchDirectionRate << ',' << r.freeGradientInf << ','
            << r.stateChangeInf << ',' << r.dbcSatisfied << ',' << r.dbcTotal << ','
            << r.stopReason << ',' << r.forceResidualInf << ',' << r.relativeForce << ',' << r.gatePass << '\n';
    out.close();
    if (!out) throw std::runtime_error("写入失败 " + path);
}

void writeStateCsv(const std::string& path, const MultiBodyScene& scene,
                   const Positions& x, const Positions& x0) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法创建 " + path);
    out << std::setprecision(17);
    out << "node,body_id,global_x,global_y,global_z,local_x,local_y,local_z,"
           "disp_x,disp_y,disp_z\n";
    for (int i = 0; i < scene.numNodes(); ++i) {
        const int b = scene.nodeBody[i];
        const Vec3 t = scene.bodies[b].translation;
        const Vec3 g = x.row(i).transpose();
        const Vec3 l = g - t;
        const Vec3 d = g - x0.row(i).transpose();
        out << i << ',' << b << ',' << g.x() << ',' << g.y() << ',' << g.z() << ','
            << l.x() << ',' << l.y() << ',' << l.z() << ','
            << d.x() << ',' << d.y() << ',' << d.z() << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("写入失败 " + path);
}

struct TrajectoryCompare {
    double maxDiff = 0.0;
    std::vector<double> perBodyMax;                 ///< [nb]
    std::vector<std::vector<double>> diff;          ///< [nb][step] 最大节点位置差
};

TrajectoryCompare compareTrajectory(const MultiBodyScene& scene,
                                    const SceneRun& global,
                                    const std::vector<IndependentRun>& indep) {
    TrajectoryCompare cmp;
    const int steps = static_cast<int>(global.traj.size());
    cmp.perBodyMax.assign(scene.numBodies(), 0.0);
    cmp.diff.assign(static_cast<size_t>(scene.numBodies()),
                    std::vector<double>(static_cast<size_t>(steps), 0.0));
    for (int b = 0; b < scene.numBodies(); ++b) {
        const BodyRecord& rec = scene.bodies[b];
        const Vec3 t = rec.translation;
        for (int s = 0; s < steps; ++s) {
            const Positions& xg = global.traj[s];
            const Positions& xl = indep[b].traj[s];
            if (xl.rows() != rec.numNodes())
                throw std::runtime_error("独立解节点数不一致");
            double localMax = 0.0;
            for (int i = 0; i < rec.numNodes(); ++i) {
                const Vec3 g = xg.row(rec.nodeBegin + i).transpose() - t;
                const Vec3 l = xl.row(i).transpose();
                localMax = std::max(localMax, (g - l).norm());
            }
            cmp.diff[b][s] = localMax;
            cmp.perBodyMax[b] = std::max(cmp.perBodyMax[b], localMax);
        }
    }
    for (double d : cmp.perBodyMax) cmp.maxDiff = std::max(cmp.maxDiff, d);
    return cmp;
}

void writeCompareCsv(const std::string& path, const MultiBodyScene& scene,
                     const TrajectoryCompare& cmp) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法创建 " + path);
    out << std::setprecision(17);
    out << "step,body_id,max_position_diff_m\n";
    for (int b = 0; b < scene.numBodies(); ++b)
        for (size_t s = 0; s < cmp.diff[b].size(); ++s)
            out << s << ',' << b << ',' << cmp.diff[b][s] << '\n';
    out.close();
    if (!out) throw std::runtime_error("写入失败 " + path);
}

std::string joinPath(const std::string& dir, const std::string& name) {
    return (std::filesystem::path(dir) / name).string();
}

} // namespace


int main(int argc,char** argv) {
 try {
  if(argc!=3)throw std::runtime_error("usage: s6a_four_main MESH OUT");
  std::string out=argv[2];std::filesystem::create_directories(out);
  const auto mesh=loadMsh(argv[1]);const auto audit=auditMesh(mesh,kDensity);
  if(!audit.positiveOrientation || audit.degenerateTets || !audit.surfaceClosedManifold)throw std::runtime_error("mesh audit failed");
  std::vector<BodyTransform> transforms(4);
  transforms[1].translation=Vec3(1.5,0,0);transforms[2].translation=Vec3(0,1.5,0);transforms[3].translation=Vec3(1.5,1.5,0);
  const auto scene=buildMultiBodyScene(mesh,fnv1aFileHash(argv[1]),transforms,kDensity,kCapZMin);
  if(scene.numNodes()!=7040 || scene.numTets()!=27404 || scene.makeDBC(kDbcK).indices.size()!=252)throw std::runtime_error("counts mismatch");
  RunSpec spec;spec.stopMode="combined";spec.gradTol=1e-6;spec.gradScale=scene.totalMass/4*9.81;
  auto writeBlock=[&](const BlockCheck& b,const std::string& name){std::ofstream f(joinPath(out,name));f<<"{\"cross_triplets\":"<<b.crossBodyTriplets<<",\"cross_entries\":"<<b.assembledCrossEntries<<",\"tets_spanning\":"<<b.tetsSpanningBodies<<",\"nonzeros\":"<<b.assembledNonzeros<<",\"pass\":"<<(b.pass?"true":"false")<<"}\n";f.close();if(!f)throw std::runtime_error("block write failed");};
  auto initial=checkBlockDiagonal(scene,scene.x0);writeBlock(initial,"block_initial.json");if(!initial.pass)return 1;
  auto result=runScene(scene,joinPath(out,"four_iterations.csv"),spec);
  writeStepsCsv(joinPath(out,"steps.csv"),result.rows);writeStateCsv(joinPath(out,"state.csv"),scene,result.traj.back(),scene.x0);
  if(result.rows.size()!=10)return 1;for(auto& r:result.rows)if(!r.gatePass)return 1;
  auto final=checkBlockDiagonal(scene,result.traj.back());writeBlock(final,"block_final.json");if(!final.pass)return 1;
  const auto one=buildMultiBodyScene(mesh,fnv1aFileHash(argv[1]),std::vector<BodyTransform>{BodyTransform{}},kDensity,kCapZMin);
  std::vector<IndependentRun> independent;
  for(int b=0;b<4;++b){
   std::string prefix="independent_"+std::to_string(b);
   auto r=runScene(one,joinPath(out,prefix+"_iterations.csv"),spec);
   writeStepsCsv(joinPath(out,prefix+"_steps.csv"),r.rows);writeStateCsv(joinPath(out,prefix+"_state.csv"),one,r.traj.back(),one.x0);
   if(r.rows.size()!=10)return 1;for(auto& row:r.rows)if(!row.gatePass)return 1;
   independent.push_back({r.traj,r.rows});
  }
  auto compare=compareTrajectory(scene,result,independent);writeCompareCsv(joinPath(out,"trajectory_compare.csv"),scene,compare);
  std::cout<<std::setprecision(17)<<"four-body max independent difference="<<compare.maxDiff<<std::endl;
  return compare.maxDiff<1e-8?0:1;
 }catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 2;}
}
