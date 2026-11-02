// =============================================================================
// dynamic_bar_check.cpp — S3 动态大挠度弯曲 / S3t 扭转 DBC 驱动（单杆）
// =============================================================================
// 预登记规格：docs/s3-s6a-execution-design-20261031.md 第 1、2 节。
//
//   S3 (--mode bend)：
//     M2 结构化方杆（20×4×4），NH E=4.4e7, nu=0.3, rho=1000；
//     左端 x=0 全端面固定；重力 (0,0,gz)（默认 -9.81）；参考形状零速度起步。
//   S3t (--mode twist)：
//     同几何/材料，无重力；x=1 端面节点绕杆轴 x（过 y=z=0）按
//     theta(t) = theta_max * smoothstep(t/T) 旋转，t>=T 后保持；每步刷新 DBC 绝对目标。
//
// 输出：
//   <out-dir>/steps.csv   每步指标（逐行 flush，崩溃保留前缀）
//   <out-dir>/state.csv   末态逐节点位置/位移
//   <out-dir>/summary.json 门槛判定 + 极值 + 资源
//
// 反力定义（与 cantilever_check 一致）：reaction = 弹性能+重力对夹持节点梯度之和，
//   即约束必须施加给杆的力。反力矩绕 x 轴：Mx = sum_{clamp} [(p_i-c) x g_i]_x，c 为杆轴点 (0,0,0)。
//   dbc_penalty_torque_x：DBC 罚力绕 x 轴的矩。IP 中罚项未乘 h^2，故物理约束力
//     F_i = -k m_i (p_i - t_i)/h^2；drive_torque_x 为驱动端材料梯度力矩（独立复算，静止时两者反号）。
// 一致性检查：reaction_torque_x + drive_torque_x ≈ 0（无重力、角动量可忽略时）。
// 能量：elastic/gravity 为未加权的物理能量（不含 h^2 IP 权重）；total = elastic+gravity（势能合计）。
// =============================================================================
#include "energies/Gravity.h"
#include "energies/NeoHookeanTetEnergy.h"
#include "io/msh_loader.h"
#include "scene/mass_spring_builder.h"
#include "time_integrator.h"

#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace fiction;

constexpr double kE = 4.40e7, kNu = 0.3, kRho = 1000.0;
constexpr double kL = 1.0, kHalf = 0.025, kVolume = 0.0025, kMass = 2.5;
constexpr double kPi = 3.14159265358979323846;

struct Options {
    std::string meshPath;
    std::string outDir;
    double dt = 0.005;
    int steps = 2000;
    double tol = 1e-4;
    int maxIter = 50;
    std::string mode = "bend";
    double twistDeg = 90.0;
    double twistTime = 2.0;
    double gravityZ = -9.81;
    double initialK = 1e8;
    bool verbose = false;
};

const char* kUsage =
    "dynamic_bar_check --mesh FILE --out-dir DIR [--dt 0.005] [--steps 2000]\n"
    "  [--tol 1e-4] [--max-iter 50] [--mode bend|twist] [--twist-deg 90]\n"
    "  [--twist-time 2.0] [--gravity-z -9.81] [--verbose]\n";

double smoothstep(double s) {
    s = std::min(1.0, std::max(0.0, s));
    return s * s * (3.0 - 2.0 * s);
}

std::string jsonEscape(const std::string& s) {
    std::string out;
    for (char c : s) {
        if (c == '"' || c == '\\') out += '\\';
        out += c;
    }
    return out;
}

long peakRssKb() {
    std::ifstream in("/proc/self/status");
    std::string key;
    while (in >> key) {
        if (key == "VmHWM:") {
            long kb = 0;
            in >> kb;
            return kb;
        }
        in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    return -1;
}

struct StepMetric {
    double tipDeflection = 0.0;
    double tipMeanSpeed = 0.0;
    double tipMaxSpeed = 0.0;
    double maxSpeed = 0.0;
    int converged = 0, iterations = 0, lineSearchIterations = 0;
    int maxIterHit = 0, lineSearchFailed = 0;
    double dbcError = 0.0;          ///< 独立复算的 DBC 目标位置误差（max L2）
    double reportDbcError = 0.0;    ///< 求解器自报值，用于交叉核对
    double minJ = 0.0, maxJ = 0.0;
    double minStretch = 0.0, maxStretch = 0.0;
    double elastic = 0.0, gravity = 0.0, total = 0.0, kinetic = 0.0;
    double reactionX = 0.0, reactionY = 0.0, reactionZ = 0.0;
    double reactionTorqueX = 0.0, penaltyTorqueX = 0.0, driveTorqueX = 0.0;
    double targetDeg = 0.0, actualDeg = 0.0;
    double searchDirectionRate = 0.0, freeGradientInf = 0.0, dbcStiffness = 0.0;
    std::string solver;
};

struct Bounds {
    double minJ = std::numeric_limits<double>::infinity();
    double dbcError = 0.0;
    double absActualDeg = 0.0;
    double minActualDeg = std::numeric_limits<double>::infinity();
    double maxActualDeg = -std::numeric_limits<double>::infinity();
};

} // namespace

int main(int argc, char** argv) {
    Options opt;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") { std::cout << kUsage; return 0; }
            if (arg == "--verbose") { opt.verbose = true; continue; }
            if (i + 1 >= argc) throw std::runtime_error("缺少参数值: " + arg);
            const std::string val = argv[++i];
            if (arg == "--mesh") opt.meshPath = val;
            else if (arg == "--out-dir") opt.outDir = val;
            else if (arg == "--dt") opt.dt = std::stod(val);
            else if (arg == "--steps") opt.steps = std::stoi(val);
            else if (arg == "--tol") opt.tol = std::stod(val);
            else if (arg == "--max-iter") opt.maxIter = std::stoi(val);
            else if (arg == "--mode") opt.mode = val;
            else if (arg == "--twist-deg") opt.twistDeg = std::stod(val);
            else if (arg == "--twist-time") opt.twistTime = std::stod(val);
            else if (arg == "--gravity-z") opt.gravityZ = std::stod(val);
            else throw std::runtime_error("未知参数: " + arg);
        }
        if (opt.meshPath.empty() || opt.outDir.empty())
            throw std::runtime_error("必须指定 --mesh 和 --out-dir");
        if (opt.mode != "bend" && opt.mode != "twist")
            throw std::runtime_error("--mode 仅支持 bend | twist");
        if (!(opt.dt > 0.0) || opt.steps <= 0 || !(opt.tol > 0.0) || opt.maxIter <= 0)
            throw std::runtime_error("--dt/--steps/--tol/--max-iter 非法");
        if (opt.mode == "twist" && (!(opt.twistTime > 0.0) || opt.twistDeg <= 0.0))
            throw std::runtime_error("--twist-deg/--twist-time 必须为正");
    } catch (const std::exception& e) {
        std::cerr << "S3 ARG ERROR: " << e.what() << '\n' << kUsage;
        return 2;
    }

    const bool twisting = (opt.mode == "twist");
    const bool useGravity = !twisting;                 // S3t 规格：无重力
    const Vec3 gravity(0.0, 0.0, useGravity ? opt.gravityZ : 0.0);
    const double thetaMax = opt.twistDeg * kPi / 180.0;
    const Vec3 axisPoint(0.0, 0.0, 0.0);               // 杆轴 x 过 y=z=0

    try {
        std::filesystem::create_directories(opt.outDir);
        const MshMesh mesh = loadMsh(opt.meshPath);
        if (mesh.empty() || mesh.tets.empty() || !mesh.x.allFinite())
            throw std::runtime_error("空网格或非有限坐标");
        const MassSpringModel model = buildMassSpringModel(mesh, kRho);

        // --- 运行前几何/质量核对（预登记规格）---
        const std::array<double, 6> box{{0.0, -kHalf, -kHalf, kL, kHalf, kHalf}};
        for (int k = 0; k < 6; ++k)
            if (std::abs(model.bbox[k] - box[k]) > 1e-12)
                throw std::runtime_error("包围盒与 S3 规格不符");
        if (std::abs(model.volume - kVolume) > 1e-12 * kVolume)
            throw std::runtime_error("体积与 0.0025 m^3 不符");
        if (std::abs(model.totalMass - kMass) > 1e-12 * kMass)
            throw std::runtime_error("总质量与 2.5 kg 不符");
        auto edges = [](const Positions& x, const std::array<int, 4>& t) {
            Eigen::Matrix3d d;
            for (int c = 0; c < 3; ++c) d.col(c) = (x.row(t[c + 1]) - x.row(t[0])).transpose();
            return d;
        };
        for (const auto& t : mesh.tets)
            if (!(edges(mesh.x, t).determinant() > 0.0))
                throw std::runtime_error("参考四面体非正定向");

        // --- 端面节点 ---
        std::vector<int> tip, clamp;
        for (int i = 0; i < mesh.x.rows(); ++i) {
            if (std::abs(mesh.x(i, 0)) < 1e-12) clamp.push_back(i);
            if (std::abs(mesh.x(i, 0) - kL) < 1e-12) tip.push_back(i);
        }
        if (clamp.empty() || tip.empty()) throw std::runtime_error("缺少端面节点");

        DBCInfo dbc;
        dbc.stiffness = opt.initialK;
        dbc.isDBC.assign(mesh.x.rows(), 0);
        for (int i : clamp) {
            dbc.indices.push_back(i);
            dbc.targetVelocities.push_back(Vec3::Zero());
            dbc.limits.push_back(mesh.x.row(i).transpose());
            dbc.isDBC[i] = 1;
        }
        const size_t clampCount = dbc.indices.size();
        if (twisting) {
            for (int i : tip) {
                dbc.indices.push_back(i);
                dbc.targetVelocities.push_back(Vec3::Zero());
                dbc.limits.push_back(mesh.x.row(i).transpose());
                dbc.isDBC[i] = 1;
            }
        }
        dbc.targetPositions.reserve(dbc.indices.size());
        for (int idx : dbc.indices) dbc.targetPositions.push_back(mesh.x.row(idx).transpose());

        // --- 能量：NH + 重力 ---
        auto elastic = std::make_shared<NeoHookeanTetEnergy>(mesh.tets, mesh.x, kE, kNu);
        std::shared_ptr<ExampleGravityEnergy> gravityEnergy;
        if (useGravity) gravityEnergy = std::make_shared<ExampleGravityEnergy>(model.masses, gravity);

        TimeIntegrator::Config cfg;
        cfg.h = opt.dt;
        cfg.tol = opt.tol;
        cfg.maxIter = opt.maxIter;
        cfg.verbose = opt.verbose;
        cfg.staticSolve = false;
        TimeIntegrator solver(cfg, model.masses);
        solver.addEnergy(elastic);
        if (gravityEnergy) solver.addEnergy(gravityEnergy);

        // --- 参考四面体 Dm^{-1}（用于 F 与奇异值）---
        std::vector<Eigen::Matrix3d> invDm(mesh.tets.size());
        for (size_t t = 0; t < mesh.tets.size(); ++t) {
            const auto& tet = mesh.tets[t];
            const Vec3 p0 = mesh.x.row(tet[0]).transpose();
            Eigen::Matrix3d Dm;
            Dm.col(0) = mesh.x.row(tet[1]).transpose() - p0;
            Dm.col(1) = mesh.x.row(tet[2]).transpose() - p0;
            Dm.col(2) = mesh.x.row(tet[3]).transpose() - p0;
            invDm[t] = Dm.inverse();
        }

        const std::string stepsPath = (std::filesystem::path(opt.outDir) / "steps.csv").string();
        std::ofstream csv(stepsPath);
        if (!csv) throw std::runtime_error("无法创建 " + stepsPath);
        csv << std::setprecision(17);
        csv << "# protocol=fiction-validation-v0.3-S3; mode=" << opt.mode
            << "; E=" << kE << ",nu=" << kNu << ",rho=" << kRho
            << ",gz=" << (useGravity ? opt.gravityZ : 0.0)
            << "; dt=" << opt.dt << ",tol=" << opt.tol << ",max_iter=" << opt.maxIter
            << ",initial_dbc_stiffness=" << opt.initialK
            << "; tip_deflection=mean_downward_z_of_x1_face_nodes\n"
            << "step,time,tip_deflection,tip_mean_speed,max_speed,converged,iterations,"
               "line_search_iterations,max_iter_hit,line_search_failed,dbc_error,min_J,max_J,"
               "min_stretch,max_stretch,elastic_energy,gravity_energy,total_energy,"
               "reaction_x,reaction_y,reaction_z,wall_seconds,"
               "tip_max_speed,report_dbc_error,search_direction_rate,free_gradient_inf,"
               "dbc_stiffness,kinetic_energy,reaction_torque_x,dbc_penalty_torque_x,"
               "drive_torque_x,target_rotation_deg,actual_rotation_deg,solver\n";

        Positions x = mesh.x;
        Velocities v = Velocities::Zero(mesh.x.rows(), 3);
        StepMetric last;
        Bounds worst;
        bool allConverged = true, allFinite = true, gatesOk = true;
        int failedStep = -1;
        double maxJump = 0.0, prevTip = 0.0;
        double totalWall = 0.0, solveWall = 0.0;
        std::vector<Vec3> prevDrivePos;
        std::string solverPath = "-";
        // 作用在杆上的外载荷力（含符号）；静止时 reaction + external 应为 0
        const Vec3 external = useGravity ? Vec3(0.0, 0.0, model.totalMass * opt.gravityZ) : Vec3::Zero();

        for (int step = 1; step <= opt.steps; ++step) {
            const double t = step * opt.dt;
            double targetDeg = 0.0;
            if (twisting) {
                targetDeg = opt.twistDeg * smoothstep(t / opt.twistTime);
                const double theta = thetaMax * smoothstep(t / opt.twistTime);
                const Eigen::Matrix3d R = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitX()).toRotationMatrix();
                for (size_t k = clampCount; k < dbc.indices.size(); ++k) {
                    const int idx = dbc.indices[k];
                    const Vec3 p = mesh.x.row(idx).transpose();
                    dbc.targetPositions[k] = axisPoint + R * (p - axisPoint);
                }
            }

            const auto begin = std::chrono::steady_clock::now();
            const TimeIntegrator::State state = solver.stepForward(x, v, dbc);
            const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count();
            totalWall += wall;
            solveWall += wall;
            x = state.x;
            v = state.v;
            const auto& r = state.report;
            solverPath = r.solver;

            StepMetric m;
            m.converged = r.converged ? 1 : 0;
            m.iterations = r.iterations;
            m.lineSearchIterations = r.lineSearchIterations;
            m.maxIterHit = r.maxIterHit ? 1 : 0;
            m.lineSearchFailed = r.lineSearchFailed ? 1 : 0;
            m.reportDbcError = r.dbcMaxError;
            m.searchDirectionRate = r.searchDirectionRate;
            m.freeGradientInf = r.freeGradientInf;
            m.dbcStiffness = r.dbcStiffness;
            m.solver = r.solver;

            // --- 尖端面指标 ---
            double sumDown = 0.0, sumSpeed = 0.0, tipMaxSpeed = 0.0;
            for (int i : tip) {
                sumDown += (mesh.x(i, 2) - x(i, 2));
                const double sp = v.row(i).norm();
                sumSpeed += sp;
                tipMaxSpeed = std::max(tipMaxSpeed, sp);
            }
            m.tipDeflection = sumDown / static_cast<double>(tip.size());
            m.tipMeanSpeed = sumSpeed / static_cast<double>(tip.size());
            m.tipMaxSpeed = tipMaxSpeed;
            m.maxSpeed = 0.0;
            m.kinetic = 0.0;
            for (int i = 0; i < x.rows(); ++i) {
                m.maxSpeed = std::max(m.maxSpeed, v.row(i).norm());
                m.kinetic += 0.5 * model.masses[i] * v.row(i).squaredNorm();
            }

            // --- DBC 目标误差（独立复算）---
            m.dbcError = 0.0;
            for (size_t k = 0; k < dbc.indices.size(); ++k) {
                const int idx = dbc.indices[k];
                const double e = (x.row(idx).transpose() - dbc.targetPositions[k]).norm();
                m.dbcError = std::max(m.dbcError, std::isfinite(e) ? e : std::numeric_limits<double>::max());
            }

            // --- F / J / 主伸长 ---
            m.minJ = std::numeric_limits<double>::infinity();
            m.maxJ = -m.minJ;
            m.minStretch = m.minJ;
            m.maxStretch = -m.minJ;
            bool finite = x.allFinite() && v.allFinite();
            for (size_t ti = 0; ti < mesh.tets.size(); ++ti) {
                const Eigen::Matrix3d F = edges(x, mesh.tets[ti]) * invDm[ti];
                const double J = F.determinant();
                finite = finite && F.allFinite() && std::isfinite(J);
                m.minJ = std::min(m.minJ, J);
                m.maxJ = std::max(m.maxJ, J);
                if (F.allFinite()) {
                    const Eigen::Vector3d s = Eigen::JacobiSVD<Eigen::Matrix3d>(F).singularValues();
                    m.minStretch = std::min(m.minStretch, s.minCoeff());
                    m.maxStretch = std::max(m.maxStretch, s.maxCoeff());
                }
            }

            // --- 能量与反力（弹性能+重力对夹持节点梯度之和）---
            m.elastic = elastic->val(x);
            m.gravity = gravityEnergy ? gravityEnergy->val(x) : 0.0;
            m.total = m.elastic + m.gravity;
            Positions grad = Positions::Zero(x.rows(), 3);
            elastic->addGrad(x, grad);
            if (gravityEnergy) gravityEnergy->addGrad(x, grad);
            Vec3 reaction = Vec3::Zero();
            for (size_t k = 0; k < clampCount; ++k) reaction += grad.row(dbc.indices[k]).transpose();
            m.reactionX = reaction.x();
            m.reactionY = reaction.y();
            m.reactionZ = reaction.z();
            // 绕 x 轴的反力矩（夹持端，弹性能+重力梯度）
            double torqueX = 0.0;
            for (size_t k = 0; k < clampCount; ++k) {
                const int idx = dbc.indices[k];
                const Vec3 p = x.row(idx).transpose() - axisPoint;
                const Vec3 g = grad.row(idx).transpose();
                const Vec3 mm = p.cross(g);
                torqueX += mm.x();
            }
            m.reactionTorqueX = torqueX;
            // DBC 罚力矩（全部 DBC 节点）：IP 中罚项未乘 h^2，故物理约束力 = -k m (p-target)/h^2
            double penTorque = 0.0;
            for (size_t k = 0; k < dbc.indices.size(); ++k) {
                const int idx = dbc.indices[k];
                const Vec3 p = x.row(idx).transpose() - axisPoint;
                const Vec3 F = -dbc.stiffness * model.masses[idx] *
                               (x.row(idx).transpose() - dbc.targetPositions[k]) / (opt.dt * opt.dt);
                penTorque += p.cross(F).x();
            }
            m.penaltyTorqueX = penTorque;
            // 驱动端材料梯度力矩（与罚力矩独立计算，静止时应互为反号）
            double driveTorque = 0.0;
            for (size_t k = clampCount; k < dbc.indices.size(); ++k) {
                const int idx = dbc.indices[k];
                const Vec3 p = x.row(idx).transpose() - axisPoint;
                driveTorque += p.cross(grad.row(idx).transpose()).x();
            }
            m.driveTorqueX = driveTorque;

            // --- 扭转角（尖端面 yz 平面相对参考的平均旋转）---
            m.targetDeg = targetDeg;
            if (twisting) {
                double sumAng = 0.0;
                int counted = 0;
                for (int i : tip) {
                    const double y0 = mesh.x(i, 1), z0 = mesh.x(i, 2);
                    const double r0 = std::hypot(y0, z0);
                    if (r0 < 1e-12) continue;
                    const double y1 = x(i, 1) - axisPoint.y(), z1 = x(i, 2) - axisPoint.z();
                    // 绕 +x 的旋转：theta = atan2(r0 × r, r0 · r)，cross 取 yz 平面
                    const double cross2 = y0 * z1 - z0 * y1;
                    const double dot2 = y0 * y1 + z0 * z1;
                    sumAng += std::atan2(cross2, dot2);
                    ++counted;
                }
                m.actualDeg = counted > 0 ? (sumAng / counted) * 180.0 / kPi : 0.0;
            } else {
                m.actualDeg = 0.0;
            }
            last = m;

            // --- 门槛累积 ---
            allConverged = allConverged && (m.converged == 1);
            allFinite = allFinite && finite && std::isfinite(m.minJ) && std::isfinite(m.maxJ) &&
                        std::isfinite(m.minStretch) && std::isfinite(m.maxStretch) &&
                        std::isfinite(m.tipDeflection) && std::isfinite(m.maxSpeed) &&
                        std::isfinite(m.dbcError) && std::isfinite(m.total) && std::isfinite(torqueX);
            worst.minJ = std::min(worst.minJ, m.minJ);
            worst.dbcError = std::max(worst.dbcError, m.dbcError);
            if (twisting) {
                worst.absActualDeg = std::max(worst.absActualDeg, std::abs(m.actualDeg));
                worst.minActualDeg = std::min(worst.minActualDeg, m.actualDeg);
                worst.maxActualDeg = std::max(worst.maxActualDeg, m.actualDeg);
            }
            double jump = 0.0;
            std::vector<Vec3> drivePos;
            if (twisting) {
                drivePos.reserve(dbc.indices.size() - clampCount);
                for (size_t k = clampCount; k < dbc.indices.size(); ++k)
                    drivePos.push_back(x.row(dbc.indices[k]).transpose());
                if (prevDrivePos.size() == drivePos.size()) {
                    for (size_t k = 0; k < drivePos.size(); ++k)
                        jump = std::max(jump, (drivePos[k] - prevDrivePos[k]).norm());
                }
                prevDrivePos = drivePos;
            } else if (step > 10) {
                jump = std::abs(m.tipDeflection - prevTip);
            }
            if (step > 10) maxJump = std::max(maxJump, jump);
            prevTip = m.tipDeflection;

            csv << step << ',' << t << ',' << m.tipDeflection << ',' << m.tipMeanSpeed << ','
                << m.maxSpeed << ',' << m.converged << ',' << m.iterations << ','
                << m.lineSearchIterations << ',' << m.maxIterHit << ',' << m.lineSearchFailed << ','
                << m.dbcError << ',' << m.minJ << ',' << m.maxJ << ',' << m.minStretch << ','
                << m.maxStretch << ',' << m.elastic << ',' << m.gravity << ',' << m.total << ','
                << m.reactionX << ',' << m.reactionY << ',' << m.reactionZ << ',' << wall << ','
                << m.tipMaxSpeed << ',' << m.reportDbcError << ',' << m.searchDirectionRate << ','
                << m.freeGradientInf << ',' << m.dbcStiffness << ',' << m.kinetic << ','
                << m.reactionTorqueX << ',' << m.penaltyTorqueX << ',' << m.driveTorqueX << ',' << m.targetDeg << ','
                << m.actualDeg << ',' << m.solver << '\n';
            csv.flush();
            if (opt.verbose)
                std::cout << "step " << step << '/' << opt.steps << " t=" << t
                          << " tip=" << m.tipDeflection << " it=" << m.iterations
                          << " conv=" << m.converged << " minJ=" << m.minJ
                          << " dbcErr=" << m.dbcError << " wall=" << wall << "s\n";

            if (!(m.converged == 1) || !finite) {
                failedStep = step;
                gatesOk = false;
                std::cerr << "S3 GATE FAIL: 步 " << step << " converged=" << m.converged
                          << " finite=" << finite << " maxIterHit=" << m.maxIterHit
                          << " lineSearchFailed=" << m.lineSearchFailed << '\n';
                break;
            }
            if (!(m.minJ > 0.3)) {
                failedStep = step;
                gatesOk = false;
                std::cerr << "S3 GATE FAIL: 步 " << step << " minJ=" << m.minJ << " <= 0.3\n";
                break;
            }
            if (!(m.dbcError < 1e-6)) {
                failedStep = step;
                gatesOk = false;
                std::cerr << "S3 GATE FAIL: 步 " << step << " dbc_error=" << m.dbcError << " >= 1e-6\n";
                break;
            }
            if (step > 10 && !(jump < 0.02)) {
                failedStep = step;
                gatesOk = false;
                std::cerr << "S3 GATE FAIL: 步 " << step << " 位移跳变=" << jump << " >= 0.02\n";
                break;
            }
        }

        csv.close();
        if (!csv) throw std::runtime_error("steps.csv 写入失败");

        // --- 末态 ---
        const std::string statePath = (std::filesystem::path(opt.outDir) / "state.csv").string();
        std::ofstream stateOut(statePath);
        if (!stateOut) throw std::runtime_error("无法创建 " + statePath);
        stateOut << std::setprecision(17);
        stateOut << "node,X,Y,Z,x,y,z,ux,uy,uz,is_dbc,is_tip\n";
        for (int i = 0; i < x.rows(); ++i) {
            stateOut << i;
            for (int d = 0; d < 3; ++d) stateOut << ',' << mesh.x(i, d);
            for (int d = 0; d < 3; ++d) stateOut << ',' << x(i, d);
            for (int d = 0; d < 3; ++d) stateOut << ',' << x(i, d) - mesh.x(i, d);
            const bool isTip = std::find(tip.begin(), tip.end(), i) != tip.end();
            stateOut << ',' << int(dbc.isDBC[i]) << ',' << (isTip ? 1 : 0) << '\n';
        }
        stateOut.close();
        if (!stateOut) throw std::runtime_error("state.csv 写入失败");

        const int stepsRun = (failedStep > 0 ? failedStep : opt.steps);
        const bool jumpOk = maxJump < 0.02;
        const bool pass = gatesOk && allConverged && allFinite && (worst.minJ > 0.3) &&
                          (worst.dbcError < 1e-6) && jumpOk && failedStep < 0;
        const bool settled = last.maxSpeed < 0.05;

        std::ostringstream js;
        js << std::setprecision(17) << "{\n"
           << "  \"protocol\": \"fiction-validation-v0.3-S3\",\n"
           << "  \"spec\": \"docs/s3-s6a-execution-design-20261031.md\",\n"
           << "  \"run_id\": \"" << jsonEscape(std::filesystem::path(opt.outDir).filename().string()) << "\",\n"
           << "  \"mode\": \"" << opt.mode << "\",\n"
           << "  \"mesh\": \"" << jsonEscape(opt.meshPath) << "\",\n"
           << "  \"params\": {\n"
           << "    \"E\": " << kE << ", \"nu\": " << kNu << ", \"rho\": " << kRho << ",\n"
           << "    \"dt\": " << opt.dt << ", \"steps\": " << opt.steps << ",\n"
           << "    \"tol\": " << opt.tol << ", \"max_iter\": " << opt.maxIter << ",\n"
           << "    \"initial_dbc_stiffness\": " << opt.initialK << ",\n"
           << "    \"gravity\": [" << gravity.x() << ", " << gravity.y() << ", " << gravity.z() << "],\n"
           << "    \"twist_deg\": " << (twisting ? opt.twistDeg : 0.0) << ",\n"
           << "    \"twist_time\": " << (twisting ? opt.twistTime : 0.0) << ",\n"
           << "    \"static_solve\": false\n  },\n"
           << "  \"mesh_check\": {\n"
           << "    \"nodes\": " << x.rows() << ", \"tets\": " << mesh.tets.size() << ",\n"
           << "    \"bbox\": [" << model.bbox[0] << ", " << model.bbox[1] << ", " << model.bbox[2]
           << ", " << model.bbox[3] << ", " << model.bbox[4] << ", " << model.bbox[5] << "],\n"
           << "    \"volume\": " << model.volume << ", \"total_mass\": " << model.totalMass << ",\n"
           << "    \"clamp_nodes\": " << clamp.size() << ", \"tip_nodes\": " << tip.size() << ",\n"
           << "    \"dbc_nodes\": " << dbc.indices.size() << ",\n"
           << "    \"all_tets_positive\": true\n  },\n"
           << "  \"steps_run\": " << stepsRun << ", \"sim_time\": " << stepsRun * opt.dt << ",\n"
           << "  \"failed_step\": " << failedStep << ",\n"
           << "  \"gates\": {\n"
           << "    \"all_steps_converged\": " << (allConverged ? "true" : "false") << ",\n"
           << "    \"finite\": " << (allFinite ? "true" : "false") << ",\n"
           << "    \"min_J\": " << worst.minJ << ", \"min_J_above_0.3\": "
           << (worst.minJ > 0.3 ? "true" : "false") << ",\n"
           << "    \"dbc_error_max\": " << worst.dbcError << ", \"dbc_error_below_1e-6\": "
           << (worst.dbcError < 1e-6 ? "true" : "false") << ",\n"
           << "    \"max_step_jump_after_10\": " << maxJump << ", \"max_step_jump_below_0.02\": "
           << (jumpOk ? "true" : "false") << ",\n"
           << "    \"max_step_jump_definition\": \""
           << (twisting ? "driven DBC node target displacement norm, consecutive steps >10"
                        : "|tip_deflection[k]-tip_deflection[k-1]|, k>10") << "\",\n"
           << "    \"end_max_speed\": " << last.maxSpeed << ",\n"
           << "    \"terminal_speed_below_0.05\": " << (settled ? "true" : "false") << ",\n"
           << "    \"pass\": " << (pass ? "true" : "false") << ",\n"
           << "    \"pass_definition\": \"S3 预登记安全门槛（收敛/有限/minJ>0.3/dbc<1e-6/jump<0.02）；"
              "终止速度 0.05 m/s 只用于判断是否可声称平衡\"\n  },\n"
           << "  \"final\": {\n"
           << "    \"tip_deflection\": " << last.tipDeflection << ",\n"
           << "    \"tip_mean_speed\": " << last.tipMeanSpeed << ",\n"
           << "    \"max_speed\": " << last.maxSpeed << ",\n"
           << "    \"min_J\": " << last.minJ << ", \"max_J\": " << last.maxJ << ",\n"
           << "    \"min_stretch\": " << last.minStretch << ", \"max_stretch\": " << last.maxStretch << ",\n"
           << "    \"elastic_energy\": " << last.elastic << ", \"gravity_energy\": " << last.gravity << ",\n"
           << "    \"total_energy\": " << last.total << ", \"kinetic_energy\": " << last.kinetic << ",\n"
           << "    \"reaction\": [" << last.reactionX << ", " << last.reactionY << ", " << last.reactionZ << "],\n"
           << "    \"reaction_torque_x\": " << last.reactionTorqueX << ",\n"
           << "    \"dbc_penalty_torque_x\": " << last.penaltyTorqueX << ",\n"
           << "    \"drive_torque_x\": " << last.driveTorqueX << ",\n"
           << "    \"target_rotation_deg\": " << last.targetDeg << ",\n"
           << "    \"actual_rotation_deg\": " << last.actualDeg << ",\n"
           << "    \"dbc_error\": " << last.dbcError << "\n  },\n"
           << "  \"rotation_range_deg\": [" << (twisting ? worst.minActualDeg : 0.0) << ", "
           << (twisting ? worst.maxActualDeg : 0.0) << "],\n"
           << "  \"external_force\": [" << external.x() << ", " << external.y() << ", " << external.z() << "],\n"
           << "  \"reaction_plus_external_norm\": "
           << std::sqrt((last.reactionX + external.x()) * (last.reactionX + external.x()) +
                        (last.reactionY + external.y()) * (last.reactionY + external.y()) +
                        (last.reactionZ + external.z()) * (last.reactionZ + external.z())) << ",\n"
           << "  \"solver_paths\": [\"" << jsonEscape(solverPath) << "\"],\n"
           << "  \"wall_seconds_total\": " << totalWall << ",\n"
           << "  \"peak_rss_kb\": " << peakRssKb() << "\n"
           << "}\n";

        const std::string summaryPath = (std::filesystem::path(opt.outDir) / "summary.json").string();
        std::ofstream summary(summaryPath);
        if (!summary) throw std::runtime_error("无法创建 " + summaryPath);
        summary << js.str();
        summary.close();
        if (!summary) throw std::runtime_error("summary.json 写入失败");

        std::cout << std::setprecision(12) << "S3[" << opt.mode << "] steps=" << stepsRun
                  << "/" << opt.steps << " tip=" << last.tipDeflection
                  << " minJ=" << worst.minJ << " maxSpeed=" << last.maxSpeed
                  << " dbcErr=" << worst.dbcError << " maxJump=" << maxJump
                  << " rot=" << last.actualDeg << "deg"
                  << " wall=" << totalWall << "s"
                  << (pass ? " PASS\n" : " FAIL\n");
        (void)solveWall;
        return pass ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "S3 ERROR: " << e.what() << '\n';
        return 3;
    }
}
