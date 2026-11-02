// S1：均匀拉伸 patch test。加载路径与误差阈值见 validation protocol v0.2。
#include "energies/NeoHookeanTetEnergy.h"
#include "io/msh_loader.h"
#include "scene/mass_spring_builder.h"
#include "time_integrator.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace fiction;
constexpr double kYoungs = 1.0e6;
constexpr double kPoisson = 0.3;
constexpr double kLength = 1.0;
constexpr double kWidth = 0.05;
constexpr double kVolume = kLength * kWidth * kWidth;
constexpr double kArea = kWidth * kWidth;
constexpr double kDensity = 1000.0;
constexpr double kTolerance = 1.0e-10;
constexpr double kDbcStiffness = 1.0e8;
constexpr double kEndpointTol = 1.0e-6 * kLength;
constexpr double kMiddleTol = 1.0e-7 * kLength;
constexpr double kFreeTol = 1.0e-5 * kLength;
constexpr double kRelativeTol = 1.0e-5;
constexpr double kJReferenceTol = 1.0e-8;
constexpr double kZeroAbsTol = 1.0e-12;
const std::array<double, 4> kLoads{{1.0, 1.01, 1.02, 1.05}};

struct Reference {
    double transverse;
    double energy;
    double reaction;
    double J;
};

Reference referenceFor(double s) {
    const double mu = kYoungs / (2.0 * (1.0 + kPoisson));
    const double lambda = kYoungs * kPoisson / ((1.0 + kPoisson) * (1.0 - 2.0 * kPoisson));
    auto residual = [=](double t) { return mu * (t*t - 1.0) + lambda * std::log(s*t*t); };
    double lo = 0.5, hi = 1.5;
    if (residual(lo) * residual(hi) > 0.0) throw std::runtime_error("横向伸长参考方程未夹根");
    for (int i = 0; i < 100; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (residual(mid) > 0.0) hi = mid;
        else lo = mid;
    }
    const double t = 0.5 * (lo + hi);
    const double J = s * t * t;
    const double logJ = std::log(J);
    const double psi = 0.5 * mu * (s*s + 2.0*t*t - 3.0) - mu * logJ
                     + 0.5 * lambda * logJ * logJ;
    const double p11 = mu * (s - 1.0/s) + lambda * logJ / s;
    return {t, kVolume * psi, kArea * p11, J};
}

Eigen::Matrix3d deformationGradient(const Positions& x, const std::array<int,4>& tet,
                                    const Positions& x0) {
    Eigen::Matrix3d dm, ds;
    const Eigen::Vector3d p0 = x0.row(tet[0]).transpose();
    const Eigen::Vector3d q0 = x.row(tet[0]).transpose();
    for (int c = 0; c < 3; ++c) {
        dm.col(c) = x0.row(tet[c+1]).transpose() - p0;
        ds.col(c) = x.row(tet[c+1]).transpose() - q0;
    }
    return ds * dm.inverse();
}

struct Metrics {
    double endpointError = 0.0;
    double middleAffineError = 0.0;
    double freeAffineError = 0.0;
    double energy = 0.0;
    double reaction = 0.0;
    double energyRelativeError = 0.0;
    double reactionRelativeError = 0.0;
    double minJ = std::numeric_limits<double>::infinity();
    double maxJ = -std::numeric_limits<double>::infinity();
    double maxJReferenceError = 0.0;
    bool allJPositive = true;
};

Metrics measure(const MshMesh& mesh, const Positions& x, double s, const std::vector<int>& left,
                const std::vector<int>& right, const std::vector<int>& middle,
                const std::vector<char>& constrained, const Reference& ref,
                const NeoHookeanTetEnergy& energy) {
    Metrics m;
    auto affineError = [&](int i) {
        const double expected[3] = {s * mesh.x(i,0), ref.transverse * mesh.x(i,1),
                                    ref.transverse * mesh.x(i,2)};
        return std::max({std::abs(x(i,0)-expected[0]), std::abs(x(i,1)-expected[1]),
                         std::abs(x(i,2)-expected[2])});
    };
    for (int i : left) m.endpointError = std::max(m.endpointError, affineError(i));
    for (int i : right) m.endpointError = std::max(m.endpointError, affineError(i));
    for (int i : middle) m.middleAffineError = std::max(m.middleAffineError, affineError(i));
    for (int i=0; i<x.rows(); ++i) if (!constrained[i])
        m.freeAffineError = std::max(m.freeAffineError, affineError(i));
    m.energy = energy.val(x);
    Positions grad = Positions::Zero(x.rows(), 3);
    energy.addGrad(x, grad);
    for (int i : right) m.reaction += grad(i,0);
    for (const auto& tet : mesh.tets) {
        const double J = deformationGradient(x, tet, mesh.x).determinant();
        m.minJ = std::min(m.minJ, J);
        m.maxJ = std::max(m.maxJ, J);
        m.allJPositive = m.allJPositive && std::isfinite(J) && J > 0.0;
        m.maxJReferenceError = std::max(m.maxJReferenceError, std::abs(J-ref.J));
    }
    if (std::abs(ref.energy) > 1e-14) m.energyRelativeError = std::abs(m.energy-ref.energy)/std::abs(ref.energy);
    else m.energyRelativeError = std::abs(m.energy);
    if (std::abs(ref.reaction) > 1e-14) m.reactionRelativeError = std::abs(m.reaction-ref.reaction)/std::abs(ref.reaction);
    else m.reactionRelativeError = std::abs(m.reaction);
    return m;
}

void writeCsvHeader(std::ofstream& csv) {
    csv << "# protocol=fiction-validation-v0.2-S1\n"
        << "# geometry=L:1,b:0.05,h:0.05,rho:1000\n"
        << "# material=neo_hookean,E:1000000,nu:0.3\n"
        << "# solver=static,tol:1e-10,max_iter:100,dbc_stiffness_initial:1e8\n"
        << "# thresholds=endpoint_abs:1e-6*L,middle_affine_abs:1e-7*L,free_affine_abs:1e-5*L,"
           "energy_rel:1e-5,reaction_rel:1e-5,max_abs_J_error:1e-8,zero_abs:1e-12\n"
        << "load,transverse_reference,reference_energy,reference_reaction,iterations,line_search_iterations,"
           "converged,max_iter_hit,line_search_failed,search_direction_rate,free_gradient_inf,dbc_max_error,"
           "dbc_stiffness,dbc_satisfied,dbc_total,solver,endpoint_error,middle_affine_error,free_affine_error,"
           "energy,reaction,energy_error,reaction_error,min_J,max_J,max_abs_J_error,all_J_positive,pass\n";
}
}

int main(int argc, char** argv) {
    std::string meshPath, csvPath = "s1_metrics.csv";
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--mesh" && i + 1 < argc) meshPath = argv[++i];
            else if (arg == "--csv" && i + 1 < argc) csvPath = argv[++i];
            else if (arg == "--help") {
                std::cout << "用法: tension_check --mesh bar.msh [--csv metrics.csv]\n";
                return 0;
            } else throw std::runtime_error("未知或缺少值的参数: " + arg);
        }
        if (meshPath.empty()) throw std::runtime_error("必须提供 --mesh <path.msh>");
    } catch (const std::exception& e) {
        std::cerr << "S1 ERROR: " << e.what() << '\n';
        return 2;
    }

    try {
        const MshMesh mesh = loadMsh(meshPath);
        if (mesh.tets.empty() || mesh.x.rows() == 0) throw std::runtime_error("网格为空");
        const auto massModel = buildMassSpringModel(mesh, kDensity);
        if (std::abs(massModel.volume-kVolume) > 1e-12*kVolume ||
            std::abs(massModel.totalMass-kDensity*kVolume) > 1e-12*kDensity*kVolume)
            throw std::runtime_error("几何体积/集中总质量与 S1 规格不符");
        const auto bbox = massModel.bbox;
        if (std::abs(bbox[0]) > 1e-12 || std::abs(bbox[1]+0.5*kWidth) > 1e-12 ||
            std::abs(bbox[2]+0.5*kWidth) > 1e-12 || std::abs(bbox[3]-kLength) > 1e-12 ||
            std::abs(bbox[4]-0.5*kWidth) > 1e-12 || std::abs(bbox[5]-0.5*kWidth) > 1e-12)
            throw std::runtime_error("网格包围盒与 S1 方杆规格不符");

        std::vector<int> left, right, middle;
        std::vector<char> constrained(mesh.x.rows(), 0);
        for (int i=0; i<mesh.x.rows(); ++i) {
            if (std::abs(mesh.x(i,0)) < 1e-12) left.push_back(i);
            if (std::abs(mesh.x(i,0)-kLength) < 1e-12) right.push_back(i);
            if (mesh.x(i,0) >= 0.2*kLength-1e-12 && mesh.x(i,0) <= 0.8*kLength+1e-12)
                middle.push_back(i);
        }
        if (left.empty() || right.empty() || middle.empty()) throw std::runtime_error("未找到完整端面或中截面节点");

        DBCInfo dbc;
        dbc.stiffness = kDbcStiffness;
        dbc.isDBC.assign(mesh.x.rows(), 0);
        for (int i : left) {
            dbc.indices.push_back(i); dbc.targetVelocities.push_back(Vec3::Zero());
            dbc.limits.push_back(Vec3::Zero()); dbc.isDBC[i] = 1;
            dbc.targetPositions.push_back(mesh.x.row(i).transpose()); constrained[i] = 1;
        }
        // 右端目标随加载步改变；左右端面完整固定在给定仿射目标上。
        for (int i : right) {
            dbc.indices.push_back(i); dbc.targetVelocities.push_back(Vec3::Zero());
            dbc.limits.push_back(Vec3::Zero()); dbc.isDBC[i] = 1; constrained[i] = 1;
            dbc.targetPositions.push_back(mesh.x.row(i).transpose());
        }

        Eigen::VectorXd masses = massModel.masses;
        TimeIntegrator::Config cfg;
        cfg.h = 1.0; cfg.tol = kTolerance; cfg.maxIter = 100; cfg.verbose = false;
        cfg.staticSolve = true;
        TimeIntegrator solver(cfg, masses);
        auto energy = std::make_shared<NeoHookeanTetEnergy>(mesh.tets, mesh.x, kYoungs, kPoisson);
        solver.addEnergy(energy);
        std::ofstream csv(csvPath);
        if (!csv) throw std::runtime_error("无法创建 CSV: " + csvPath);
        csv << std::setprecision(17)
            << "# mesh_path=" << meshPath << ",nodes=" << mesh.x.rows()
            << ",tets=" << mesh.tets.size() << ",volume=" << massModel.volume
            << ",total_mass=" << massModel.totalMass << '\n';
        writeCsvHeader(csv);
        if (!csv) throw std::runtime_error("写入 CSV header 失败: " + csvPath);

        Positions x = mesh.x;
        Velocities v = Velocities::Zero(mesh.x.rows(),3);
        std::cout << std::setprecision(10) << "S1 mesh=" << meshPath << " nodes=" << mesh.x.rows()
                  << " tets=" << mesh.tets.size() << " volume=" << massModel.volume
                  << " mass=" << massModel.totalMass << '\n';
        for (double s : kLoads) {
            const Reference ref = referenceFor(s);
            for (size_t q=0; q<left.size(); ++q) {
                const int i = left[q];
                dbc.targetPositions[q] = Vec3(mesh.x(i,0), ref.transverse*mesh.x(i,1),
                                              ref.transverse*mesh.x(i,2));
            }
            for (size_t q=0; q<right.size(); ++q) {
                const int i = right[q];
                dbc.targetPositions[left.size()+q] = Vec3(s*mesh.x(i,0), ref.transverse*mesh.x(i,1),
                                                           ref.transverse*mesh.x(i,2));
            }
            const auto state = solver.stepForward(x, v, dbc);
            x = state.x;
            v = state.v;

            Metrics m = measure(mesh, x, s, left, right, middle, constrained, ref, *energy);
            double maxDisplacement = 0.0;
            for (int i=0; i<x.rows(); ++i)
                maxDisplacement = std::max(maxDisplacement, (x.row(i)-mesh.x.row(i)).cwiseAbs().maxCoeff());

            const bool zero = (s == 1.0);
            const bool metricsFinite = std::isfinite(ref.transverse) && std::isfinite(ref.energy) &&
                std::isfinite(ref.reaction) && std::isfinite(m.endpointError) &&
                std::isfinite(m.middleAffineError) && std::isfinite(m.freeAffineError) &&
                std::isfinite(m.energy) && std::isfinite(m.reaction) &&
                std::isfinite(m.energyRelativeError) && std::isfinite(m.reactionRelativeError) &&
                std::isfinite(m.minJ) && std::isfinite(m.maxJ) && std::isfinite(m.maxJReferenceError) &&
                std::isfinite(state.report.searchDirectionRate) && std::isfinite(state.report.freeGradientInf) &&
                std::isfinite(state.report.dbcMaxError) && std::isfinite(state.report.dbcStiffness) &&
                state.x.allFinite() && state.v.allFinite();
            const bool pass = metricsFinite && state.report.converged && !state.report.maxIterHit &&
                !state.report.lineSearchFailed && m.endpointError < kEndpointTol &&
                state.report.dbcSatisfied == state.report.dbcTotal &&
                m.middleAffineError < kMiddleTol && m.freeAffineError < kFreeTol &&
                (zero ? (maxDisplacement < kZeroAbsTol && std::abs(m.energy) < kZeroAbsTol && std::abs(m.reaction) < kZeroAbsTol)
                      : (m.energyRelativeError < kRelativeTol && m.reactionRelativeError < kRelativeTol)) &&
                m.allJPositive && m.maxJReferenceError < kJReferenceTol && state.x.allFinite();

            csv << s << ',' << ref.transverse << ',' << ref.energy << ',' << ref.reaction << ','
                << state.report.iterations << ',' << state.report.lineSearchIterations << ','
                << state.report.converged << ',' << state.report.maxIterHit << ',' << state.report.lineSearchFailed << ','
                << state.report.searchDirectionRate << ',' << state.report.freeGradientInf << ','
                << state.report.dbcMaxError << ',' << state.report.dbcStiffness << ','
                << state.report.dbcSatisfied << ',' << state.report.dbcTotal << ',' << state.report.solver << ','
                << m.endpointError << ',' << m.middleAffineError << ',' << m.freeAffineError << ','
                << m.energy << ',' << m.reaction << ',' << m.energyRelativeError << ','
                << m.reactionRelativeError << ',' << m.minJ << ',' << m.maxJ << ','
                << m.maxJReferenceError << ',' << m.allJPositive << ',' << pass << '\n';
            csv.flush();
            if (!csv) throw std::runtime_error("写入 CSV 失败: " + csvPath);
            std::cout << "load=" << s << " t=" << ref.transverse << " energy=" << m.energy
                      << " reaction=" << m.reaction << " affine_free=" << m.freeAffineError
                      << " J=[" << m.minJ << ',' << m.maxJ << "] iter=" << state.report.iterations
                      << " converged=" << state.report.converged << " solver=" << state.report.solver
                      << (pass ? " PASS\n" : " FAIL\n");
            if (!pass) return 1; // 第一处失败即停，不自动调参或继续载荷。
        }
        std::cout << "S1 PASS; metrics=" << csvPath << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "S1 ERROR: " << e.what() << '\n';
        return 2;
    }
}
