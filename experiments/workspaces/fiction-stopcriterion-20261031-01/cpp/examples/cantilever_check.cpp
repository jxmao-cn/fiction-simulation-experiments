// S2：小变形自重悬臂。只判定单网格数值求解，精度及单调性由批次汇总判断。
#include "energies/Gravity.h"
#include "energies/NeoHookeanTetEnergy.h"
#include "io/msh_loader.h"
#include "scene/mass_spring_builder.h"
#include "time_integrator.h"
#include <Eigen/SVD>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
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
constexpr double E = 5.90e8, nu = 0.3, rho = 1000.0, g = 9.81;
constexpr double L = 1.0, b = 0.05, h = 0.05, A = b*h, V = L*A;
constexpr double tol = 1e-10, initialK = 1e8;
Eigen::Matrix3d edges(const Positions& x, const std::array<int,4>& t) {
    Eigen::Matrix3d d;
    for (int c=0; c<3; ++c) d.col(c) = (x.row(t[c+1])-x.row(t[0])).transpose();
    return d;
}
void checkedClose(std::ofstream& out, const std::string& path) {
    out.close();
    if (!out) throw std::runtime_error("写入失败: " + path);
}
}

int main(int argc, char** argv) {
    try {
        std::string meshPath, csvPath, statePath, traceDir;
        // 停止判据诊断参数（默认 step，行为与原版一致）
        std::string stopMode = "step";
        double gradTol = 1e-7;
        double gradScale = std::numeric_limits<double>::quiet_NaN(); // NaN -> 用总重 24.525 N
        double actualStepTol = 1e-12;
        int actualStepCount = 2;
        for (int i=1; i<argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help") {
                std::cout << "cantilever_check --mesh FILE --csv FILE --state FILE [--trace-dir DIR]\n"
                             "  [--stop-mode step|grad|actual|combined] [--grad-tol 1e-7]\n"
                             "  [--grad-scale 24.525] [--actual-step-tol 1e-12] [--actual-step-count 2]\n";
                return 0;
            }
            if (i+1 >= argc) throw std::runtime_error("缺少参数值: " + arg);
            if (arg == "--mesh") meshPath = argv[++i];
            else if (arg == "--csv") csvPath = argv[++i];
            else if (arg == "--state") statePath = argv[++i];
            else if (arg == "--trace-dir") traceDir = argv[++i];
            else if (arg == "--stop-mode") stopMode = argv[++i];
            else if (arg == "--grad-tol") gradTol = std::stod(argv[++i]);
            else if (arg == "--grad-scale") gradScale = std::stod(argv[++i]);
            else if (arg == "--actual-step-tol") actualStepTol = std::stod(argv[++i]);
            else if (arg == "--actual-step-count") actualStepCount = std::stoi(argv[++i]);
            else throw std::runtime_error("未知参数: " + arg);
        }
        if (meshPath.empty() || csvPath.empty() || statePath.empty())
            throw std::runtime_error("必须指定 --mesh、--csv 和 --state");
        if (std::filesystem::path(csvPath).lexically_normal() == std::filesystem::path(statePath).lexically_normal())
            throw std::runtime_error("CSV 和状态路径必须不同");
        const MshMesh mesh = loadMsh(meshPath);
        if (mesh.empty() || mesh.tets.empty() || !mesh.x.allFinite()) throw std::runtime_error("空网格或非有限坐标");
        const auto model = buildMassSpringModel(mesh, rho);
        const std::array<double,6> box{{0.0,-b/2,-h/2,L,b/2,h/2}};
        for (int k=0; k<6; ++k) if (std::abs(model.bbox[k]-box[k]) > 1e-12)
            throw std::runtime_error("包围盒与 S2 规格不符");
        if (std::abs(model.volume-V) > 1e-12*V || std::abs(model.totalMass-rho*V) > 1e-12*rho*V)
            throw std::runtime_error("体积或质量不符");
        for (const auto& t : mesh.tets) if (!(edges(mesh.x,t).determinant() > 0.0))
            throw std::runtime_error("参考四面体非正定向");

        std::vector<int> tip;
        DBCInfo dbc;
        dbc.stiffness = initialK;
        dbc.isDBC.assign(mesh.x.rows(), 0);
        for (int i=0; i<mesh.x.rows(); ++i) {
            if (std::abs(mesh.x(i,0)) < 1e-12) {
                dbc.indices.push_back(i); dbc.targetPositions.push_back(mesh.x.row(i).transpose());
                dbc.targetVelocities.push_back(Vec3::Zero()); dbc.limits.push_back(Vec3::Zero());
                dbc.isDBC[i] = 1;
            }
            if (std::abs(mesh.x(i,0)-L) < 1e-12) tip.push_back(i);
        }
        if (dbc.indices.empty() || tip.empty()) throw std::runtime_error("缺少端面节点");
        // 参考端面上，P1 位移的面积积分恰等于三角形三节点平均乘面积。
        Eigen::VectorXd tipWeights = Eigen::VectorXd::Zero(mesh.x.rows());
        double tipArea = 0.0;
        for (const auto& f : mesh.surface) {
            bool end = true;
            for (int i : f) end = end && std::abs(mesh.x(i,0)-L) < 1e-12;
            if (!end) continue;
            const Vec3 u = (mesh.x.row(f[1])-mesh.x.row(f[0])).transpose();
            const Vec3 v = (mesh.x.row(f[2])-mesh.x.row(f[0])).transpose();
            const double area = 0.5*u.cross(v).norm();
            if (!(area > 0.0) || !std::isfinite(area)) throw std::runtime_error("端面三角形退化");
            tipArea += area;
            for (int i : f) tipWeights[i] += area/3.0;
        }
        if (std::abs(tipArea-A) > 1e-12*A) throw std::runtime_error("端面面积不符");
        for (int i : tip) if (!(tipWeights[i] > 0.0)) throw std::runtime_error("端节点未被面积积分覆盖");

        std::ofstream csv(csvPath), stateOut(statePath);
        if (!csv || !stateOut) throw std::runtime_error("无法创建输出文件");
        csv << std::setprecision(17);
        stateOut << std::setprecision(17);
        TimeIntegrator::Config cfg;
        cfg.h = 1.0; cfg.tol = tol; cfg.maxIter = 100; cfg.staticSolve = true; cfg.verbose = false;
        // 停止判据诊断：gradScale 缺省取总重（= model.totalMass * g = 24.525 N）
        if (!std::isfinite(gradScale)) gradScale = model.totalMass * g;
        cfg.stopMode = stopMode;
        cfg.gradTolRel = gradTol;
        cfg.gradScale = gradScale;
        cfg.actualStepTol = actualStepTol;
        cfg.actualStepCount = actualStepCount;
        if (!traceDir.empty()) {
            std::filesystem::create_directories(traceDir);
            cfg.traceCsvPath = (std::filesystem::path(traceDir)/"iterations.csv").string();
            cfg.traceLabel = "S2";
        }
        auto elastic = std::make_shared<NeoHookeanTetEnergy>(mesh.tets,mesh.x,E,nu);
        auto gravity = std::make_shared<ExampleGravityEnergy>(model.masses,Vec3(0,0,-g));
        TimeIntegrator solver(cfg,model.masses);
        solver.addEnergy(elastic); solver.addEnergy(gravity);
        const auto begin = std::chrono::steady_clock::now();
        const auto state = solver.stepForward(mesh.x,Velocities::Zero(mesh.x.rows(),3),dbc);
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
        const auto& x = state.x;
        const auto& r = state.report;
        double displacement = 0.0, arithmetic = 0.0;
        double tipMin = std::numeric_limits<double>::infinity(), tipMax = -tipMin;
        for (int i : tip) {
            const double down = mesh.x(i,2)-x(i,2);
            displacement += tipWeights[i]*down/tipArea;
            arithmetic += down/static_cast<double>(tip.size());
            tipMin = std::min(tipMin,down); tipMax = std::max(tipMax,down);
        }
        double dbcError = 0.0;
        for (int i : dbc.indices) dbcError = std::max(dbcError,(x.row(i)-mesh.x.row(i)).cwiseAbs().maxCoeff());
        double minJ = std::numeric_limits<double>::infinity(), maxJ = -minJ;
        double minStretch = minJ, maxStretch = -minJ;
        bool finite = x.allFinite() && state.v.allFinite();
        for (const auto& t : mesh.tets) {
            const Eigen::Matrix3d F = edges(x,t)*edges(mesh.x,t).inverse();
            const double J = F.determinant();
            finite = finite && F.allFinite() && std::isfinite(J);
            minJ = std::min(minJ,J); maxJ = std::max(maxJ,J);
            if (F.allFinite()) {
                const Eigen::Vector3d stretch = Eigen::JacobiSVD<Eigen::Matrix3d>(F).singularValues();
                minStretch = std::min(minStretch,stretch.minCoeff());
                maxStretch = std::max(maxStretch,stretch.maxCoeff());
            }
        }
        Positions grad = Positions::Zero(mesh.x.rows(),3);
        elastic->addGrad(x,grad); gravity->addGrad(x,grad);
        // 总势能对夹持节点的导数即维持该约束所需的外部支反力，包含夹持节点自重。
        Vec3 reaction = Vec3::Zero();
        for (int i : dbc.indices) reaction += grad.row(i).transpose();
        const Vec3 external(0,0,-model.totalMass*g);
        const double I = b*h*h*h/12.0, G = E/(2.0*(1.0+nu)), kappa = 5.0/6.0, w = rho*A*g;
        const double bend = w*std::pow(L,4)/(8.0*E*I);
        const double shear = w*L*L/(2.0*kappa*G*A), reference = bend+shear;
        const double error = std::abs(arithmetic-reference)/reference;
        const double elasticEnergy = elastic->val(x), gravityEnergy = gravity->val(x);
        for (double value : {displacement,arithmetic,tipMin,tipMax,dbcError,minJ,maxJ,minStretch,maxStretch,
                             error,elasticEnergy,gravityEnergy,wall,r.searchDirectionRate,r.freeGradientInf,
                             r.dbcMaxError,r.dbcStiffness}) finite = finite && std::isfinite(value);
        finite = finite && grad.allFinite() && reaction.allFinite();
        const bool pass = finite && r.converged && !r.maxIterHit && !r.lineSearchFailed && minJ > 0.9 &&
                          r.dbcSatisfied == r.dbcTotal && dbcError < 1e-6;
        csv << "# protocol=fiction-validation-v0.2-S2; E=590000000,nu=0.3,rho=1000,gz=-9.81\n"
            << "# static=true,tol=1e-10,max_iter=100,initial_dbc_stiffness=1e8; tip=tip_face_node_arithmetic_mean_downward_displacement\n"
            << "# stopcriterion-diagnostic: stop_mode=" << stopMode << ",grad_tol_rel=" << gradTol
            << ",grad_scale=" << gradScale << ",actual_step_tol=" << actualStepTol
            << ",actual_step_count=" << actualStepCount << "\n"
            << "nodes,tets,volume,total_mass,clamp_nodes,tip_nodes,tip_area,reference_bending,reference_shear,reference_tip,tip_area_weighted,tip_deflection,tip_min,tip_max,relative_error,iterations,line_search_iterations,converged,max_iter_hit,line_search_failed,search_direction_rate,free_gradient_inf,dbc_max_error,measured_dbc_error,dbc_stiffness,dbc_satisfied,dbc_total,solver,min_J,max_J,min_stretch,max_stretch,elastic_energy,gravity_energy,external_fx,external_fy,external_fz,reaction_x,reaction_y,reaction_z,force_balance_norm,solve_wall_seconds,numerical_pass,stop_mode,stop_reason,grad_scale,grad_tol_rel,free_gradient_relative,actual_step_tol,actual_step_count\n"
            << x.rows() << ',' << mesh.tets.size() << ',' << model.volume << ',' << model.totalMass << ','
            << dbc.indices.size() << ',' << tip.size() << ',' << tipArea << ',' << bend << ',' << shear << ','
            << reference << ',' << displacement << ',' << arithmetic << ',' << tipMin << ',' << tipMax << ','
            << error << ',' << r.iterations << ',' << r.lineSearchIterations << ',' << r.converged << ','
            << r.maxIterHit << ',' << r.lineSearchFailed << ',' << r.searchDirectionRate << ',' << r.freeGradientInf << ','
            << r.dbcMaxError << ',' << dbcError << ',' << r.dbcStiffness << ',' << r.dbcSatisfied << ',' << r.dbcTotal << ','
            << r.solver << ',' << minJ << ',' << maxJ << ',' << minStretch << ',' << maxStretch << ','
            << elasticEnergy << ',' << gravityEnergy << ',' << external.x() << ',' << external.y() << ',' << external.z() << ','
            << reaction.x() << ',' << reaction.y() << ',' << reaction.z() << ',' << (reaction+external).norm() << ',' << wall << ',' << pass
            << ',' << stopMode << ',' << r.stopReason << ',' << gradScale << ',' << gradTol << ','
            << (r.freeGradientInf / gradScale) << ',' << actualStepTol << ',' << actualStepCount << '\n';
        stateOut << "node,X,Y,Z,x,y,z,ux,uy,uz,is_dbc,tip_area_weight\n";
        for (int i=0; i<x.rows(); ++i) {
            stateOut << i;
            for (int d=0; d<3; ++d) stateOut << ',' << mesh.x(i,d);
            for (int d=0; d<3; ++d) stateOut << ',' << x(i,d);
            for (int d=0; d<3; ++d) stateOut << ',' << x(i,d)-mesh.x(i,d);
            stateOut << ',' << int(dbc.isDBC[i]) << ',' << tipWeights[i] << '\n';
        }
        checkedClose(csv,csvPath); checkedClose(stateOut,statePath);
        std::cout << std::setprecision(12) << "S2 nodes=" << x.rows() << " tets=" << mesh.tets.size()
                  << " tip=" << arithmetic << " reference=" << reference << " error=" << error
                  << " minJ=" << minJ << " iterations=" << r.iterations << " solver=" << r.solver
                  << " stop=" << r.stopReason
                  << (pass ? " NUMERICAL_PASS\n" : " NUMERICAL_FAIL\n");
        return pass ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "S2 ERROR: " << e.what() << '\n';
        return 2;
    }
}
