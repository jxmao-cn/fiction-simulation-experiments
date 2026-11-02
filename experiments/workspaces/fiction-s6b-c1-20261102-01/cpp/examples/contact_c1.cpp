// =============================================================================
// contact_c1.cpp — S6b C1：sphere1K 对解析平面的无摩擦接触 pilot/正式入口
// =============================================================================
// 隔离实验入口；不修改生产默认路径。数值配置由 CLI 冻结并写入 config.json。
// 指标：接触数、最小距离、CCD/非候选步长、Newton/线搜索、min J、接触反力、
//       重力/弹性能、墙钟与峰值 RSS。
// =============================================================================

#include "energies/ContactEnergy.h"
#include "energies/Gravity.h"
#include "energies/NeoHookeanTetEnergy.h"
#include "io/msh_loader.h"
#include "io/obj_writer.h"
#include "scene/mass_spring_builder.h"
#include "time_integrator.h"
#include "types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace fiction;

namespace {

struct Options {
    std::string mesh;
    std::string out;
    int steps = 200;
    double dt = 2e-3;
    double dhat = 1.7e-3;
    double tol = 1e-4;
    int maxIter = 50;
    double youngs = 1e6;
    double poisson = 0.3;
    double density = 1000.0;
    double kappa = 0.0;      // <=0: t=0 heuristic
    double gravity = -9.81;
    int exportEvery = 0;
    bool verbose = false;
};

double argDouble(const std::string& s) {
    size_t p = 0;
    const double v = std::stod(s, &p);
    if (p != s.size()) throw std::runtime_error("invalid number: " + s);
    return v;
}

int argInt(const std::string& s) {
    size_t p = 0;
    const long v = std::stol(s, &p);
    if (p != s.size() || v < 0) throw std::runtime_error("invalid int: " + s);
    return static_cast<int>(v);
}

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&](const char* flag) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string(flag) + " missing value");
            return argv[++i];
        };
        if (a == "--mesh") o.mesh = next("--mesh");
        else if (a == "--out") o.out = next("--out");
        else if (a == "--steps") o.steps = argInt(next("--steps"));
        else if (a == "--dt") o.dt = argDouble(next("--dt"));
        else if (a == "--dhat") o.dhat = argDouble(next("--dhat"));
        else if (a == "--tol") o.tol = argDouble(next("--tol"));
        else if (a == "--max-iter") o.maxIter = argInt(next("--max-iter"));
        else if (a == "--youngs") o.youngs = argDouble(next("--youngs"));
        else if (a == "--poisson") o.poisson = argDouble(next("--poisson"));
        else if (a == "--density") o.density = argDouble(next("--density"));
        else if (a == "--kappa") o.kappa = argDouble(next("--kappa"));
        else if (a == "--gravity") o.gravity = argDouble(next("--gravity"));
        else if (a == "--export-every") o.exportEvery = argInt(next("--export-every"));
        else if (a == "--verbose") o.verbose = true;
        else throw std::runtime_error("unknown option: " + a);
    }
    if (o.mesh.empty() || o.out.empty()) throw std::runtime_error("--mesh and --out are required");
    if (o.steps <= 0 || !(o.dt > 0) || !(o.dhat > 0) || !(o.tol > 0) || o.maxIter < 0)
        throw std::runtime_error("invalid numeric option");
    return o;
}

Eigen::MatrixXi surfaceMatrix(const MshMesh& mesh) {
    Eigen::MatrixXi F(mesh.numSurface(), 3);
    for (int i = 0; i < mesh.numSurface(); ++i)
        for (int d = 0; d < 3; ++d) F(i, d) = mesh.surface[i][d];
    return F;
}

struct TetJacobian {
    Eigen::Matrix3d invDm;
    bool valid = false;
};

std::vector<TetJacobian> buildTetJacobians(const MshMesh& mesh, const Positions& rest) {
    std::vector<TetJacobian> out(mesh.tets.size());
    for (size_t t = 0; t < mesh.tets.size(); ++t) {
        const auto& tet = mesh.tets[t];
        const Vec3 p0 = rest.row(tet[0]).transpose();
        Eigen::Matrix3d Dm;
        Dm.col(0) = rest.row(tet[1]).transpose() - p0;
        Dm.col(1) = rest.row(tet[2]).transpose() - p0;
        Dm.col(2) = rest.row(tet[3]).transpose() - p0;
        if (std::abs(Dm.determinant()) > 1e-18) {
            out[t].invDm = Dm.inverse();
            out[t].valid = true;
        }
    }
    return out;
}

std::pair<double, double> minMaxJ(const MshMesh& mesh, const Positions& x,
                                  const std::vector<TetJacobian>& jac) {
    double lo = std::numeric_limits<double>::infinity();
    double hi = -std::numeric_limits<double>::infinity();
    for (size_t t = 0; t < mesh.tets.size(); ++t) {
        if (!jac[t].valid) continue;
        const auto& tet = mesh.tets[t];
        const Vec3 p0 = x.row(tet[0]).transpose();
        Eigen::Matrix3d Ds;
        Ds.col(0) = x.row(tet[1]).transpose() - p0;
        Ds.col(1) = x.row(tet[2]).transpose() - p0;
        Ds.col(2) = x.row(tet[3]).transpose() - p0;
        const double J = (Ds * jac[t].invDm).determinant();
        lo = std::min(lo, J);
        hi = std::max(hi, J);
    }
    return {lo, hi};
}

long peakRssKiB() {
    std::ifstream in("/proc/self/status");
    std::string key;
    while (in >> key) {
        if (key == "VmHWM:") {
            long v = 0; in >> v; return v;
        }
        std::getline(in, key);
    }
    return -1;
}

double mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0; for (double x : v) s += x; return s / v.size();
}

double stddev(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    const double m = mean(v);
    double s = 0.0; for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / (v.size() - 1));
}

} // namespace

int main(int argc, char** argv) {
    Options opt;
    try {
        opt = parse(argc, argv);
        std::filesystem::create_directories(opt.out);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 2;
    }

    try {
        MshMesh mesh = loadMsh(opt.mesh);
        if (mesh.numTets() == 0 || mesh.numSurface() == 0)
            throw std::runtime_error("mesh must contain tets and surface triangles");

        MassSpringModel model = buildMassSpringModel(mesh, opt.density);
        const int N = model.numNodes();
        const Eigen::MatrixXi faces = surfaceMatrix(mesh);

        // 初始几何：最低表面点位于 2*dhat；整体平移不改变材料/接触能量。
        double minZ = std::numeric_limits<double>::infinity();
        for (int i = 0; i < faces.rows(); ++i)
            for (int d = 0; d < 3; ++d) minZ = std::min(minZ, mesh.x(faces(i, d), 2));
        const double dz = 2.0 * opt.dhat - minZ;
        Positions rest = model.x0;
        rest.col(2).array() += dz;
        mesh.x.col(2).array() += dz;

        const double bboxDiag = std::hypot(
            model.bbox[3] - model.bbox[0],
            std::hypot(model.bbox[4] - model.bbox[1], model.bbox[5] - model.bbox[2]));

        auto gravity = std::make_shared<ExampleGravityEnergy>(
            model.masses, Vec3(0.0, 0.0, opt.gravity));
        auto material = std::make_shared<NeoHookeanTetEnergy>(
            mesh.tets, rest, opt.youngs, opt.poisson);

        ContactEnergy::Options cop;
        cop.dhat = opt.dhat;
        cop.dmin = 0.0;
        cop.kappa = opt.kappa;
        cop.useAreaWeighting = false;
        auto contact = std::make_shared<ContactEnergy>(
            rest, faces, Eigen::VectorXi::Zero(0), cop);
        contact->addPlane(Eigen::Vector3d(0, 0, 1), 0.0);

        // t=0 的弹性+重力梯度；接触初始无活跃对，grad_barrier 为 0。
        Positions ge = Positions::Zero(N, kDim);
        gravity->addGrad(rest, ge);
        material->addGrad(rest, ge);
        const double kappa = contact->initKappa(rest, bboxDiag, model.masses.mean(), flatten(ge));

        const auto jac = buildTetJacobians(mesh, rest);
        auto [minJ0, maxJ0] = minMaxJ(mesh, rest, jac);

        TimeIntegrator::Config cfg;
        cfg.h = opt.dt;
        cfg.tol = opt.tol;
        cfg.maxIter = opt.maxIter;
        cfg.verbose = opt.verbose;
        cfg.stopMode = "step";
        TimeIntegrator integrator(cfg, model.masses);
        integrator.addEnergy(gravity);
        integrator.addEnergy(material);
        integrator.addEnergy(contact);

        DBCInfo dbc; // 无 DBC
        Positions x = rest;
        Velocities v = Velocities::Zero(N, kDim);

        const double totalMass = model.masses.sum();
        const double weight = totalMass * std::abs(opt.gravity);

        std::ostringstream csv;
        csv << "step,t,contacts,candidates,d_min,alpha_ccd,alpha_nc,iterations,line_search,"
               "converged,stop_reason,min_J,max_J,Fz,E_gravity,E_elastic,E_contact,wall_s\n";
        csv << std::setprecision(17);

        double wallTotal = 0.0;
        int gateFailed = -1;
        std::string gateReason;
        for (int s = 1; s <= opt.steps; ++s) {
            const auto t0 = std::chrono::steady_clock::now();
            TimeIntegrator::State st = integrator.stepForward(x, v, dbc);
            const double wall = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
            wallTotal += wall;
            x = st.x; v = st.v;

            const auto [minJ, maxJ] = minMaxJ(mesh, x, jac);
            const Eigen::VectorXd f = contact->forces(x);
            double Fz = 0.0;
            for (int i = 0; i < N; ++i) Fz += f(3 * i + 2);

            csv << s << ',' << s * opt.dt << ',' << contact->contactCount() << ','
                << contact->candidateCount() << ',' << contact->minDistance() << ','
                << contact->lastAlphaCcd() << ',' << contact->lastAlphaNc() << ','
                << st.report.iterations << ',' << st.report.lineSearchIterations << ','
                << (st.report.converged ? 1 : 0) << ',' << st.report.stopReason << ','
                << minJ << ',' << maxJ << ',' << Fz << ','
                << gravity->val(x) << ',' << material->val(x) << ',' << contact->val(x) << ','
                << wall << '\n';

            if (!st.report.converged) {
                gateFailed = s; gateReason = "not_converged:" + st.report.stopReason;
                std::cerr << "step " << s << " failed: reason=" << st.report.stopReason
                          << " iter=" << st.report.iterations
                          << " ls=" << st.report.lineSearchIterations
                          << " alpha_ccd=" << contact->lastAlphaCcd()
                          << " alpha_nc=" << contact->lastAlphaNc()
                          << " dmin=" << contact->minDistance()
                          << " contacts=" << contact->contactCount()
                          << " minJ=" << minJ << "\n";
                break;
            }
            if (!x.allFinite() || !v.allFinite()) { gateFailed = s; gateReason = "nonfinite"; break; }
            if (!(minJ > 0.3)) { gateFailed = s; gateReason = "minJ<=0.3"; break; }
            if (!(contact->minDistance() > 0.0)) { gateFailed = s; gateReason = "d<=0"; break; }
            if (contact->contactCount() > 2000) { gateFailed = s; gateReason = "contacts>2000"; break; }
            if (opt.exportEvery > 0 && s % opt.exportEvery == 0) {
                std::ostringstream p;
                p << opt.out << "/frame_" << std::setw(5) << std::setfill('0') << s << ".obj";
                writeSurfaceObj(p.str(), x, mesh.surface);
            }
        }

        {
            std::ofstream out(opt.out + "/steps.csv");
            out << csv.str();
        }

        const Positions finalX = x;
        const auto [minJf, maxJf] = minMaxJ(mesh, finalX, jac);
        writeSurfaceObj(opt.out + "/final.obj", finalX, mesh.surface);

        // 从 steps.csv 重读指标，避免循环内提前停止时统计口径不一致。
        std::vector<double> dmin, fz;
        int lastStep = 0;
        bool allConv = true;
        {
            std::ifstream in(opt.out + "/steps.csv");
            std::string line;
            std::getline(in, line);
            while (std::getline(in, line)) {
                std::stringstream ss(line);
                std::string tok;
                std::vector<std::string> cols;
                while (std::getline(ss, tok, ',')) cols.push_back(tok);
                if (cols.size() < 18) continue;
                const int step = std::stoi(cols[0]);
                lastStep = step;
                allConv = allConv && (std::stoi(cols[9]) == 1);
                dmin.push_back(std::stod(cols[4]));
                fz.push_back(std::stod(cols[13]));
            }
        }
        const int steadyN = std::max(1, static_cast<int>(0.2 / opt.dt));
        const int start = std::max(0, static_cast<int>(dmin.size()) - steadyN);
        const std::vector<double> dminSteady(dmin.begin() + start, dmin.end());
        const std::vector<double> fzSteady(fz.begin() + start, fz.end());
        const double dminMean = mean(dminSteady);
        const double fzMean = mean(fzSteady);
        const double fzRel = std::abs(fzMean - weight) / weight;
        const double fzStdRel = mean(fzSteady) > 0 ? stddev(fzSteady) / std::abs(mean(fzSteady)) : 1.0;
        bool dPositive = true;
        for (double d : dmin) dPositive = dPositive && d > 0.0;

        const bool pass = (gateFailed < 0) && allConv && (lastStep == opt.steps) &&
                          dPositive && (dminMean > 0.1 * opt.dhat) &&
                          (fzRel < 0.05) && (fzStdRel < 0.05) &&
                          (minJf > 0.3) && (dminMean < opt.dhat) &&
                          (contact->contactCount() <= 2000);

        std::ofstream meta(opt.out + "/config.json");
        meta << std::setprecision(17)
             << "{\n"
             << "  \"run_id\": \"s6b-c1-plane-20261102-01\",\n"
             << "  \"mesh\": \"" << opt.mesh << "\",\n"
             << "  \"nodes\": " << N << ",\n"
             << "  \"tets\": " << mesh.numTets() << ",\n"
             << "  \"surface_faces\": " << mesh.numSurface() << ",\n"
             << "  \"collision_vertices\": " << contact->collisionVertices() << ",\n"
             << "  \"surface_edges\": " << contact->surfaceEdges() << ",\n"
             << "  \"steps\": " << opt.steps << ",\n"
             << "  \"dt\": " << opt.dt << ",\n"
             << "  \"dhat\": " << opt.dhat << ",\n"
             << "  \"dmin\": 0.0,\n"
             << "  \"tol\": " << opt.tol << ",\n"
             << "  \"max_iter\": " << opt.maxIter << ",\n"
             << "  \"youngs\": " << opt.youngs << ",\n"
             << "  \"poisson\": " << opt.poisson << ",\n"
             << "  \"density\": " << opt.density << ",\n"
             << "  \"gravity\": " << opt.gravity << ",\n"
             << "  \"initial_shift_z\": " << dz << ",\n"
             << "  \"initial_min_J\": " << minJ0 << ",\n"
             << "  \"initial_max_J\": " << maxJ0 << ",\n"
             << "  \"kappa\": " << kappa << ",\n"
             << "  \"kappa_max\": " << contact->kappaMax() << ",\n"
             << "  \"total_mass\": " << totalMass << ",\n"
             << "  \"weight\": " << weight << ",\n"
             << "  \"wall_total_s\": " << wallTotal << ",\n"
             << "  \"peak_rss_kib\": " << peakRssKiB() << "\n"
             << "}\n";

        std::ofstream summary(opt.out + "/summary.json");
        summary << std::setprecision(17)
                << "{\n"
                << "  \"status\": \"" << (pass ? "passed" : "failed") << "\",\n"
                << "  \"steps_completed\": " << lastStep << ",\n"
                << "  \"gate_failed_step\": " << gateFailed << ",\n"
                << "  \"gate_reason\": \"" << gateReason << "\",\n"
                << "  \"all_converged\": " << (allConv ? "true" : "false") << ",\n"
                << "  \"d_min_positive\": " << (dPositive ? "true" : "false") << ",\n"
                << "  \"d_min_steady_mean\": " << dminMean << ",\n"
                << "  \"d_min_steady_gt_0.1dhat\": " << (dminMean > 0.1 * opt.dhat ? "true" : "false") << ",\n"
                << "  \"d_min_steady_lt_dhat\": " << (dminMean < opt.dhat ? "true" : "false") << ",\n"
                << "  \"Fz_steady_mean\": " << fzMean << ",\n"
                << "  \"weight\": " << weight << ",\n"
                << "  \"Fz_rel_error\": " << fzRel << ",\n"
                << "  \"Fz_std_over_mean\": " << fzStdRel << ",\n"
                << "  \"min_J_final\": " << minJf << ",\n"
                << "  \"contact_count_final\": " << contact->contactCount() << ",\n"
                << "  \"wall_total_s\": " << wallTotal << ",\n"
                << "  \"peak_rss_kib\": " << peakRssKiB() << "\n"
                << "}\n";

        std::cout << "C1 " << (pass ? "PASS" : "FAIL")
                  << " steps=" << lastStep << "/" << opt.steps
                  << " dmin_mean=" << dminMean
                  << " Fz_mean=" << fzMean << " weight=" << weight
                  << " minJ=" << minJf << " kappa=" << kappa
                  << " wall=" << wallTotal << "s\n";
        if (!pass) std::cout << "  reason=" << gateReason << "\n";
        return pass ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 3;
    }
}
