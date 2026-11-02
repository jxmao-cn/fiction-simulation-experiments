// =============================================================================
// examples/s5_scale_check.cpp — S5 单体规模 pilot（Armadillo13K / 28K）
// =============================================================================
// 职责（对应 docs/s3-s6a-execution-design-20261031.md 第 4 节预登记规格）：
//   1. 只读资产审计：节点/四面体/表面、bbox、体积、质量、正定向、表面闭合
//   2. NH tet FEM + 重力 + pin-box DBC，隐式欧拉 + projected Newton 推进
//   3. 每步记录：wall、iterations、line_search_iterations、converged、maxIterHit、
//      lineSearchFailed、minJ、maxJ、max_disp、max_speed、solver
//   4. 参考装配与首次分解耗时、H 非零数（驱动侧复算，与 TimeIntegrator 同项）
//
// 用法：
//   s5_scale_check --mesh FILE --steps N --pin-box 6f --csv FILE --state FILE
//                  [--audit FILE] [--meta FILE] [--trace FILE.step1.csv]
//                  [--youngs 1e7] [--poisson 0.3] [--density 1000]
//                  [--dt 0.005] [--tol 1e-2] [--max-iter 50] [--gravity 0 0 -9.81]
//
// 退出码：0 全部步收敛且有限；3 出现不收敛/线搜索失败/非有限；2 异常。
// =============================================================================

#include "energies/Gravity.h"
#include "energies/NeoHookeanTetEnergy.h"
#include "io/msh_loader.h"
#include "scene/mass_spring_builder.h"
#include "time_integrator.h"
#include "types.h"

// 与 src/time_integrator.cpp 相同的 CHOLMOD 探测（用于单独计首次分解耗时）
#if defined(FICTION_USE_CHOLMOD) && __has_include(<Eigen/CholmodSupport>)
#  if __has_include(<cholmod.h>) || __has_include(<suitesparse/cholmod.h>)
#    include <Eigen/CholmodSupport>
#    define S5_HAS_CHOLMOD 1
#  endif
#endif

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

using namespace fiction;

double nowSeconds() {
    return std::chrono::duration<double>(
               std::chrono::steady_clock::now().time_since_epoch()).count();
}

/// 单四面体参考行列式（Dm = [x1-x0, x2-x0, x3-x0] 的列）
double refDet(const Positions& x0, const std::array<int, 4>& tet) {
    const Vec3 p0 = x0.row(tet[0]).transpose();
    Eigen::Matrix3d Dm;
    Dm.col(0) = x0.row(tet[1]).transpose() - p0;
    Dm.col(1) = x0.row(tet[2]).transpose() - p0;
    Dm.col(2) = x0.row(tet[3]).transpose() - p0;
    return Dm.determinant();
}

std::uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint32_t>(b);
}

std::string boxJson(const std::array<double, 6>& b) {
    std::ostringstream os;
    os << std::setprecision(17) << "[" << b[0] << ", " << b[1] << ", " << b[2] << ", "
       << b[3] << ", " << b[4] << ", " << b[5] << "]";
    return os.str();
}

struct Audit {
    int nodes = 0;
    int tets = 0;
    int surface = 0;
    std::array<double, 6> bbox{};
    double volume = 0.0;
    double mass = 0.0;
    double detMin = 0.0;
    double detMax = 0.0;
    int nonPositive = 0;
    int degenerate = 0;
    int orphanNodes = 0;
    bool surfaceClosed = false;
    int nonManifoldEdges = 0;
    int boundaryEdgesOnce = 0;
    int boundaryFaceMismatch = 0;
};

Audit auditMesh(const MshMesh& mesh, double rho) {
    Audit a;
    a.nodes = static_cast<int>(mesh.x.rows());
    a.tets = mesh.numTets();
    a.surface = mesh.numSurface();
    a.bbox = boundingBox(mesh.x);
    a.detMin = std::numeric_limits<double>::infinity();
    a.detMax = -a.detMin;
    std::vector<char> used(static_cast<size_t>(a.nodes), 0);
    for (const auto& t : mesh.tets) {
        const double d = refDet(mesh.x, t);
        a.detMin = std::min(a.detMin, d);
        a.detMax = std::max(a.detMax, d);
        if (!(d > 0.0)) ++a.nonPositive;
        if (std::abs(d) / 6.0 <= 1e-14) ++a.degenerate;
        a.volume += tetVolume(mesh.x, t);
        for (int i : t) if (i >= 0 && i < a.nodes) used[static_cast<size_t>(i)] = 1;
    }
    for (char u : used) if (!u) ++a.orphanNodes;
    a.mass = rho * a.volume;

    // 表面闭合：每条无向边恰好被 2 个表面三角形共享
    std::unordered_map<std::uint64_t, int> edgeUse;
    edgeUse.reserve(static_cast<size_t>(a.surface) * 3);
    std::unordered_set<std::uint64_t> surfaceFaces;
    surfaceFaces.reserve(static_cast<size_t>(a.surface) * 2);
    for (const auto& f : mesh.surface) {
        std::array<int, 3> s{{f[0], f[1], f[2]}};
        std::sort(s.begin(), s.end());
        surfaceFaces.insert((static_cast<std::uint64_t>(s[0]) << 40) ^
                            (static_cast<std::uint64_t>(s[1]) << 20) ^
                            static_cast<std::uint64_t>(s[2]));
        for (int k = 0; k < 3; ++k) ++edgeUse[edgeKey(f[k], f[(k + 1) % 3])];
    }
    for (const auto& kv : edgeUse) {
        if (kv.second == 1) ++a.boundaryEdgesOnce;
        else if (kv.second != 2) ++a.nonManifoldEdges;
    }
    a.surfaceClosed = (a.boundaryEdgesOnce == 0 && a.nonManifoldEdges == 0);

    // 表面三角面集合应等于四面体边界（出现一次的面）集合
    std::unordered_map<std::uint64_t, int> faceUse;
    faceUse.reserve(static_cast<size_t>(a.tets) * 2);
    for (const auto& t : mesh.tets) {
        static const int kF[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
        for (const auto& f : kF) {
            std::array<int, 3> s{{t[f[0]], t[f[1]], t[f[2]]}};
            std::sort(s.begin(), s.end());
            ++faceUse[(static_cast<std::uint64_t>(s[0]) << 40) ^
                      (static_cast<std::uint64_t>(s[1]) << 20) ^
                      static_cast<std::uint64_t>(s[2])];
        }
    }
    for (const auto& kv : faceUse)
        if (kv.second == 1 && !surfaceFaces.count(kv.first)) ++a.boundaryFaceMismatch;
    for (std::uint64_t key : surfaceFaces)
        if (faceUse.find(key) == faceUse.end() || faceUse[key] != 1) ++a.boundaryFaceMismatch;
    return a;
}

void writeAuditJson(const std::string& path, const Audit& a,
                    const std::string& meshPath, double rho) {
    std::ofstream os(path);
    if (!os) throw std::runtime_error("无法创建审计 JSON: " + path);
    os << std::setprecision(17);
    os << "{\n  \"mesh\": \"" << meshPath << "\",\n"
       << "  \"nodes\": " << a.nodes << ",\n"
       << "  \"tets\": " << a.tets << ",\n"
       << "  \"surface_triangles\": " << a.surface << ",\n"
       << "  \"bbox\": " << boxJson(a.bbox) << ",\n"
       << "  \"volume\": " << a.volume << ",\n"
       << "  \"density\": " << rho << ",\n"
       << "  \"total_mass\": " << a.mass << ",\n"
       << "  \"reference_det_min\": " << a.detMin << ",\n"
       << "  \"reference_det_max\": " << a.detMax << ",\n"
       << "  \"non_positive_tets\": " << a.nonPositive << ",\n"
       << "  \"degenerate_tets\": " << a.degenerate << ",\n"
       << "  \"orphan_nodes\": " << a.orphanNodes << ",\n"
       << "  \"surface_closed\": " << (a.surfaceClosed ? "true" : "false") << ",\n"
       << "  \"surface_boundary_edges_once\": " << a.boundaryEdgesOnce << ",\n"
       << "  \"surface_nonmanifold_edges\": " << a.nonManifoldEdges << ",\n"
       << "  \"surface_tet_face_mismatch\": " << a.boundaryFaceMismatch << "\n}\n";
    if (!os) throw std::runtime_error("写审计 JSON 失败: " + path);
}

} // namespace

int main(int argc, char** argv) {
    try {
        // ------------------------- CLI -------------------------
        std::string meshPath, csvPath, statePath, auditPath, metaPath, tracePath, label = "S5";
        int steps = 1, maxIter = 50;
        double youngs = 1e7, poisson = 0.3, rho = 1000.0, dt = 0.005, tol = 1e-2;
        Vec3 gravity(0.0, 0.0, -9.81);
        std::array<double, 6> pinBox{{1e300, 1e300, 1e300, -1e300, -1e300, -1e300}};

        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            auto need = [&](const char* name) {
                if (i + 1 >= argc) throw std::runtime_error(std::string("缺少参数值: ") + name);
            };
            if (a == "--help") {
                std::cout << "s5_scale_check --mesh F --steps N --pin-box 6f --csv F --state F "
                             "[--audit F] [--meta F] [--trace F] [--youngs E] [--poisson nu] "
                             "[--density rho] [--dt h] [--tol t] [--max-iter n] "
                             "[--gravity gx gy gz] [--label L]\n";
                return 0;
            } else if (a == "--mesh") { need("--mesh"); meshPath = argv[++i]; }
            else if (a == "--csv") { need("--csv"); csvPath = argv[++i]; }
            else if (a == "--state") { need("--state"); statePath = argv[++i]; }
            else if (a == "--audit") { need("--audit"); auditPath = argv[++i]; }
            else if (a == "--meta") { need("--meta"); metaPath = argv[++i]; }
            else if (a == "--trace") { need("--trace"); tracePath = argv[++i]; }
            else if (a == "--label") { need("--label"); label = argv[++i]; }
            else if (a == "--steps") { need("--steps"); steps = std::stoi(argv[++i]); }
            else if (a == "--max-iter") { need("--max-iter"); maxIter = std::stoi(argv[++i]); }
            else if (a == "--youngs") { need("--youngs"); youngs = std::stod(argv[++i]); }
            else if (a == "--poisson") { need("--poisson"); poisson = std::stod(argv[++i]); }
            else if (a == "--density") { need("--density"); rho = std::stod(argv[++i]); }
            else if (a == "--dt") { need("--dt"); dt = std::stod(argv[++i]); }
            else if (a == "--tol") { need("--tol"); tol = std::stod(argv[++i]); }
            else if (a == "--gravity") {
                need("--gravity");
                if (i + 3 >= argc) throw std::runtime_error("--gravity 需要 3 个参数");
                gravity = Vec3(std::stod(argv[i + 1]), std::stod(argv[i + 2]), std::stod(argv[i + 3]));
                i += 3;
            } else if (a == "--pin-box") {
                need("--pin-box");
                if (i + 6 >= argc) throw std::runtime_error("--pin-box 需要 6 个参数");
                for (int d = 0; d < 6; ++d) pinBox[d] = std::stod(argv[++i]);
            } else {
                throw std::runtime_error("未知参数: " + a);
            }
        }
        if (meshPath.empty() || csvPath.empty() || statePath.empty())
            throw std::runtime_error("必须指定 --mesh、--csv、--state");
        if (steps < 1) throw std::runtime_error("--steps 必须 >= 1");
        for (const auto& p : {csvPath, statePath, auditPath, metaPath, tracePath}) {
            if (p.empty()) continue;
            const auto parent = std::filesystem::path(p).parent_path();
            if (!parent.empty()) std::filesystem::create_directories(parent);
        }

        // ------------------------- 资产审计（只读） -------------------------
        const double tLoad0 = nowSeconds();
        const MshMesh mesh = loadMsh(meshPath);
        const double tLoad = nowSeconds() - tLoad0;
        if (mesh.empty() || mesh.tets.empty() || !mesh.x.allFinite())
            throw std::runtime_error("空网格或非有限坐标");
        const Audit audit = auditMesh(mesh, rho);
        if (!auditPath.empty()) writeAuditJson(auditPath, audit, meshPath, rho);
        std::cout << std::setprecision(12)
                  << "S5 审计: nodes=" << audit.nodes << " tets=" << audit.tets
                  << " surface=" << audit.surface << " bbox=" << boxJson(audit.bbox)
                  << " volume=" << audit.volume << " mass=" << audit.mass
                  << " det_min=" << audit.detMin << " nonpos=" << audit.nonPositive
                  << " closed=" << (audit.surfaceClosed ? 1 : 0)
                  << " face_mismatch=" << audit.boundaryFaceMismatch << '\n';

        // ------------------------- 模型与边界 -------------------------
        const double tModel0 = nowSeconds();
        const MassSpringModel model = buildMassSpringModel(mesh, rho);
        const double tModel = nowSeconds() - tModel0;

        const std::vector<int> pinned = nodesInBox(model.x0, pinBox);
        if (pinned.empty()) throw std::runtime_error("pin-box 未选中任何节点");
        const int N = model.numNodes();
        DBCInfo dbc;
        dbc.stiffness = 1e8;
        dbc.isDBC.assign(static_cast<size_t>(N), 0);
        for (int i : pinned) {
            dbc.indices.push_back(i);
            dbc.targetVelocities.push_back(Vec3::Zero());
            dbc.limits.push_back(model.x0.row(i).transpose());
            dbc.isDBC[static_cast<size_t>(i)] = 1;
        }
        const int dof = N * kDim;
        const int freeDof = (N - static_cast<int>(pinned.size())) * kDim;

        // ------------------------- 能量组装 -------------------------
        const double tEnergy0 = nowSeconds();
        auto elastic = std::make_shared<NeoHookeanTetEnergy>(mesh.tets, model.x0, youngs, poisson);
        auto grav = std::make_shared<ExampleGravityEnergy>(model.masses, gravity);
        const double tEnergy = nowSeconds() - tEnergy0;

        // ------------------------- 参考装配 + 首次分解（驱动侧同项复算） ---------
        const double tAsm0 = nowSeconds();
        std::vector<Triplet> trips;
        trips.reserve(static_cast<size_t>(dof) * 4 + static_cast<size_t>(mesh.tets.size()) * 144);
        for (int i = 0; i < N; ++i)
            for (int d = 0; d < kDim; ++d) {
                const int idx = kDim * i + d;
                trips.emplace_back(idx, idx, model.masses[i]);
            }
        const size_t before = trips.size();
        elastic->addHess(model.x0, trips);
        const double wElastic = elastic->hessianWeight(dt);
        for (size_t i = before; i < trips.size(); ++i)
            trips[i] = Triplet(trips[i].row(), trips[i].col(), trips[i].value() * wElastic);
        for (int i : pinned)
            for (int d = 0; d < kDim; ++d) {
                const int idx = kDim * i + d;
                trips.emplace_back(idx, idx, dbc.stiffness * model.masses[i]);
            }
        SparseMatrix H(dof, dof);
        H.setFromTriplets(trips.begin(), trips.end());
        trips.clear();
        trips.shrink_to_fit();
        const double tAsm = nowSeconds() - tAsm0;
        const Eigen::Index nnz = H.nonZeros();

        std::string firstFactorSolver = "none";
        double tFactor = -1.0;
        {
            const auto f0 = nowSeconds();
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dof);
#ifdef S5_HAS_CHOLMOD
            Eigen::CholmodSupernodalLLT<SparseMatrix> chol;
            chol.compute(H);
            if (chol.info() == Eigen::Success && chol.solve(rhs).allFinite()) {
                firstFactorSolver = "CHOLMOD";
                tFactor = nowSeconds() - f0;
            }
#endif
            if (tFactor < 0.0) {
                Eigen::SimplicialLDLT<SparseMatrix> ldlt;
                ldlt.compute(H);
                if (ldlt.info() == Eigen::Success) {
                    firstFactorSolver = "SimplicialLDLT";
                    tFactor = nowSeconds() - f0;
                } else {
                    firstFactorSolver = "failed";
                }
            }
        }

        // ------------------------- 逐步推进 -------------------------
        std::vector<double> detRef(static_cast<size_t>(mesh.tets.size()), 0.0);
        for (size_t t = 0; t < mesh.tets.size(); ++t) detRef[t] = refDet(model.x0, mesh.tets[t]);

        std::ofstream csv(csvPath);
        if (!csv) throw std::runtime_error("无法创建 CSV: " + csvPath);
        csv << std::setprecision(17)
            << "step,wall_seconds,iterations,line_search_iterations,converged,maxIterHit,"
               "lineSearchFailed,minJ,maxJ,max_disp,max_speed,solver,dbc_max_error,free_gradient_inf\n";

        Positions x = model.x0;
        Velocities v = Velocities::Zero(N, kDim);
        bool allOk = true;
        bool allConverged = true;
        double totalWall = 0.0;
        int stepsDone = 0;

        for (int step = 1; step <= steps; ++step) {
            TimeIntegrator::Config cfg;
            cfg.h = dt;
            cfg.tol = tol;
            cfg.maxIter = maxIter;
            cfg.verbose = false;
            if (step == 1 && !tracePath.empty()) {
                cfg.traceCsvPath = tracePath;
                cfg.traceLabel = label;
            }
            TimeIntegrator solver(cfg, model.masses);
            solver.addEnergy(grav);
            solver.addEnergy(elastic);
            const auto s0 = std::chrono::steady_clock::now();
            const auto state = solver.stepForward(x, v, dbc);
            const double wall = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - s0).count();
            x = state.x;
            v = state.v;
            const auto& r = state.report;
            totalWall += wall;

            // 每步诊断：J = det(Ds) / det(Dm)，位移与速度取逐节点范数
            double minJ = std::numeric_limits<double>::infinity();
            double maxJ = -minJ;
            bool finite = x.allFinite() && v.allFinite();
            for (size_t t = 0; t < mesh.tets.size() && finite; ++t) {
                const auto& tet = mesh.tets[t];
                const Vec3 p0 = x.row(tet[0]).transpose();
                Eigen::Matrix3d Ds;
                Ds.col(0) = x.row(tet[1]).transpose() - p0;
                Ds.col(1) = x.row(tet[2]).transpose() - p0;
                Ds.col(2) = x.row(tet[3]).transpose() - p0;
                const double J = Ds.determinant() / detRef[t];
                if (!std::isfinite(J)) { finite = false; break; }
                minJ = std::min(minJ, J);
                maxJ = std::max(maxJ, J);
            }
            double maxDisp = 0.0, maxSpeed = 0.0;
            for (int i = 0; i < N; ++i) {
                maxDisp = std::max(maxDisp, (x.row(i) - model.x0.row(i)).norm());
                maxSpeed = std::max(maxSpeed, v.row(i).norm());
            }
            allOk = allOk && finite;
            allConverged = allConverged && r.converged;
            ++stepsDone;

            csv << step << ',' << wall << ',' << r.iterations << ',' << r.lineSearchIterations
                << ',' << (r.converged ? 1 : 0) << ',' << (r.maxIterHit ? 1 : 0) << ','
                << (r.lineSearchFailed ? 1 : 0) << ',' << minJ << ',' << maxJ << ','
                << maxDisp << ',' << maxSpeed << ',' << r.solver << ','
                << r.dbcMaxError << ',' << r.freeGradientInf << '\n';
            csv.flush();

            std::cout << std::setprecision(6) << "step " << step << ": wall=" << wall
                      << " s iter=" << r.iterations << " ls=" << r.lineSearchIterations
                      << " conv=" << (r.converged ? 1 : 0) << " maxIterHit=" << (r.maxIterHit ? 1 : 0)
                      << " lsFail=" << (r.lineSearchFailed ? 1 : 0) << " minJ=" << minJ
                      << " maxJ=" << maxJ << " maxDisp=" << maxDisp
                      << " maxSpeed=" << maxSpeed << " solver=" << r.solver << '\n';

            if (!finite || r.maxIterHit || r.lineSearchFailed || !r.converged) {
                std::cout << "[STOP] 第 " << step << " 步未通过（有限=" << (finite ? 1 : 0)
                          << " conv=" << (r.converged ? 1 : 0) << " maxIterHit=" << r.maxIterHit
                          << " lsFail=" << r.lineSearchFailed << "）\n";
                break;
            }
        }
        csv.close();

        // ------------------------- 末态与元数据 -------------------------
        {
            std::ofstream st(statePath);
            if (!st) throw std::runtime_error("无法创建状态文件: " + statePath);
            st << std::setprecision(17) << "node,X,Y,Z,x,y,z,ux,uy,uz,is_dbc\n";
            for (int i = 0; i < N; ++i) {
                st << i;
                for (int d = 0; d < kDim; ++d) st << ',' << model.x0(i, d);
                for (int d = 0; d < kDim; ++d) st << ',' << x(i, d);
                for (int d = 0; d < kDim; ++d) st << ',' << x(i, d) - model.x0(i, d);
                st << ',' << (dbc.isDBC[static_cast<size_t>(i)] ? 1 : 0) << '\n';
            }
            if (!st) throw std::runtime_error("写状态文件失败: " + statePath);
        }
        if (!metaPath.empty()) {
            std::ofstream meta(metaPath);
            if (!meta) throw std::runtime_error("无法创建元数据 JSON: " + metaPath);
            meta << std::setprecision(17)
                 << "{\n  \"label\": \"" << label << "\",\n"
                 << "  \"mesh\": \"" << meshPath << "\",\n"
                 << "  \"nodes\": " << N << ",\n"
                 << "  \"tets\": " << mesh.numTets() << ",\n"
                 << "  \"surface_triangles\": " << mesh.numSurface() << ",\n"
                 << "  \"dof\": " << dof << ",\n"
                 << "  \"free_dof\": " << freeDof << ",\n"
                 << "  \"pin_box\": " << boxJson(pinBox) << ",\n"
                 << "  \"pin_nodes\": " << pinned.size() << ",\n"
                 << "  \"volume\": " << model.volume << ",\n"
                 << "  \"total_mass\": " << model.totalMass << ",\n"
                 << "  \"avg_edge_length\": " << model.avgEdgeLength << ",\n"
                 << "  \"youngs\": " << youngs << ",\n"
                 << "  \"poisson\": " << poisson << ",\n"
                 << "  \"density\": " << rho << ",\n"
                 << "  \"gravity\": [" << gravity.x() << ", " << gravity.y() << ", "
                 << gravity.z() << "],\n"
                 << "  \"dt\": " << dt << ",\n"
                 << "  \"tol\": " << tol << ",\n"
                 << "  \"max_iter\": " << maxIter << ",\n"
                 << "  \"dbc_stiffness_initial\": " << dbc.stiffness << ",\n"
                 << "  \"hessian_nnz\": " << nnz << ",\n"
                 << "  \"load_seconds\": " << tLoad << ",\n"
                 << "  \"model_build_seconds\": " << tModel << ",\n"
                 << "  \"energy_init_seconds\": " << tEnergy << ",\n"
                 << "  \"reference_assembly_seconds\": " << tAsm << ",\n"
                 << "  \"first_factorization_seconds\": " << tFactor << ",\n"
                 << "  \"first_factorization_solver\": \"" << firstFactorSolver << "\",\n"
                 << "  \"steps_requested\": " << steps << ",\n"
                 << "  \"steps_completed\": " << stepsDone << ",\n"
                 << "  \"total_step_wall_seconds\": " << totalWall << ",\n"
                 << "  \"all_steps_converged\": " << (allConverged ? "true" : "false") << ",\n"
                 << "  \"all_finite\": " << (allOk ? "true" : "false") << "\n}\n";
            if (!meta) throw std::runtime_error("写元数据 JSON 失败: " + metaPath);
        }

        const bool pass = allOk && allConverged && stepsDone == steps;
        std::cout << std::setprecision(6) << "S5 汇总: steps=" << stepsDone << '/' << steps
                  << " total_step_wall=" << totalWall << " s nnz_H=" << nnz
                  << " first_factor=" << tFactor << " s (" << firstFactorSolver << ")"
                  << " pin_nodes=" << pinned.size() << " free_dof=" << freeDof << '\n'
                  << (pass ? "[PASS] S5 pilot 短段完成" : "[FAIL] S5 pilot 未完成或未收敛") << '\n';
        return pass ? 0 : 3;
    } catch (const std::exception& e) {
        std::cerr << "S5 ERROR: " << e.what() << '\n';
        return 2;
    }
}
