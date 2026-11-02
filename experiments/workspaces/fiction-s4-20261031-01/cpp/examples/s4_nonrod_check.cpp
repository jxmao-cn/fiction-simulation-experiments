// =============================================================================
// examples/s4_nonrod_check.cpp — S4 非杆几何（torus）NH tet FEM 验证
// =============================================================================
// 预登记见 docs/s3-s6a-execution-design-20261031.md 第 3 节：
//   - 几何：IPC torus.msh（x/y∈[-0.5,0.5]，z∈[-0.095,0.095]）
//   - 材料：Neo-Hookean tet FEM，E=1e6 Pa，nu=0.3，rho=1000 kg/m^3
//   - 重力：(0,0,-9.81)；固定 y>=0.45 的节点，其余自由；无接触
//   - 先静态满载（dt=1, staticSolve, tol=1e-8, maxIter=100, DBC k=1e8）；
//     静态通过后从参考形状、零速度跑 2 s 动态（dt=0.005, 400 步,
//     tol=1e-4, maxIter=50），用于关键帧与短片
//   - 每步检查：有限、min J>0.3、DBC 误差<1e-6；失败即停，不调参重跑
//
// 输出：summary.txt / static_steps.csv / static_state.csv /
//       dynamic_steps.csv / final_state.csv / keyframe_*.png /
//       frames/*.png / video.mp4
// =============================================================================

#include "energies/Gravity.h"
#include "energies/NeoHookeanTetEnergy.h"
#include "io/image_writer.h"
#include "io/msh_loader.h"
#include "io/obj_writer.h"
#include "io/rasterizer.h"
#include "io/video_encoder.h"
#include "scene/mass_spring_builder.h"
#include "time_integrator.h"
#include "types.h"

#include <Eigen/SVD>
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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace fiction;

constexpr double kE = 1.0e6;
constexpr double kNu = 0.3;
constexpr double kRho = 1000.0;
constexpr double kGravity = 9.81;
constexpr double kFixedY = 0.45;
constexpr double kInitialDbcK = 1.0e8;

constexpr double kStaticTol = 1e-8;
constexpr int kStaticMaxIter = 100;
constexpr double kDynamicDt = 0.005;
constexpr int kDynamicSteps = 400;
constexpr double kDynamicTol = 1e-4;
constexpr int kDynamicMaxIter = 50;

constexpr int kMinJGateNumerator = 3;  // min J > 3/10
constexpr int kMinJGateDenominator = 10;
constexpr double kDbcErrorGate = 1e-6;

Eigen::Matrix3d edges(const Positions& x, const std::array<int, 4>& t) {
    Eigen::Matrix3d d;
    for (int c = 0; c < 3; ++c) d.col(c) = (x.row(t[c + 1]) - x.row(t[0])).transpose();
    return d;
}

struct Metrics {
    double minJ = 0.0, maxJ = 0.0;
    double minStretch = 0.0, maxStretch = 0.0;
    double dbcError = 0.0, maxDisp = 0.0, maxSpeed = 0.0;
    double elasticEnergy = 0.0, gravityEnergy = 0.0, potentialEnergy = 0.0;
    bool finite = true;
};

Metrics measure(const MshMesh& mesh, const Positions& x, const Velocities& v,
                const DBCInfo& dbc, const NeoHookeanTetEnergy& elastic,
                const ExampleGravityEnergy& gravity) {
    Metrics m;
    m.minJ = std::numeric_limits<double>::infinity();
    m.maxJ = -m.minJ;
    m.minStretch = m.minJ;
    m.maxStretch = -m.minJ;
    m.finite = x.allFinite() && v.allFinite();
    for (const auto& t : mesh.tets) {
        const Eigen::Matrix3d F = edges(x, t) * edges(mesh.x, t).inverse();
        const double J = F.determinant();
        m.finite = m.finite && F.allFinite() && std::isfinite(J);
        m.minJ = std::min(m.minJ, J);
        m.maxJ = std::max(m.maxJ, J);
        if (F.allFinite()) {
            const Eigen::Vector3d s = Eigen::JacobiSVD<Eigen::Matrix3d>(F).singularValues();
            m.minStretch = std::min(m.minStretch, s.minCoeff());
            m.maxStretch = std::max(m.maxStretch, s.maxCoeff());
        }
    }
    for (int i : dbc.indices) {
        m.dbcError = std::max(m.dbcError, (x.row(i) - mesh.x.row(i)).cwiseAbs().maxCoeff());
    }
    for (int i = 0; i < x.rows(); ++i) {
        m.maxDisp = std::max(m.maxDisp, (x.row(i) - mesh.x.row(i)).norm());
    }
    m.maxSpeed = v.cwiseAbs().maxCoeff();
    m.elasticEnergy = elastic.val(x);
    m.gravityEnergy = gravity.val(x);
    m.potentialEnergy = m.elasticEnergy + m.gravityEnergy;
    for (double value : {m.minJ, m.maxJ, m.minStretch, m.maxStretch, m.dbcError, m.maxDisp,
                         m.maxSpeed, m.elasticEnergy, m.gravityEnergy, m.potentialEnergy}) {
        m.finite = m.finite && std::isfinite(value);
    }
    return m;
}

bool staticPass(const TimeIntegrator::StepReport& r, const Metrics& m) {
    return m.finite && r.converged && !r.maxIterHit && !r.lineSearchFailed &&
           m.minJ > static_cast<double>(kMinJGateNumerator) / kMinJGateDenominator &&
           m.dbcError < kDbcErrorGate;
}

bool dynamicStepPass(const TimeIntegrator::StepReport& r, const Metrics& m) {
    return m.finite && r.converged && !r.maxIterHit && !r.lineSearchFailed &&
           m.minJ > static_cast<double>(kMinJGateNumerator) / kMinJGateDenominator &&
           m.dbcError < kDbcErrorGate;
}

const char* kStepHeader =
    "step,t,dt,static_solve,converged,iterations,line_search_iterations,max_iter_hit,"
    "line_search_failed,search_direction_rate,free_gradient_inf,min_J,max_J,min_stretch,"
    "max_stretch,dbc_max_error,measured_dbc_error,dbc_stiffness,dbc_satisfied,dbc_total,"
    "max_disp,max_speed,elastic_energy,gravity_energy,total_potential,wall_seconds,solver,"
    "gate_pass\n";

void writeStepRow(std::ofstream& out, int step, double t, double dt, bool staticSolve,
                  const TimeIntegrator::StepReport& r, const Metrics& m, double wall,
                  bool pass) {
    out << step << ',' << t << ',' << dt << ',' << int(staticSolve) << ',' << r.converged << ','
        << r.iterations << ',' << r.lineSearchIterations << ',' << r.maxIterHit << ','
        << r.lineSearchFailed << ',' << r.searchDirectionRate << ',' << r.freeGradientInf << ','
        << m.minJ << ',' << m.maxJ << ',' << m.minStretch << ',' << m.maxStretch << ','
        << r.dbcMaxError << ',' << m.dbcError << ',' << r.dbcStiffness << ',' << r.dbcSatisfied
        << ',' << r.dbcTotal << ',' << m.maxDisp << ',' << m.maxSpeed << ',' << m.elasticEnergy
        << ',' << m.gravityEnergy << ',' << m.potentialEnergy << ',' << wall << ',' << r.solver
        << ',' << pass << '\n';
}

void writeState(const std::string& path, const Positions& x0, const Positions& x,
                const std::vector<char>& isDbc) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("无法创建状态文件: " + path);
    out << std::setprecision(17);
    out << "node,X,Y,Z,x,y,z,ux,uy,uz,is_dbc\n";
    for (int i = 0; i < x.rows(); ++i) {
        out << i;
        for (int d = 0; d < 3; ++d) out << ',' << x0(i, d);
        for (int d = 0; d < 3; ++d) out << ',' << x(i, d);
        for (int d = 0; d < 3; ++d) out << ',' << x(i, d) - x0(i, d);
        out << ',' << int(isDbc[i]) << '\n';
    }
    out.close();
    if (!out) throw std::runtime_error("状态文件写入失败: " + path);
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::string meshPath, outDir;
        int steps = kDynamicSteps;
        int exportEvery = 5;
        double fps = 25.0;
        int videoW = 960, videoH = 720;
        bool render = true;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--help") {
                std::cout << "s4_nonrod_check --mesh FILE --out DIR [--steps N] [--export-every N]\n"
                             "  [--fps F] [--video-size W H] [--no-render]\n";
                return 0;
            }
            if (i + 1 >= argc) throw std::runtime_error("缺少参数值: " + a);
            if (a == "--mesh") meshPath = argv[++i];
            else if (a == "--out") outDir = argv[++i];
            else if (a == "--steps") steps = std::stoi(argv[++i]);
            else if (a == "--export-every") exportEvery = std::stoi(argv[++i]);
            else if (a == "--fps") fps = std::stod(argv[++i]);
            else if (a == "--no-render") { render = false; --i; }
            else if (a == "--video-size") {
                videoW = std::stoi(argv[++i]);
                if (i + 1 >= argc) throw std::runtime_error("--video-size 需要宽高");
                videoH = std::stoi(argv[++i]);
            } else throw std::runtime_error("未知参数: " + a);
        }
        if (meshPath.empty() || outDir.empty())
            throw std::runtime_error("必须指定 --mesh 与 --out");
        if (steps <= 0 || exportEvery <= 0) throw std::runtime_error("steps/export-every 必须为正");

        const MshMesh mesh = loadMsh(meshPath);
        if (mesh.empty() || mesh.tets.empty() || !mesh.x.allFinite())
            throw std::runtime_error("空网格或非有限坐标");
        const auto model = buildMassSpringModel(mesh, kRho);

        // 参考定向（与审计脚本一致）
        for (const auto& t : mesh.tets)
            if (!(edges(mesh.x, t).determinant() > 0.0))
                throw std::runtime_error("参考四面体非正定向");

        std::vector<int> fixedNodes;
        DBCInfo dbc;
        dbc.stiffness = kInitialDbcK;
        dbc.isDBC.assign(mesh.x.rows(), 0);
        for (int i = 0; i < mesh.x.rows(); ++i) {
            if (mesh.x(i, 1) >= kFixedY) {
                fixedNodes.push_back(i);
                dbc.indices.push_back(i);
                dbc.targetPositions.push_back(mesh.x.row(i).transpose());
                dbc.targetVelocities.push_back(Vec3::Zero());
                dbc.limits.push_back(Vec3::Zero());
                dbc.isDBC[i] = 1;
            }
        }
        if (fixedNodes.empty()) throw std::runtime_error("y>=0.45 未选中任何节点");
        std::array<double, 6> fixedBox{{1e300, 1e300, 1e300, -1e300, -1e300, -1e300}};
        for (int i : fixedNodes)
            for (int d = 0; d < 3; ++d) {
                fixedBox[d] = std::min(fixedBox[d], mesh.x(i, d));
                fixedBox[3 + d] = std::max(fixedBox[3 + d], mesh.x(i, d));
            }

        std::filesystem::create_directories(outDir);
        if (render) {
            std::error_code ec;
            std::filesystem::remove_all(std::filesystem::path(outDir) / "frames", ec);
            std::filesystem::create_directories(std::filesystem::path(outDir) / "frames");
        }

        auto elastic = std::make_shared<NeoHookeanTetEnergy>(mesh.tets, mesh.x, kE, kNu);
        auto gravity = std::make_shared<ExampleGravityEnergy>(model.masses, Vec3(0.0, 0.0, -kGravity));

        std::ofstream summary(std::filesystem::path(outDir) / "summary.txt");
        if (!summary) throw std::runtime_error("无法创建 summary.txt");
        summary << std::setprecision(17);
        summary << "run_id=s4-nonrod-20261031-01\n"
                << "mesh=" << meshPath << "\n"
                << "nodes=" << mesh.x.rows() << "\n"
                << "tets=" << mesh.tets.size() << "\n"
                << "surface_triangles=" << mesh.surface.size() << "\n"
                << "volume=" << model.volume << "\n"
                << "total_mass=" << model.totalMass << "\n"
                << "bbox_min=" << model.bbox[0] << ',' << model.bbox[1] << ',' << model.bbox[2] << "\n"
                << "bbox_max=" << model.bbox[3] << ',' << model.bbox[4] << ',' << model.bbox[5] << "\n"
                << "youngs=" << kE << "\npoisson=" << kNu << "\ndensity=" << kRho << "\n"
                << "gravity=0,0," << -kGravity << "\n"
                << "fixed_rule=y>=0.45\n"
                << "fixed_nodes=" << fixedNodes.size() << "\n"
                << "fixed_bbox_min=" << fixedBox[0] << ',' << fixedBox[1] << ',' << fixedBox[2] << "\n"
                << "fixed_bbox_max=" << fixedBox[3] << ',' << fixedBox[4] << ',' << fixedBox[5] << "\n"
                << "free_nodes=" << (mesh.x.rows() - static_cast<int>(fixedNodes.size())) << "\n"
                << "static_dt=1\nstatic_tol=" << kStaticTol << "\nstatic_max_iter=" << kStaticMaxIter
                << "\ndynamic_dt=" << kDynamicDt << "\ndynamic_steps=" << steps
                << "\ndynamic_tol=" << kDynamicTol << "\ndynamic_max_iter=" << kDynamicMaxIter
                << "\ninitial_dbc_stiffness=" << kInitialDbcK << "\n"
                << "dynamic_initial_condition=reference_shape_zero_velocity\n"
                << "gate_min_J_gt=0.3\ngate_dbc_error_lt=" << kDbcErrorGate << "\n";
        summary.close();
        if (!summary) throw std::runtime_error("summary.txt 写入失败");

        // ---------------- 静态满载 ----------------
        TimeIntegrator::Config scfg;
        scfg.h = 1.0;
        scfg.tol = kStaticTol;
        scfg.maxIter = kStaticMaxIter;
        scfg.staticSolve = true;
        scfg.verbose = false;
        TimeIntegrator staticSolver(scfg, model.masses);
        staticSolver.addEnergy(elastic);
        staticSolver.addEnergy(gravity);
        const auto t0 = std::chrono::steady_clock::now();
        const auto sstate = staticSolver.stepForward(mesh.x, Velocities::Zero(mesh.x.rows(), 3), dbc);
        const double staticWall =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const Metrics sm = measure(mesh, sstate.x, sstate.v, dbc, *elastic, *gravity);
        const bool sPass = staticPass(sstate.report, sm);

        std::ofstream staticCsv(std::filesystem::path(outDir) / "static_steps.csv");
        staticCsv << std::setprecision(17) << kStepHeader;
        writeStepRow(staticCsv, 0, 0.0, 1.0, true, sstate.report, sm, staticWall, sPass);
        staticCsv.close();
        if (!staticCsv) throw std::runtime_error("static_steps.csv 写入失败");
        writeState((std::filesystem::path(outDir) / "static_state.csv").string(), mesh.x, sstate.x,
                   dbc.isDBC);

        std::cout << std::setprecision(12) << "S4 static: converged=" << sstate.report.converged
                  << " iters=" << sstate.report.iterations << " minJ=" << sm.minJ
                  << " maxDisp=" << sm.maxDisp << " dbcErr=" << sm.dbcError
                  << " wall=" << staticWall << "s pass=" << sPass << '\n';

        if (!sPass) {
            std::ofstream fail(std::filesystem::path(outDir) / "FAILURE.txt");
            fail << "stage=static\nconverged=" << sstate.report.converged
                 << "\nmax_iter_hit=" << sstate.report.maxIterHit
                 << "\nline_search_failed=" << sstate.report.lineSearchFailed
                 << "\nfinite=" << sm.finite << "\nmin_J=" << sm.minJ
                 << "\ndbc_error=" << sm.dbcError << "\niterations=" << sstate.report.iterations
                 << "\nnote=静态未通过，未运行动态（不调参重跑）\n";
            fail.close();
            std::cerr << "S4 STATIC_FAIL\n";
            return 1;
        }

        // ---------------- 动态 2 s ----------------
        const int videoWEven = std::max(2, videoW - (videoW % 2));
        const int videoHEven = std::max(2, videoH - (videoH % 2));
        const Vec3 center = boxCenter(model.bbox);
        const double diag = std::max(boxDiagonal(model.bbox), 1e-6);
        const std::string ext = imageExtension();
        SoftRenderer renderer(render ? videoWEven : 2, render ? videoHEven : 2, 2);
        std::string framesDir = (std::filesystem::path(outDir) / "frames").string();
        int frameIndex = 0;

        auto renderFrame = [&](const Positions& x, const std::string& path) {
            SoftCamera cam;
            cam.target = center;
            cam.radius = 1.7 * diag;
            cam.yawDeg = 45.0;
            cam.pitchDeg = 25.0;
            renderer.clear(0.08f, 0.09f, 0.11f);
            renderer.setCamera(cam);
            const Vec3 groundOrigin(center.x(), center.y(), model.bbox[2]);
            renderer.drawGrid(groundOrigin, Vec3(0.0, 0.0, 1.0), 0.9 * diag, 0.9 * diag / 10.0,
                              0.35f, 0.38f, 0.42f);
            renderer.drawTriangles(x, mesh.surface, 0.45f, 0.62f, 0.90f);
            renderer.drawNodes(x, fixedNodes, 3, 1.0f, 0.25f, 0.2f);
            std::vector<unsigned char> rgb;
            renderer.downsample(rgb);
            if (!writeImage(path, videoWEven, videoHEven, rgb))
                std::cerr << "[s4] 帧写出失败: " << path << "\n";
        };
        auto frameName = [&](int idx) {
            std::ostringstream os;
            os << framesDir << "/frame_" << std::setw(6) << std::setfill('0') << idx << "." << ext;
            return os.str();
        };
        auto keyframe = [&](const Positions& x, int ms) {
            std::ostringstream os;
            os << outDir << "/keyframe_t" << std::setw(4) << std::setfill('0') << ms;
            const std::string base = os.str();
            writeSurfaceObj(base + ".obj", x, mesh.surface);
            if (render) renderFrame(x, base + "." + ext);
        };

        TimeIntegrator::Config dcfg;
        dcfg.h = kDynamicDt;
        dcfg.tol = kDynamicTol;
        dcfg.maxIter = kDynamicMaxIter;
        dcfg.staticSolve = false;
        dcfg.verbose = false;
        TimeIntegrator dynSolver(dcfg, model.masses);
        dynSolver.addEnergy(elastic);
        dynSolver.addEnergy(gravity);

        Positions x = mesh.x;
        Velocities v = Velocities::Zero(mesh.x.rows(), 3);
        std::ofstream dynCsv(std::filesystem::path(outDir) / "dynamic_steps.csv");
        if (!dynCsv) throw std::runtime_error("无法创建 dynamic_steps.csv");
        dynCsv << std::setprecision(17) << kStepHeader;

        keyframe(x, 0);
        if (render) renderFrame(x, frameName(frameIndex));
        ++frameIndex;

        const int midStep = steps / 2;
        const int midMs = static_cast<int>(std::lround(midStep * kDynamicDt * 1000.0));
        bool midKeyframeWritten = false;
        int failStep = -1;
        bool dynamicPass = true;
        double lastMaxDisp = 0.0, lastMinJ = sm.minJ;
        int lastIterations = 0;
        const auto d0 = std::chrono::steady_clock::now();
        for (int step = 1; step <= steps; ++step) {
            const auto b = std::chrono::steady_clock::now();
            const auto state = dynSolver.stepForward(x, v, dbc);
            const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - b).count();
            x = state.x;
            v = state.v;
            const Metrics m = measure(mesh, x, v, dbc, *elastic, *gravity);
            const bool pass = dynamicStepPass(state.report, m);
            writeStepRow(dynCsv, step, step * kDynamicDt, kDynamicDt, false, state.report, m, wall, pass);
            lastMaxDisp = m.maxDisp;
            lastMinJ = m.minJ;
            lastIterations = state.report.iterations;
            if (step == midStep) {
                keyframe(x, midMs);
                midKeyframeWritten = true;
            }
            if (render && step % exportEvery == 0) {
                renderFrame(x, frameName(frameIndex));
                ++frameIndex;
            }
            if (!pass) {
                failStep = step;
                dynamicPass = false;
                std::cerr << "S4 dynamic step " << step << " FAILED: converged="
                          << state.report.converged << " iters=" << state.report.iterations
                          << " minJ=" << m.minJ << " dbcErr=" << m.dbcError
                          << " finite=" << m.finite << '\n';
                break;
            }
            if (step % 50 == 0) {
                std::cout << "  step " << step << "/" << steps << " minJ=" << m.minJ
                          << " maxDisp=" << m.maxDisp << " iters=" << state.report.iterations
                          << " wall=" << wall << "s\n";
            }
        }
        dynCsv.close();
        if (!dynCsv) throw std::runtime_error("dynamic_steps.csv 写入失败");
        const double dynWall =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - d0).count();
        const int finalMs =
            static_cast<int>(std::lround((failStep < 0 ? steps : failStep) * kDynamicDt * 1000.0));
        keyframe(x, finalMs);
        writeState((std::filesystem::path(outDir) / "final_state.csv").string(), mesh.x, x, dbc.isDBC);

        // 视频编码（帧数 = steps/export-every + 1）
        std::string videoPath = (std::filesystem::path(outDir) / "video.mp4").string();
        std::ostringstream videoNote;
        if (render && frameIndex > 1) {
            const std::string pattern = ffmpegFramePattern(framesDir, ext);
            if (encodeVideo(pattern, fps, videoPath)) {
                videoNote << "video=" << videoPath << "\nvideo_frames=" << frameIndex
                          << "\nvideo_fps=" << fps
                          << "\nvideo_container_seconds=" << frameIndex / fps
                          << "\nvideo_physical_seconds=" << (frameIndex - 1) * exportEvery * kDynamicDt
                          << "\n";
            } else {
                videoNote << "video=ENCODE_FAILED\nframes_dir=" << framesDir << "\n";
            }
        } else {
            videoNote << "video=SKIPPED\n";
        }

        std::ostringstream midName;
        if (midKeyframeWritten) {
            midName << outDir << "/keyframe_t" << std::setw(4) << std::setfill('0') << midMs
                    << "." << ext;
        } else {
            midName << "NOT_WRITTEN";
        }
        std::ofstream sum2(std::filesystem::path(outDir) / "summary.txt", std::ios::app);
        sum2 << std::setprecision(17)
             << "static_converged=" << sstate.report.converged << "\n"
             << "static_iterations=" << sstate.report.iterations << "\n"
             << "static_min_J=" << sm.minJ << "\nstatic_max_J=" << sm.maxJ << "\n"
             << "static_min_stretch=" << sm.minStretch << "\nstatic_max_stretch=" << sm.maxStretch
             << "\nstatic_max_disp=" << sm.maxDisp << "\nstatic_dbc_error=" << sm.dbcError
             << "\nstatic_wall_seconds=" << staticWall << "\n"
             << "static_gate_pass=" << sPass << "\n"
             << "dynamic_steps_completed=" << (failStep < 0 ? steps : failStep - 1) << "\n"
             << "dynamic_fail_step=" << failStep << "\n"
             << "dynamic_gate_pass=" << dynamicPass << "\n"
             << "dynamic_final_min_J=" << lastMinJ << "\n"
             << "dynamic_final_max_disp=" << lastMaxDisp << "\n"
             << "dynamic_final_iterations=" << lastIterations << "\n"
             << "dynamic_wall_seconds=" << dynWall << "\n"
             << videoNote.str()
             << "keyframe_initial=" << outDir << "/keyframe_t0000." << ext << "\n"
             << "keyframe_middle=" << midName.str() << "\n"
             << "keyframe_final=" << outDir << "/keyframe_t"
             << std::setw(4) << std::setfill('0') << finalMs << "." << ext << "\n";
        sum2.close();

        std::cout << std::setprecision(12) << "S4 dynamic: completed="
                  << (failStep < 0 ? steps : failStep - 1) << "/" << steps
                  << " failStep=" << failStep << " minJ=" << lastMinJ
                  << " maxDisp=" << lastMaxDisp << " wall=" << dynWall << "s\n";
        if (!dynamicPass) {
            std::ofstream fail(std::filesystem::path(outDir) / "FAILURE.txt");
            fail << "stage=dynamic\nfail_step=" << failStep << "\nmin_J=" << lastMinJ << "\n"
                 << "note=动态某步未通过门槛，已停止（不调参重跑）\n";
            fail.close();
            return 1;
        }
        std::cout << "S4 PASS (static + dynamic gates, no analytical reference)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "S4 ERROR: " << e.what() << '\n';
        return 2;
    }
}
