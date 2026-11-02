// =============================================================================
// examples/ipc_mesh_demo.cpp — IPC .msh 模型控制台仿真（无可视化）
// =============================================================================
// 职责：
//   - 加载 IPC input/tetMeshes/*.msh（MshIO 4.1 ASCII 四面体网格）
//   - 装配 mass-spring（tet 唯一边）+ 集中质量 + pin-box DBC + force-box 外力
//   - 隐式欧拉 + projected Newton 推进，可选导出 OBJ 帧
//   - 可选导出 mp4 视频：CPU 软光栅化（无 OpenGL 依赖）+ ffmpeg，SSH 远程可用
//   - 结束时自检结果有限（无 NaN/Inf）；退出码 0/1 可供 ctest 捕获
//
// 用法：./ipc_mesh_demo --mesh <path.msh> [选项...]（详见 --help）
// =============================================================================

#include "io/image_writer.h"
#include "io/obj_writer.h"
#include "io/rasterizer.h"
#include "io/video_encoder.h"
#include "scene/ipc_scene.h"
#include "time_integrator.h"
#include "types.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

/// 帧文件名 frame_000000.png / .ppm
std::string framePath(const std::string& dir, int index, const std::string& ext) {
    std::ostringstream os;
    os << dir << "/frame_" << std::setw(6) << std::setfill('0') << index << "." << ext;
    return os.str();
}

} // namespace

int main(int argc, char** argv) {
    const std::vector<std::string> args(argv + 1, argv + argc);
    if (fiction::ipcOptionsWantHelp(args)) {
        fiction::printIpcSceneUsage(argv[0], std::cout);
        return 0;
    }

    fiction::IpcSceneOptions opt;
    std::string err;
    if (fiction::parseIpcSceneOptions(args, opt, err) != fiction::IpcParseResult::Ok) {
        std::cerr << "[ipc_mesh_demo] 参数错误: " << err << "\n\n";
        fiction::printIpcSceneUsage(argv[0], std::cerr);
        return 1;
    }

    try {
        fiction::IpcScene scene = fiction::buildIpcScene(opt);
        std::cout << "===== fiction IPC 模型仿真（mass-spring）=====\n"
                  << fiction::ipcSceneSummary(scene, opt) << "\n";

        fiction::TimeIntegrator::Config cfg;
        cfg.h = opt.dt;
        cfg.tol = opt.tol;
        cfg.maxIter = opt.maxIter;
        cfg.verbose = opt.verbose;
        fiction::TimeIntegrator integrator(cfg, scene.model.masses);
        for (auto& e : fiction::makeIpcSceneEnergies(scene, opt)) {
            integrator.addEnergy(std::move(e));
        }

        // 导出目录准备（OBJ 帧 0 = 初始状态）
        const bool exporting = !opt.exportDir.empty();
        if (exporting) {
            std::error_code ec;
            std::filesystem::create_directories(opt.exportDir, ec);
            if (ec) {
                std::cerr << "[ipc_mesh_demo] 无法创建导出目录 " << opt.exportDir
                          << ": " << ec.message() << "\n";
                return 1;
            }
            fiction::writeSurfaceObj(opt.exportDir + "/0.obj", scene.model.x0, scene.mesh.surface);
        }

        // 视频模式准备（CPU 软光栅化，无显示环境可用）
        const bool rendering = !opt.videoPath.empty();
        int videoW = 0, videoH = 0;
        std::string framesDir;
        fiction::SoftRenderer renderer(1, 1); // 占位，rendering 时重建
        if (rendering) {
            videoW = opt.videoWidth - (opt.videoWidth % 2);
            videoH = opt.videoHeight - (opt.videoHeight % 2);
            videoW = std::max(videoW, 2);
            videoH = std::max(videoH, 2);
            framesDir = opt.videoPath + ".frames";
            std::error_code ec;
            std::filesystem::create_directories(framesDir, ec);
            if (ec) {
                std::cerr << "[ipc_mesh_demo] 无法创建视频帧目录 " << framesDir
                          << ": " << ec.message() << "\n";
                return 1;
            }
            const std::filesystem::path videoFs(opt.videoPath);
            if (!videoFs.parent_path().empty()) {
                std::filesystem::create_directories(videoFs.parent_path(), ec);
                if (ec) {
                    std::cerr << "[ipc_mesh_demo] 无法创建视频输出目录 "
                              << videoFs.parent_path() << ": " << ec.message() << "\n";
                    return 1;
                }
            }
            renderer = fiction::SoftRenderer(videoW, videoH, 2);
        }

        // 相机与地面（按初始包围盒取景；环绕速度 18 deg/s 仿真时间）
        const fiction::Vec3 center = fiction::boxCenter(scene.model.bbox);
        const double diag = std::max(fiction::boxDiagonal(scene.model.bbox), 1e-6);
        const fiction::Vec3 groundOrigin(center.x(), center.y(), scene.model.bbox[2]);
        const double gridExtent = 0.9 * diag;
        const std::string ext = fiction::imageExtension();
        int frameIndex = 0;

        auto renderFrame = [&](const fiction::Positions& x, double simTime) {
            if (!rendering) return;
            fiction::SoftCamera cam;
            cam.target = center;
            cam.radius = 1.7 * diag;
            cam.yawDeg = 45.0 + (opt.videoOrbit ? 18.0 * simTime : 0.0);
            cam.pitchDeg = 25.0;
            renderer.clear(0.08f, 0.09f, 0.11f);
            renderer.setCamera(cam);
            renderer.drawGrid(groundOrigin, fiction::Vec3(0.0, 0.0, 1.0),
                              gridExtent, gridExtent / 10.0, 0.35f, 0.38f, 0.42f);
            renderer.drawTriangles(x, scene.mesh.surface, 0.45f, 0.62f, 0.90f);
            if (!scene.pinnedNodes.empty())
                renderer.drawNodes(x, scene.pinnedNodes, 3, 1.0f, 0.25f, 0.2f);
            for (const auto& nodes : scene.forceNodes)
                renderer.drawNodes(x, nodes, 3, 1.0f, 0.75f, 0.15f);
            std::vector<unsigned char> rgb;
            renderer.downsample(rgb);
            if (!fiction::writeImage(framePath(framesDir, frameIndex, ext), videoW, videoH, rgb)) {
                std::cerr << "[ipc_mesh_demo] 帧写出失败: " << framePath(framesDir, frameIndex, ext)
                          << "\n";
            }
            ++frameIndex;
        };

        fiction::Positions x = scene.model.x0;
        fiction::Velocities v = fiction::Velocities::Zero(scene.model.numNodes(), fiction::kDim);
        if (opt.initVel.squaredNorm() > 0.0) v.rowwise() = opt.initVel.transpose();

        renderFrame(x, 0.0); // 初始帧

        const auto t0 = std::chrono::steady_clock::now();
        for (int step = 1; step <= opt.steps; ++step) {
            fiction::refreshIpcDbcTargets(scene, opt, step * opt.dt); // 扭转 DBC 轨迹
            const auto state = integrator.stepForward(x, v, scene.dbc);
            x = state.x;
            v = state.v;
            if (step % opt.exportEvery == 0) {
                if (exporting) {
                    fiction::writeSurfaceObj(opt.exportDir + "/" + std::to_string(step) + ".obj",
                                             x, scene.mesh.surface);
                }
                if (rendering) renderFrame(x, step * opt.dt);
            }
        }
        const auto t1 = std::chrono::steady_clock::now();
        const double secs = std::chrono::duration<double>(t1 - t0).count();

        // 位移与速度统计
        const fiction::Positions dx = x - scene.model.x0;
        double maxDisp = 0.0;
        int argMaxDisp = 0;
        for (int i = 0; i < x.rows(); ++i) {
            const double d = dx.row(i).norm();
            if (d > maxDisp) {
                maxDisp = d;
                argMaxDisp = i;
            }
        }
        const double maxSpeed = v.cwiseAbs().maxCoeff();
        const bool finite = x.allFinite() && v.allFinite();

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "仿真完成 : " << opt.steps << " 步, 用时 " << secs << " s ("
                  << secs / std::max(opt.steps, 1) * 1000.0 << " ms/步)\n"
                  << "最大位移 : " << maxDisp << " (node " << argMaxDisp << ")\n"
                  << "最大速度 : " << maxSpeed << "\n";
        if (exporting) {
            std::cout << "OBJ 导出 : " << opt.exportDir << "/（每 " << opt.exportEvery
                      << " 步一帧，共 " << opt.steps / opt.exportEvery + 1 << " 帧）\n";
        }

        // 视频编码
        if (rendering) {
            const std::string pattern = fiction::ffmpegFramePattern(framesDir, ext);
            std::cout << "视频渲染 : " << videoW << "x" << videoH << ", " << frameIndex
                      << " 帧, " << opt.videoFps << " fps\n";
            if (frameIndex <= 1) {
                std::cout << "视频帧数不足（<2），跳过编码，保留帧目录: " << framesDir << "\n";
            } else if (fiction::encodeVideo(pattern, opt.videoFps, opt.videoPath)) {
                std::cout << "视频输出 : " << opt.videoPath << "\n";
                if (!opt.videoKeepFrames) {
                    std::error_code ec;
                    std::filesystem::remove_all(framesDir, ec);
                }
            } else {
                std::cerr << "[ipc_mesh_demo] ffmpeg 编码失败（或未安装），保留帧目录: " << framesDir
                          << "\n  手工编码: ffmpeg -y -framerate " << opt.videoFps << " -i \""
                          << pattern << "\" -c:v libx264 -pix_fmt yuv420p -crf 20 \""
                          << opt.videoPath << "\"\n";
            }
        }

        std::cout << (finite ? "[PASS] " : "[FAIL] ") << "结果有限（无 NaN/Inf）\n";
        return finite ? 0 : 1;
    } catch (const std::exception& ex) {
        std::cerr << "[ipc_mesh_demo] 错误: " << ex.what() << "\n";
        return 1;
    }
}
