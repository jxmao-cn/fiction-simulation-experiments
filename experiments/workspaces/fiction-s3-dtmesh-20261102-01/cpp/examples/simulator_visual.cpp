// =============================================================================
// examples/simulator_visual.cpp — 交互式可视化 demo（GLFW + OpenGL）
// =============================================================================
// 职责：
//   - 对标 examples/simulator.cpp 的控制台验证，改为实时窗口可视化
//   - 仅使用现有能量，不新增物理内容：
//       1) 一组带随机初速的质点 + ExampleGravityEnergy（对标 Demo 2 重力）
//       2) 一个红色 DBC 驱动节点（对标 Demo 3：目标速度向下，限位后悬停）
//       3) 地面网格仅作视觉参考（C++ 侧尚无接触/Barrier 能量，质点会穿过）
//
// 操作：
//   左键拖拽 = 旋转视角   滚轮 = 缩放   右键(或 Shift+左键)拖拽 = 平移
//   Space = 暂停/继续     R = 复位相机   Up/Down = 每帧仿真步数 ±1   ESC = 退出
// =============================================================================

#include "viewer.h"
#include "time_integrator.h"
#include "energies/Gravity.h"
#include "types.h"

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

using fiction::Vec3;

namespace {

// -----------------------------------------------------------------------------
// 场景参数
// -----------------------------------------------------------------------------
constexpr int kParticles = 200;      // 重力质点数
constexpr double kBoxHalf = 1.5;     // 初始位置立方体半宽
constexpr double kMaxSpeed = 1.5;    // 初始最大速度
constexpr double kGroundZ = 0.0;     // 地面高度（仅视觉参考）
constexpr double kH = 0.01;          // 时间步长
constexpr double kTol = 1e-2;        // Newton 容差（实时演示放宽）
constexpr int kDbcNode = kParticles; // DBC 节点在数组末尾

struct Scene {
    fiction::Positions x;
    fiction::Velocities v;
    Eigen::VectorXd masses;
    std::vector<char> isDBC;
    fiction::DBCInfo dbc;
    std::unique_ptr<fiction::TimeIntegrator> integrator;
};

Scene makeScene() {
    Scene s;
    const int N = kParticles + 1;
    s.x.setZero(N, 3);
    s.v.setZero(N, 3);
    s.masses.setOnes(N);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> pos(-kBoxHalf, kBoxHalf);
    std::uniform_real_distribution<double> vel(-kMaxSpeed, kMaxSpeed);
    for (int i = 0; i < kParticles; ++i) {
        s.x(i, 0) = pos(rng);
        s.x(i, 1) = pos(rng);
        s.x(i, 2) = kGroundZ + 2.0 + pos(rng) * 0.5;  // 初始悬浮在地面之上
        s.v(i, 0) = vel(rng);
        s.v(i, 1) = vel(rng);
        s.v(i, 2) = vel(rng);
    }

    // DBC 驱动节点：从 (0, 0, 3) 以 (0,0,-1) 下落，到 z=1.0 悬停
    s.x(kDbcNode, 0) = 0.0;
    s.x(kDbcNode, 1) = 0.0;
    s.x(kDbcNode, 2) = 3.0;
    s.isDBC.assign(N, 0);
    s.isDBC[kDbcNode] = 1;
    s.dbc.indices = {kDbcNode};
    s.dbc.targetVelocities = {fiction::Vec3(0.0, 0.0, -1.0)};
    s.dbc.limits = {fiction::Vec3(0.0, 0.0, 1.0)};
    s.dbc.stiffness = 1e6;
    s.dbc.isDBC = s.isDBC;

    fiction::TimeIntegrator::Config cfg;
    cfg.h = kH;
    cfg.tol = kTol;
    cfg.verbose = false;  // 实时窗口不刷屏
    s.integrator = std::make_unique<fiction::TimeIntegrator>(cfg, s.masses);
    s.integrator->addEnergy(
        std::make_shared<fiction::ExampleGravityEnergy>(s.masses, fiction::Vec3(0.0, 0.0, -9.81)));
    return s;
}

} // anonymous namespace

int main() {
    fiction::Viewer::Config vcfg;
    vcfg.title = "fiction visual demo (Space=暂停 R=复位相机 Up/Down=步数 ESC=退出)";
    fiction::Viewer viewer(vcfg, fiction::Vec3(0.0, 0.0, 0.8), 7.0);
    if (!viewer.ok()) {
        std::cerr << "[simulator_visual] Viewer 初始化失败\n";
        return 1;
    }

    Scene scene = makeScene();
    std::cout << "\n===== fiction 可视化 demo =====\n"
              << "蓝色质点: 重力作用 (ExampleGravity, g=(0,0,-9.81))\n"
              << "红色节点: DBC 驱动 (vel=(0,0,-1), 限位 z=1.0)\n"
              << "地面仅作视觉参考 (尚无 Barrier/接触能量, 质点会穿过)\n"
              << "左键旋转 | 滚轮缩放 | 右键平移 | Space 暂停 | R 复位相机"
              << " | Up/Down 每帧步数 | ESC 退出\n\n";

    const Vec3 groundOrigin(0.0, 0.0, kGroundZ);
    const Vec3 groundNormal(0.0, 0.0, 1.0);
    int frame = 0;

    while (!viewer.shouldClose()) {
        viewer.pollEvents();

        viewer.beginFrame();
        viewer.drawGround(groundOrigin, groundNormal, 5.0, 0.5,
                          0.35f, 0.38f, 0.42f);
        viewer.drawPoints(scene.x, 3.0f, 0.35f, 0.55f, 1.0f);      // 全部质点
        viewer.drawPoints(scene.x, scene.isDBC, 9.0f, 1.0f, 0.25f, 0.2f);  // DBC 高亮
        viewer.endFrame();

        if (!viewer.paused() && viewer.stepsPerFrame() > 0) {
            for (int k = 0; k < viewer.stepsPerFrame(); ++k) {
                auto state = scene.integrator->stepForward(scene.x, scene.v, scene.dbc);
                scene.x = state.x;
                scene.v = state.v;
            }
        }
        ++frame;
    }
    return 0;
}
