// =============================================================================
// examples/barrier_check.cpp — BarrierEnergy 数值验证（纯 3D，无可视化）
// =============================================================================
// 对标 Python/BarrierEnergy.py，验证四件事：
//   1) 有限差分校验梯度（floor 固定平面 + ceiling 节点锚定平面同时激活）
//   2) 有限差分校验海森矩阵（含 ceiling 的 (i, driver) 6×6 耦合块）
//   3) initStepSize 数值 CCD 行为（靠近/远离/锚点相对运动）
//   4) 落地仿真：重力 + Barrier，验证全程不穿透且最终停在 dhat 接触带内
// 退出码：0 = 全部 PASS，1 = 存在 FAIL（可供 ctest 捕获）
// =============================================================================

#include "time_integrator.h"
#include "energies/BarrierEnergy.h"
#include "energies/Gravity.h"
#include "types.h"
#include "utils.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

namespace {

int g_failCount = 0;

void check(bool ok, const std::string& msg) {
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << msg << "\n";
    if (!ok) ++g_failCount;
}

// -----------------------------------------------------------------------------
// 测试场景：node3 为 ceiling 锚定节点（对标 Python x[-1]，n=(0,-1,0)）
//   node0: y=0.005 → floor 激活（d=0.005 < dhat=0.01）
//   node1: y=0.020 → 双平面均不激活
//   node2: y=0.044 → ceiling 激活（d = 0.05-0.044 = 0.006）
//   node3: y=0.050 → floor 不激活；ceiling 中作为 driver 被跳过，仅受反作用
// 所有激活对的 d 距离边界 {0, dhat} ≥ 0.004，远离有限差分扰动半径
// -----------------------------------------------------------------------------
fiction::Positions makeTestPositions() {
    fiction::Positions x(4, 3);
    x << 0.1, 0.005, 0.3,
         0.4, 0.020, 0.1,
         0.2, 0.044, 0.5,
         0.7, 0.050, 0.2;
    return x;
}

Eigen::VectorXd makeTestContactArea() {
    Eigen::VectorXd A(4);
    A << 1.0, 2.0, 0.5, 1.5;
    return A;
}

std::vector<fiction::BarrierEnergy::Plane> makeTestPlanes() {
    return {
        {fiction::Vec3(0.0, 1.0, 0.0),  fiction::Vec3::Zero(), -1}, // floor：固定平面
        {fiction::Vec3(0.0, -1.0, 0.0), fiction::Vec3::Zero(),  3}, // ceiling：锚定 node3
    };
}

// -----------------------------------------------------------------------------
// 1) 有限差分校验梯度
// -----------------------------------------------------------------------------
void fdCheckGrad() {
    std::cout << "\n===== 检验 1: 有限差分 vs addGrad =====\n";
    const fiction::Positions x = makeTestPositions();
    fiction::BarrierEnergy barrier(makeTestContactArea(), makeTestPlanes());

    fiction::Positions g = fiction::Positions::Zero(4, 3);
    barrier.addGrad(x, g);

    const double eps = 1e-7;
    double maxErr = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int d = 0; d < 3; ++d) {
            fiction::Positions xp = x, xm = x;
            xp(i, d) += eps;
            xm(i, d) -= eps;
            const double fd = (barrier.val(xp) - barrier.val(xm)) / (2.0 * eps);
            maxErr = std::max(maxErr, std::abs(fd - g(i, d)));
        }
    }
    const double scale = std::max(g.cwiseAbs().maxCoeff(), 1.0);
    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "  |grad|_inf = " << g.cwiseAbs().maxCoeff()
              << "  max FD err = " << maxErr << "\n";
    check(maxErr / scale < 1e-6, "梯度相对误差 < 1e-6");
}

// -----------------------------------------------------------------------------
// 2) 有限差分校验海森矩阵
// -----------------------------------------------------------------------------
void fdCheckHess() {
    std::cout << "\n===== 检验 2: 有限差分(梯度) vs addHess =====\n";
    const fiction::Positions x = makeTestPositions();
    fiction::BarrierEnergy barrier(makeTestContactArea(), makeTestPlanes());

    // 精确海森：Triplet 组装为稠密（同一下标累加，与 setFromTriplets 一致）
    std::vector<fiction::Triplet> trips;
    barrier.addHess(x, trips);
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(12, 12);
    for (const auto& t : trips) H(t.row(), t.col()) += t.value();

    // FD 海森：对梯度做中心差分 H_fd(:,k) = (g(x+eps*e_k) - g(x-eps*e_k)) / (2*eps)
    const double eps = 1e-6;
    Eigen::MatrixXd Hfd = Eigen::MatrixXd::Zero(12, 12);
    for (int k = 0; k < 12; ++k) {
        fiction::Positions xp = x, xm = x;
        xp(k / 3, k % 3) += eps;
        xm(k / 3, k % 3) -= eps;
        fiction::Positions gp = fiction::Positions::Zero(4, 3);
        fiction::Positions gm = fiction::Positions::Zero(4, 3);
        barrier.addGrad(xp, gp);
        barrier.addGrad(xm, gm);
        // flatten 为 RowMajor 顺序 [x0,y0,z0, x1,...]，与 DOF 索引 3*node+d 一致
        Hfd.col(k) = (fiction::flatten(gp) - fiction::flatten(gm)) / (2.0 * eps);
    }

    const double denom = std::max(H.cwiseAbs().maxCoeff(), 1.0);
    const double err = (H - Hfd).cwiseAbs().maxCoeff();
    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "  |H|_inf = " << H.cwiseAbs().maxCoeff()
              << "  max FD err = " << err << "\n";
    check(err / denom < 1e-6, "海森相对误差 < 1e-6");
    // ceiling 耦合块存在性：H(3*2+1, 3*3+1) 与 H(3*3+1, 3*2+1) 非零
    check(H(3 * 2 + 1, 3 * 3 + 1) != 0.0 && H(3 * 3 + 1, 3 * 2 + 1) != 0.0,
          "ceiling (node2, driver) 6x6 耦合块非零");
}

// -----------------------------------------------------------------------------
// 3) initStepSize 数值 CCD 行为
// -----------------------------------------------------------------------------
void checkInitStepSize() {
    std::cout << "\n===== 检验 3: initStepSize（数值 CCD） =====\n";
    const fiction::Positions x = makeTestPositions();
    fiction::BarrierEnergy barrier(makeTestContactArea(), makeTestPlanes());

    // 场景 A：node0 高度 0.005 向下运动 → alpha ≈ 0.9*0.005/1 = 0.0045
    {
        fiction::Positions p = fiction::Positions::Zero(4, 3);
        p.row(0) << 0.0, -1.0, 0.0;
        const double alpha = barrier.initStepSize(x, p);
        const double expect = 0.9 * 0.005 / 1.0;
        check(std::abs(alpha - expect) < 1e-12,
              "朝地面运动 alpha=" + std::to_string(alpha) + " ≈ " + std::to_string(expect));
        check((x + alpha * p)(0, 1) > 0.0, "步后 node0 高度 > 0（不穿透）");
    }
    // 场景 B/D：仅 floor 平面（floor+ceiling 对置时，任何竖直运动必朝向其一，
    //   CCD 恒会限步——这是对标 Python 的保守防隧穿行为，无 dhat 过滤）
    {
        const std::vector<fiction::BarrierEnergy::Plane> floorPlane = {
            {fiction::Vec3(0.0, 1.0, 0.0), fiction::Vec3::Zero(), -1}
        };
        fiction::BarrierEnergy floorBarrier(makeTestContactArea(), floorPlane);
        // 场景 B：node0 向上远离地面 → alpha = 1
        fiction::Positions p = fiction::Positions::Zero(4, 3);
        p.row(0) << 0.0, 1.0, 0.0;
        check(floorBarrier.initStepSize(x, p) == 1.0, "远离平面 alpha = 1");
        // 场景 D：远处（d=0.044 > dhat）向地面运动同样受限（保守防隧穿）
        fiction::Positions pFar = fiction::Positions::Zero(4, 3);
        pFar.row(2) << 0.0, -1.0, 0.0;
        const double alphaFar = floorBarrier.initStepSize(x, pFar);
        check(std::abs(alphaFar - 0.9 * 0.044) < 1e-12,
              "远处逼近也受限 alpha=" + std::to_string(alphaFar)
                  + " ≈ " + std::to_string(0.9 * 0.044));
    }
    // 场景 C：ceiling 锚点 node3 以 (0,-2,0) 下压 node2（d=0.006）
    //   相对速度 n·(p2 - p3) = (0,-1,0)·(0,2,0) = -2 → alpha ≈ 0.9*0.006/2 = 0.0027
    {
        fiction::Positions p = fiction::Positions::Zero(4, 3);
        p.row(3) << 0.0, -2.0, 0.0;
        const double alpha = barrier.initStepSize(x, p);
        const double expect = 0.9 * 0.006 / 2.0;
        check(std::abs(alpha - expect) < 1e-12,
              "ceiling 相对逼近 alpha=" + std::to_string(alpha) + " ≈ " + std::to_string(expect));
        const fiction::Positions xt = x + alpha * p;
        const double dAfter = -(xt(2, 1) - xt(3, 1));
        check(dAfter > 0.0, "步后 node2-锚点间距 > 0");
    }
}

// -----------------------------------------------------------------------------
// 4) 落地仿真：重力 + Barrier，全程不穿透且停在 dhat 接触带内
// -----------------------------------------------------------------------------
void dropSimulation() {
    std::cout << "\n===== 检验 4: 落地仿真（重力 + Barrier，300 步） =====\n";
    const int N = 2;
    fiction::Positions x(N, 3), v(N, 3);
    x << 0.0, 0.30, 0.0,
         0.5, 0.50, 0.0;
    v.setZero();
    const Eigen::VectorXd masses = Eigen::VectorXd::Constant(N, 1.0);
    const Eigen::VectorXd A = Eigen::VectorXd::Constant(N, 1.0);

    fiction::TimeIntegrator::Config cfg;
    cfg.h = 0.01;
    cfg.tol = 1e-6;
    cfg.verbose = false;
    fiction::TimeIntegrator integrator(cfg, masses);
    integrator.addEnergy(std::make_shared<fiction::ExampleGravityEnergy>(
        masses, fiction::Vec3(0.0, -9.81, 0.0)));
    std::vector<fiction::BarrierEnergy::Plane> planes = {
        {fiction::Vec3(0.0, 1.0, 0.0), fiction::Vec3::Zero(), -1} // 地面 y=0
    };
    integrator.addEnergy(std::make_shared<fiction::BarrierEnergy>(A, planes));

    double minH = 1e9;
    fiction::Positions curX = x, curV = v;
    for (int s = 0; s < 300; ++s) {
        const auto st = integrator.stepForward(curX, curV, fiction::DBCInfo{});
        curX = st.x;
        curV = st.v;
        for (int i = 0; i < N; ++i) minH = std::min(minH, curX(i, 1));
    }

    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "  最终高度: node0 y=" << curX(0, 1) << "  node1 y=" << curX(1, 1)
              << "  全程最小高度=" << minH << "\n";
    check(minH > 0.0, "全程未穿透（min height > 0）");
    const double finalMinH = std::min(curX(0, 1), curX(1, 1));
    check(finalMinH > 0.0 && finalMinH < 0.01,
          "最终停留在 dhat=0.01 接触带内");
    const double vmax = curV.cwiseAbs().maxCoeff();
    std::cout << "  最终 max|v| = " << vmax << "\n";
    check(vmax < 0.1, "最终速度趋近于零（稳定接触）");
}

} // namespace

int main() {
    std::cout << "fiction C++ BarrierEnergy 数值验证 (纯 3D, 对标 Python/BarrierEnergy.py)\n";
    fdCheckGrad();
    fdCheckHess();
    checkInitStepSize();
    dropSimulation();
    std::cout << "\n===== BarrierEnergy 验证 " << (g_failCount == 0 ? "全部 PASS" : "存在 FAIL")
              << " =====\n";
    return g_failCount == 0 ? 0 : 1;
}
