// =============================================================================
// examples/simulator.cpp — 控制台验证（纯 3D 版）
// =============================================================================
// 职责：
//   - 演示 TimeIntegrator + Inertia 内置 + IEnergy 扩展（纯 3D，固定 [N,3]）
//   - 仅控制台打印，不导出 .obj，不依赖可视化
//   - 展示：
//       1) 纯惯性 3D（验证 Inertia + 线性求解 + 固定 3D 索引）
//       2) 重力 3D（验证 z 方向）
//       3) DBC 3D 驱动（Vec3）
//       4) 自定义二次势（纯 3D）
//       5) 2D 退化回归：upgrade2DTo3D + lockZ 后结果与原 2D 逻辑一致，z 保持 0
// =============================================================================

#include "time_integrator.h"
#include "energies/Gravity.h"
#include "energies/Gravitation.h"
#include "energies/zlock_energy.h"
#include "types.h"
#include "utils.h"

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>

// -----------------------------------------------------------------------------
// 辅助：打印 [N,3] 位置与速度（截断大 N）
// -----------------------------------------------------------------------------
static void printState(const fiction::Positions& x,
                       const fiction::Velocities& v,
                       const std::string& label,
                       int maxRows = 6) {
    const int N = static_cast<int>(x.rows());
    std::cout << label << " N=" << N << " Dim=" << fiction::kDim << "\n";
    const int rows = std::min(N, maxRows);
    std::cout << std::setprecision(6) << std::fixed;
    for (int i = 0; i < rows; ++i) {
        std::cout << "  node " << i << ": x=(" << x(i,0) << ", " << x(i,1) << ", " << x(i,2)
                  << ") v=(" << v(i,0) << ", " << v(i,1) << ", " << v(i,2) << ")\n";
    }
    if (N > rows) std::cout << "  ... (" << N - rows << " more nodes)\n";
}

// -----------------------------------------------------------------------------
// 场景 1：纯惯性 3D
// -----------------------------------------------------------------------------
static void demoPureInertia3D() {
    std::cout << "\n===== Demo 1: 纯惯性 3D (固定 [N,3]) =====\n";
    const int N = 2;
    fiction::Positions x(N,3), v(N,3);
    x << 0,0,0,
         1,0,1;
    v << 1,0,0.5,
         0,1,-0.5;
    Eigen::VectorXd masses(N); masses << 1.0, 2.0;

    fiction::TimeIntegrator::Config cfg; cfg.h = 0.01; cfg.tol = 1e-8; cfg.verbose = true;
    fiction::TimeIntegrator integrator(cfg, masses);
    fiction::DBCInfo dbc;

    std::cout << "初态:\n"; printState(x, v, " before");
    auto state = integrator.stepForward(x, v, dbc);
    std::cout << "单步后 (期望 x ≈ x+v*h):\n"; printState(state.x, state.v, " after");
    fiction::Positions expected = x + v * cfg.h;
    double err = (state.x - expected).lpNorm<Eigen::Infinity>();
    std::cout << "max |x_new - (x+v*h)| = " << err << " " << (err<1e-10?"[PASS]":"[FAIL]") << "\n";
}

// -----------------------------------------------------------------------------
// 场景 2：重力 3D
// -----------------------------------------------------------------------------
static void demoGravity3D() {
    std::cout << "\n===== Demo 2: 重力 3D (ExampleGravityEnergy, g=(0,0,-9.81)) =====\n";
    const int N = 3;
    fiction::Positions x(N,3), v(N,3);
    x << 0,1,0,
         1,1,0,
         0.5,2,0;
    v.setZero();
    Eigen::VectorXd masses(N); masses << 1.0, 1.0, 1.0;

    fiction::TimeIntegrator::Config cfg; cfg.h = 0.01; cfg.tol = 1e-7; cfg.verbose = true;
    fiction::TimeIntegrator integrator(cfg, masses);
    auto gravity = std::make_shared<fiction::ExampleGravityEnergy>(masses, fiction::Vec3(0,0,-9.81));
    integrator.addEnergy(gravity);
    std::cout << "已注册: " << gravity->name() << " 重量 h*h=" << gravity->hessianWeight(cfg.h) << "\n";

    fiction::DBCInfo dbc;
    std::cout << "初态:\n"; printState(x, v, " before");
    auto state = integrator.stepForward(x, v, dbc);
    std::cout << "单步后:\n"; printState(state.x, state.v, " after");
    double dz = state.x(0,2) - x(0,2);
    double dx = state.x(0,0) - x(0,0);
    double dy = state.x(0,1) - x(0,1);
    std::cout << "node0 dz=" << dz << " dx=" << dx << " dy=" << dy
              << " dz<0 && |dx|,|dy|≈0 " << ((dz<0 && std::abs(dx)<1e-9 && std::abs(dy)<1e-9)?"[PASS]":"[FAIL]") << "\n";
}

// -----------------------------------------------------------------------------
// 场景 2.5：万有引力 3D
// -----------------------------------------------------------------------------
static void demoGravitation3D() {
    std::cout << "\n===== Demo 2.5: 万有引力 3D (GravitationEnergy, 地球绕太阳) =====\n";
    const int N = 1;
    fiction::Positions x(N,3), v(N,3);
    x << 1.496e8, 0, 0;   // 日地距离 (m)
    v.setZero();

    const double mSun   = 1.989e30;   // 太阳质量 (kg)
    const double mEarth = 5.965e24;   // 地球质量 (kg)
    Eigen::VectorXd masses(N); masses << mEarth;   // 惯性能量节点质量

    // 注意：坐标为 1.5e8 m 量级，双精度 ulp≈3e-8，位置更新无法低于该精度，
    // 故收敛容差取 1e-4（残余 1e-6 远小于它，2~3 次迭代即收敛）
    fiction::TimeIntegrator::Config cfg; cfg.h = 0.01; cfg.tol = 1e-4; cfg.verbose = true;
    fiction::TimeIntegrator integrator(cfg, masses);

    // 太阳固定于原点，作为点质量吸引子（质量 + 位置）
    fiction::GravitationEnergy::Attractor sun{mSun, fiction::Vec3(0,0,0)};
    auto grav = std::make_shared<fiction::GravitationEnergy>(
        masses, std::vector<fiction::GravitationEnergy::Attractor>{sun});
    integrator.addEnergy(grav);
    std::cout << "已注册: " << grav->name() << " 权重 h*h=" << grav->hessianWeight(cfg.h) << "\n";

    fiction::DBCInfo dbc;
    std::cout << "初态:\n"; printState(x, v, " before");
    auto state = integrator.stepForward(x, v, dbc);
    std::cout << "单步后:\n"; printState(state.x, state.v, " after");

    // 校验：仅受向心引力，地球应向原点（-x 方向）移动，y/z 保持 0
    double dx = state.x(0,0) - x(0,0);
    double dy = state.x(0,1) - x(0,1);
    double dz = state.x(0,2) - x(0,2);
    // 隐式欧拉 v=0：x_new = x + h*v_new, v_new = h*a  =>  dx ≈ -h^2 * G*M_sun/r^2
    const double accel = 6.674e-11 * mSun / x.row(0).squaredNorm();
    const double expectDx = -accel * cfg.h * cfg.h;
    std::cout << "node0 dx=" << dx << " (期望 ≈ " << expectDx << ") dy=" << dy << " dz=" << dz << "\n";
    std::cout << "向原点 dx<0 && |dy|,|dz|≈0 "
              << (dx < 0 && std::abs(dy) < 1e-6 && std::abs(dz) < 1e-6 ? "[PASS]" : "[FAIL]") << "\n";
    std::cout << "|dx - expectDx| < 1e-6 "
              << (std::abs(dx - expectDx) < 1e-6 ? "[PASS]" : "[FAIL]") << "\n";
}

// -----------------------------------------------------------------------------
// 场景 3：DBC 3D 驱动
// -----------------------------------------------------------------------------
static void demoDBC3D() {
    std::cout << "\n===== Demo 3: DBC 3D 驱动 (Vec3) =====\n";
    const int N = 2;
    fiction::Positions x(N,3), v(N,3);
    x << 0,0,0,
         0,1,0;
    v.setZero();
    Eigen::VectorXd masses(N); masses << 1.0, 1.0;

    fiction::TimeIntegrator::Config cfg; cfg.h = 0.01; cfg.tol = 1e-2; cfg.verbose = true;
    fiction::TimeIntegrator integrator(cfg, masses);

    fiction::DBCInfo dbc;
    dbc.indices = {1};
    dbc.isDBC = std::vector<char>(N, 0); dbc.isDBC[1]=1;
    dbc.stiffness = 2000;
    dbc.targetVelocities = { fiction::Vec3(0,0,-1.0) };
    dbc.limits           = { fiction::Vec3(0,0,-2.0) };

    std::cout << "初态 x[DBC=1]=(0,1,0) DBC vel=(0,0,-1)\n";
    printState(x, v, " before");

    fiction::Positions curX = x; fiction::Velocities curV = v;
    for (int step=0; step<3; ++step) {
        std::cout << "\n--- Step " << step << " ---\n";
        auto state = integrator.stepForward(curX, curV, dbc);
        curX = state.x; curV = state.v;
        std::cout << " after step " << step << ": DBC node x=(" << curX(1,0) << "," << curX(1,1)
                  << "," << curX(1,2) << ") v=(" << curV(1,0) << "," << curV(1,1) << "," << curV(1,2) << ")\n";
    }
}

// -----------------------------------------------------------------------------
// 场景 4：自定义能量（纯 3D）
// -----------------------------------------------------------------------------
static void demoCustomEnergy() {
    std::cout << "\n===== Demo 4: 自定义能量接口（纯 3D） =====\n";
    struct QuadraticEnergy final : public fiction::IEnergy {
        double w; int N;
        QuadraticEnergy(double w_, int N_): w(w_), N(N_) {}
        const char* name() const override { return "QuadraticDemo"; }
        double hessianWeight(double h) const override { return h*h; }
        double val(const fiction::Positions& x) const override {
            double E=0;
            for(int i=0;i<N && i<(int)x.rows();++i) {
                fiction::Vec3 xi(x(i,0),x(i,1),x(i,2));
                E+= 0.5*w*xi.squaredNorm();
            }
            return E;
        }
        void addGrad(const fiction::Positions& x, fiction::Positions& grad) const override {
            (void)x;
            for(int i=0;i<N && i<(int)grad.rows();++i)
                grad.row(i) += w * x.row(i);
        }
        void addHess(const fiction::Positions& x, std::vector<fiction::Triplet>& trips) const override {
            (void)x;
            for(int i=0;i<N;++i) for(int d=0; d<3; ++d) trips.emplace_back(3*i+d, 3*i+d, w);
        }
    };

    const int N=2;
    fiction::Positions x(N,3), v(N,3);
    x << 1,0,0, 0,1,1; v.setZero();
    Eigen::VectorXd masses(N); masses<<1,1;
    fiction::TimeIntegrator::Config cfg; cfg.h=0.01; cfg.tol=1e-8; cfg.verbose=false;
    fiction::TimeIntegrator integrator(cfg, masses);
    integrator.addEnergy(std::make_shared<QuadraticEnergy>(10.0, N));
    fiction::DBCInfo dbc;
    auto state = integrator.stepForward(x, v, dbc);
    std::cout << " 3D Quadratic: |x| " << x.norm() << " -> " << state.x.norm()
              << (state.x.norm() < x.norm() ? " [PASS]" : " [FAIL]")
              << " Dim=" << state.x.cols() << "\n";
}

// -----------------------------------------------------------------------------
// 场景 5：2D 退化回归（upgrade2DTo3D + ZLockEnergy 锁 z）
//   与原 2D 算例等价：z 恒为 0，x/y 位移与原 2D 结果一致
// -----------------------------------------------------------------------------
static void demo2DCompat() {
    std::cout << "\n===== Demo 5: 2D 退化兼容 (upgrade2DTo3D + ZLockEnergy) =====\n";
    // 原 2D 场景：重力 2D（g=(0,-9.81)），N=3
    const int N = 3;
    Eigen::MatrixXd x2(N,2);
    x2 << 0,1,
          1,1,
          0.5,2;
    // 2D → 3D 补零
    fiction::Positions x3 = fiction::upgrade2DTo3D(x2);
    fiction::Velocities v3 = fiction::Velocities::Zero(N, 3);
    Eigen::VectorXd masses(N); masses << 1.0, 1.0, 1.0;

    fiction::TimeIntegrator::Config cfg; cfg.h = 0.01; cfg.tol = 1e-7; cfg.verbose = true;
    fiction::TimeIntegrator integrator(cfg, masses);

    // 2D 重力：y 方向 -9.81（在 3D 中映射到 y 轴，与 old 2D 一致）
    auto gravity = std::make_shared<fiction::ExampleGravityEnergy>(masses, fiction::Vec3(0,-9.81,0));
    integrator.addEnergy(gravity);
    // 锁 z 轴：仅罚 z 分量，不影响 x/y 运动
    auto zlock = fiction::makeZLockEnergy(masses);
    integrator.addEnergy(zlock);
    std::cout << "已注册: " << gravity->name() << " + " << zlock->name() << "\n";

    std::cout << "初态 (z=0):\n"; printState(x3, v3, " before");
    fiction::DBCInfo dbc;
    auto state = integrator.stepForward(x3, v3, dbc);
    std::cout << "单步后:\n"; printState(state.x, state.v, " after");

    // 校验 1：z 保持 0
    double maxAbsZ = state.x.col(2).cwiseAbs().maxCoeff();
    // 校验 2：y 方向位移 ≈ g*h²（隐式欧拉：x_new = x + h*v_new, v_new = v + h*g）
    double dy = state.x(0,1) - x3(0,1);
    double expectDy = (-9.81) * cfg.h * cfg.h;
    std::cout << "max|z| = " << maxAbsZ << " (应≈0) "
              << (maxAbsZ < 1e-9 ? "[PASS]" : "[FAIL]") << "\n";
    std::cout << "node0 dy = " << dy << " (期望 ≈ " << expectDy << ") "
              << (std::abs(dy - expectDy) < 1e-6 ? "[PASS]" : "[FAIL]") << "\n";
}

int main() {
    std::cout << "fiction C++ TimeIntegrator 控制台验证 (纯 3D, Eigen+CHOLMOD)\n";
    std::cout << "编译信息: C++" << __cplusplus << "  Eigen " << EIGEN_WORLD_VERSION << "."
              << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << "\n";
#ifdef FICTION_HAS_CHOLMOD
    std::cout << "求解器: Eigen::CholmodSupernodalLLT + SimplicialLDLT 回退\n";
#else
    std::cout << "求解器: Eigen::SimplicialLDLT (纯 Eigen)\n";
#endif
    std::cout << "框架：Inertia 内置，IEnergy 纯 3D [N,3]，索引 3*node+d\n";
    std::cout << "2D 适配：upgrade2DTo3D 补零 + lockZAxis 锁 z\n";

    demoPureInertia3D();
    demoGravity3D();
    demoGravitation3D();
    demoDBC3D();
    demoCustomEnergy();
    demo2DCompat();

    std::cout << "\n===== 全部演示完成 (纯 3D + 2D 退化兼容) =====\n";
    std::cout << "下一步：\n";
    std::cout << "  - 新能量继承 IEnergy，Positions 固定 [N,3]，索引 = 3*node+d\n";
    std::cout << "  - 边单元 H_local 6x6（2 节点×3 维），四面体 12x12\n";
    std::cout << "  - DBC 使用 Vec3，2D 算例用 upgrade2DTo3D + makeZLockEnergy\n";
    return 0;
}
