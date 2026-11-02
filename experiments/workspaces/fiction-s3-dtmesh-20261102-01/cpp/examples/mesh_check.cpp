// =============================================================================
// examples/mesh_check.cpp — IPC .msh 加载 + mass-spring 数值验证（纯 3D）
// =============================================================================
// 验证内容：
//   1) 解析器：cube.msh 规模/体积/包围盒；无 $Surface 时的边界提取（封闭性）
//   2) 解析器错误路径：二进制、非 4.1 版本、缺失 $Nodes 均应报错
//   3) 自由落体刚性平移：重力 + 弹簧、无 DBC，单步位移应精确 = g*h^2
//   4) 有限差分：ExampleSpringEnergy 的梯度/海森（真实网格扰动位形）
//   5) 有限差分：RegionForceEnergy 的梯度与外力的零海森
//   6) 仿真：pin-box 固定 + 区域外力拉伸 / pin-box 固定 + 重力下垂
// 退出码：0 = 全部 PASS，1 = 存在 FAIL（供 ctest 捕获）
// =============================================================================

#include "energies/NeoHookeanTetEnergy.h"
#include "energies/RegionForce.h"
#include "energies/SpringEnergy.h"
#include "io/msh_loader.h"
#include "scene/mass_spring_builder.h"
#include "scene/ipc_scene.h"
#include "time_integrator.h"
#include "types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

namespace {

int g_failCount = 0;

void check(bool ok, const std::string& msg) {
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << msg << "\n";
    if (!ok) ++g_failCount;
}

/// 测试数据目录（CMake 注入绝对路径；手工编译时回退相对路径）
std::string dataPath(const std::string& name) {
#ifdef FICTION_TEST_DATA_DIR
    return std::string(FICTION_TEST_DATA_DIR) + "/" + name;
#else
    return "tests/data/" + name;
#endif
}

/// 期望某段代码抛出 std::exception
template <typename Fn>
bool throwsException(Fn&& fn, std::string& what) {
    try {
        fn();
    } catch (const std::exception& ex) {
        what = ex.what();
        return true;
    }
    return false;
}

// -----------------------------------------------------------------------------
// 1) 解析器：正常加载
// -----------------------------------------------------------------------------
void checkLoaderBasics() {
    std::cout << "\n===== 检验 1: msh 解析器（cube.msh） =====\n";
    const fiction::MshMesh mesh = fiction::loadMsh(dataPath("cube.msh"));
    check(mesh.x.rows() == 8, "节点数 = 8（实际 " + std::to_string(mesh.x.rows()) + "）");
    check(mesh.numTets() == 6, "四面体数 = 6（实际 " + std::to_string(mesh.numTets()) + "）");
    check(mesh.numSurface() == 12, "表面三角面数 = 12（实际 " + std::to_string(mesh.numSurface()) + "）");

    const fiction::MassSpringModel model = fiction::buildMassSpringModel(mesh, 1000.0, 1e5);
    check(std::abs(model.volume - 1.0) < 1e-12,
          "体积 = 1.0（实际 " + std::to_string(model.volume) + "）");
    // 该四面体分解下 cube 的唯一边为 19 条（非几何上全部 28 条，部分面对角线不在弹簧网中）
    check(model.numEdges() == 19, "唯一边数 = 19（实际 " + std::to_string(model.numEdges()) + "）");
    const std::array<double, 6> expectBox = {0, 0, 0, 1, 1, 1};
    bool boxOk = true;
    for (int d = 0; d < 6; ++d)
        boxOk = boxOk && std::abs(model.bbox[static_cast<size_t>(d)] - expectBox[static_cast<size_t>(d)]) < 1e-12;
    check(boxOk, "包围盒 = [0,0,0]-[1,1,1]");
    check(std::abs(model.totalMass - 1000.0) < 1e-9,
          "总质量 = 1000 kg（实际 " + std::to_string(model.totalMass) + "）");
}

// -----------------------------------------------------------------------------
// 2) 解析器：无 $Surface 时边界提取 + 错误路径
// -----------------------------------------------------------------------------
void checkLoaderFallbackAndErrors() {
    std::cout << "\n===== 检验 2: $Surface 缺失回退与错误路径 =====\n";
    const fiction::MshMesh mesh = fiction::loadMsh(dataPath("cube_no_surface.msh"));
    check(mesh.numSurface() == 12, "边界提取面数 = 12（实际 " + std::to_string(mesh.numSurface()) + "）");

    // 封闭性：每条无向边恰好被 2 个边界三角面共享
    std::map<std::pair<int, int>, int> edgeCount;
    for (const auto& f : mesh.surface) {
        for (int k = 0; k < 3; ++k) {
            int a = f[k], b = f[(k + 1) % 3];
            if (a > b) std::swap(a, b);
            ++edgeCount[{a, b}];
        }
    }
    bool closed = true;
    for (const auto& [e, c] : edgeCount) closed = closed && (c == 2);
    check(closed, "边界提取结果封闭（每条边被 2 个面共享）");

    // 错误路径：二进制 flag
    {
        std::istringstream bad("$MeshFormat\n4.1 1 8\n$EndMeshFormat\n");
        std::string what;
        check(throwsException([&] { fiction::loadMsh(bad, "<binary>"); }, what), "二进制 MSH 报错");
    }
    // 错误路径：非 4.1 版本
    {
        std::istringstream bad("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
        std::string what;
        check(throwsException([&] { fiction::loadMsh(bad, "<v2.2>"); }, what), "非 4.1 版本报错");
    }
    // 错误路径：缺少 $Elements（无四面体）
    {
        std::istringstream bad(
            "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
            "$Nodes\n1 1 1 1\n3 0 0 1\n1\n0 0 0\n$EndNodes\n");
        std::string what;
        check(throwsException([&] { fiction::loadMsh(bad, "<no-elements>"); }, what),
              "缺少 $Elements 报错");
    }
}

// -----------------------------------------------------------------------------
// 3) 自由落体刚性平移：x_new = x0 + g*h^2（弹簧在平移零空间上不受力）
// -----------------------------------------------------------------------------
void checkFreeFall() {
    std::cout << "\n===== 检验 3: 自由落体刚性平移 =====\n";
    const fiction::MshMesh mesh = fiction::loadMsh(dataPath("cube.msh"));
    const fiction::MassSpringModel model = fiction::buildMassSpringModel(mesh, 1000.0, 0.0);

    fiction::TimeIntegrator::Config cfg;
    cfg.h = 0.01;
    cfg.tol = 1e-10;
    cfg.maxIter = 50;
    cfg.verbose = false;
    fiction::TimeIntegrator integrator(cfg, model.masses);
    const fiction::Vec3 g(0.0, 0.0, -9.81);
    integrator.addEnergy(std::make_shared<fiction::ExampleGravityEnergy>(model.masses, g));
    integrator.addEnergy(std::make_shared<fiction::ExampleSpringEnergy>(
        model.edges, model.l2, model.k));

    fiction::Positions x = model.x0;
    const fiction::Velocities v = fiction::Velocities::Zero(model.numNodes(), fiction::kDim);
    const auto state = integrator.stepForward(x, v, fiction::DBCInfo{});

    const fiction::Positions expected = model.x0.rowwise() + (cfg.h * cfg.h * g).transpose();
    const double err = (state.x - expected).cwiseAbs().maxCoeff();
    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "  max |x_new - (x0 + g*h^2)| = " << err << "\n";
    check(err < 1e-12, "单步位移精确等于 g*h^2");
}

// -----------------------------------------------------------------------------
// 4) ExampleSpringEnergy 有限差分（真实网格 + 拉伸位形）
//    注意：ExampleSpringEnergy 对压缩态使用 makePSD（投影牛顿设计），
//    精确海森在压缩时非 PSD，FD 无法对照；这里整体放大 1.1 倍保证所有弹簧
//    处于拉伸态（精确海森 PSD，makePSD 为恒等），FD 与解析式应完全一致。
// -----------------------------------------------------------------------------
void checkSpringFiniteDifference() {
    std::cout << "\n===== 检验 4: SpringEnergy 有限差分（cube 网格拉伸位形） =====\n";
    const fiction::MshMesh mesh = fiction::loadMsh(dataPath("cube.msh"));
    const fiction::MassSpringModel model = fiction::buildMassSpringModel(mesh, 1000.0, 1e4);
    fiction::ExampleSpringEnergy spring(model.edges, model.l2, model.k);

    // 以包围盒中心为基准整体放大 1.1 倍 + 小扰动
    const fiction::Vec3 center = fiction::boxCenter(model.bbox);
    fiction::Positions x = model.x0;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < x.rows(); ++i) {
        x.row(i) = (center + 1.1 * (x.row(i).transpose() - center)).transpose();
        for (int d = 0; d < fiction::kDim; ++d)
            x(i, d) += 0.01 * model.avgEdgeLength * dist(rng);
    }

    const int N = model.numNodes();
    const int dof = N * fiction::kDim;

    // 解析梯度
    fiction::Positions grad = fiction::Positions::Zero(N, fiction::kDim);
    spring.addGrad(x, grad);
    const fiction::FlatVector gAnalytic = fiction::flatten(grad);

    // 解析海森
    std::vector<fiction::Triplet> trips;
    spring.addHess(x, trips);
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(dof, dof);
    for (const auto& t : trips) H(t.row(), t.col()) += t.value();

    // FD 梯度 + FD 海森
    const double eps = 1e-7;
    fiction::FlatVector gFd(dof);
    Eigen::MatrixXd Hfd = Eigen::MatrixXd::Zero(dof, dof);
    for (int k = 0; k < dof; ++k) {
        fiction::Positions xp = x, xm = x;
        xp(k / fiction::kDim, k % fiction::kDim) += eps;
        xm(k / fiction::kDim, k % fiction::kDim) -= eps;
        const double fp = spring.val(xp), fm = spring.val(xm);
        gFd[k] = (fp - fm) / (2.0 * eps);

        fiction::Positions gp = fiction::Positions::Zero(N, fiction::kDim);
        fiction::Positions gm = fiction::Positions::Zero(N, fiction::kDim);
        spring.addGrad(xp, gp);
        spring.addGrad(xm, gm);
        Hfd.col(k) = (fiction::flatten(gp) - fiction::flatten(gm)) / (2.0 * eps);
    }

    const double gradErr = (gAnalytic - gFd).cwiseAbs().maxCoeff();
    const double gradScale = std::max(gAnalytic.cwiseAbs().maxCoeff(), 1.0);
    const double hessErr = (H - Hfd).cwiseAbs().maxCoeff();
    const double hessScale = std::max(H.cwiseAbs().maxCoeff(), 1.0);

    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "  |grad|_inf=" << gradScale << "  max FD err=" << gradErr << "\n"
              << "  |H|_inf=" << hessScale << "  max FD err=" << hessErr << "\n";
    check(gradErr / gradScale < 1e-6, "梯度相对误差 < 1e-6");
    check(hessErr / hessScale < 1e-6, "海森相对误差 < 1e-6");
    check((H - H.transpose()).cwiseAbs().maxCoeff() / hessScale < 1e-10, "海森对称");
}

// -----------------------------------------------------------------------------
// 5) RegionForceEnergy 有限差分 + 零海森
// -----------------------------------------------------------------------------
void checkRegionForceFiniteDifference() {
    std::cout << "\n===== 检验 5: RegionForceEnergy 有限差分 =====\n";
    const fiction::MshMesh mesh = fiction::loadMsh(dataPath("cube.msh"));
    const fiction::MassSpringModel model = fiction::buildMassSpringModel(mesh, 1000.0, 0.0);

    // 选 x > 0.5 的节点区域，加速度 (2.5, -1.0, 0.75)
    std::vector<int> nodes;
    for (int i = 0; i < model.numNodes(); ++i)
        if (model.x0(i, 0) > 0.5) nodes.push_back(i);
    check(!nodes.empty(), "区域选点非空（实际 " + std::to_string(nodes.size()) + " 节点）");

    const fiction::Vec3 accel(2.5, -1.0, 0.75);
    fiction::RegionForceEnergy force(model.masses, nodes, accel);

    const int N = model.numNodes();
    fiction::Positions x = model.x0;
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    for (int i = 0; i < N; ++i)
        for (int d = 0; d < fiction::kDim; ++d) x(i, d) += 0.1 * dist(rng);

    fiction::Positions grad = fiction::Positions::Zero(N, fiction::kDim);
    force.addGrad(x, grad);
    const fiction::FlatVector gA = fiction::flatten(grad);

    // 解析：grad_i = -m_i * a（区域外为 0）
    bool gradExact = true;
    for (int i = 0; i < N; ++i) {
        const bool inRegion = std::find(nodes.begin(), nodes.end(), i) != nodes.end();
        const fiction::Vec3 applied = -model.masses[i] * accel;
        const fiction::Vec3 expect = inRegion ? applied : fiction::Vec3::Zero();
        gradExact = gradExact && (expect - grad.row(i).transpose()).norm() < 1e-14;
    }
    check(gradExact, "梯度 = -m_i * a（区域外为 0）");

    // FD 梯度（对若干随机 DOF；val 量级 ~ m*a*x，用相对误差判定）
    const double eps = 1e-7;
    double maxErr = 0.0;
    for (int k = 0; k < N * fiction::kDim; k += 5) {
        fiction::Positions xp = x, xm = x;
        xp(k / fiction::kDim, k % fiction::kDim) += eps;
        xm(k / fiction::kDim, k % fiction::kDim) -= eps;
        const double fd = (force.val(xp) - force.val(xm)) / (2.0 * eps);
        maxErr = std::max(maxErr, std::abs(fd - gA[k]));
    }
    const double gradScale = std::max(gA.cwiseAbs().maxCoeff(), 1.0);
    std::cout << std::setprecision(6) << std::fixed
              << "  |grad|_inf=" << gradScale << "  max FD err=" << maxErr << "\n";
    check(maxErr / gradScale < 1e-6, "梯度与有限差分一致（相对误差 < 1e-6）");

    std::vector<fiction::Triplet> trips;
    force.addHess(x, trips);
    check(trips.empty(), "外力海森为空（线性势能）");
}

// -----------------------------------------------------------------------------
// 5.5) NeoHookeanTetEnergy：等参零能量 / 旋转不变性 / FD 梯度&海森 / ν 单调性
// -----------------------------------------------------------------------------
void checkNeoHookean() {
    std::cout << "\n===== 检验 5.5: NeoHookeanTetEnergy =====\n";
    // 参考四面体（单位直角 tet：Dm = I，V0 = 1/6）
    fiction::Positions rest(4, 3);
    rest << 0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1;
    std::vector<std::array<int, 4>> tets = {{0, 1, 2, 3}};
    const double E_young = 2e5, nu = 0.32;

    auto makeEnergy = [&](const fiction::Positions& x0) {
        return std::make_shared<fiction::NeoHookeanTetEnergy>(tets, x0, E_young, nu);
    };

    // (a) 等参位形：零能量/零梯度；Hessian 非零但应对称且 PSD（rest 处 NH 是刚硬的）
    {
        auto nh = makeEnergy(rest);
        const double E0 = nh->val(rest);
        fiction::Positions g = fiction::Positions::Zero(4, 3);
        nh->addGrad(rest, g);
        std::vector<fiction::Triplet> trips;
        nh->addHess(rest, trips);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(12, 12);
        for (const auto& t : trips) H(t.row(), t.col()) += t.value();
        const double hScale = std::max(H.cwiseAbs().maxCoeff(), 1e-30);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
        const double minEig = es.info() == Eigen::Success ? es.eigenvalues().minCoeff() : 1.0;
        std::cout << std::setprecision(3) << std::fixed;
        std::cout << "  等参: E=" << E0 << ", |grad|=" << g.cwiseAbs().maxCoeff()
                  << ", |H|=" << hScale << ", min_eig=" << minEig << "\n";
        check(std::abs(E0) < 1e-10 && g.cwiseAbs().maxCoeff() < 1e-10,
              "等参位形 E=0, grad=0");
        check((H - H.transpose()).cwiseAbs().maxCoeff() / hScale < 1e-10
                  && minEig > -1e-9 * std::max(1.0, hScale),
              "等参 Hessian 对称且 PSD（makePSD 生效）");
    }

    // (b) 旋转客观性：固定原参考网格，仅旋转当前位形，能量仍为 0
    {
        const double ang = 0.7;
        Eigen::Matrix3d R;
        R = Eigen::AngleAxisd(ang, Eigen::Vector3d(1, 2, 3).normalized())
            * Eigen::AngleAxisd(0.21, Eigen::Vector3d(1, 0, 0));
        fiction::Positions rot = rest * R.transpose();
        rot.rowwise() += Eigen::RowVector3d(0.4, -1.2, 2.0); // 平移不应影响能量
        auto nh = makeEnergy(rest); // 参考构形保持不动，不能同步旋转 rest
        const double E = nh->val(rot);
        fiction::Positions g = fiction::Positions::Zero(4, 3);
        nh->addGrad(rot, g);
        std::cout << "  固定参考的旋转位形 E=" << E << ", |grad|=" << g.norm() << "\n";
        check(std::abs(E) < 1e-10 && g.norm() < 1e-9,
              "固定参考网格、仅旋转当前位形时零能量/零弹性力");
    }

    // (c) 三配置 FD：梯度与 Hessian（rest 固定，形变位形评估）
    //     注意：与弹簧检验同理，只有精确 Hessian 已 PSD 的位形才能把「投影后」结果
    //     拿去和 FD 对照；压缩/剪切位形改用 addHessRaw（精确未投影）做 FD。
    {
        auto nh = makeEnergy(rest); // rest 才是参考位形，评估点独立取变形位形
        const std::vector<std::pair<std::string, Eigen::Matrix3d>> cases = {
            {"轴向拉伸 x1.15", (Eigen::Matrix3d)Eigen::DiagonalMatrix<double, 3, 3>(1.15, 1.0, 1.0)},
            {"纯剪切", [] {
                 Eigen::Matrix3d S = Eigen::Matrix3d::Identity();
                 S(0, 1) = 0.3;
                 return S;
             }()},
            {"体积压缩 x0.8", (Eigen::Matrix3d)Eigen::DiagonalMatrix<double, 3, 3>(0.8, 0.8, 0.9)},
        };
        for (const auto& [label, A] : cases) {
            fiction::Positions def = rest * A.transpose(); // F = A
            const int dof = 12;
            fiction::Positions grad = fiction::Positions::Zero(4, 3);
            nh->addGrad(def, grad);
            const fiction::FlatVector gA = fiction::flatten(grad);

            // 精确（未投影）Hessian
            std::vector<fiction::Triplet> trips;
            nh->addHessRaw(def, trips);
            Eigen::MatrixXd H = Eigen::MatrixXd::Zero(dof, dof);
            for (const auto& t : trips) H(t.row(), t.col()) += t.value();
            const double hScale = std::max(H.cwiseAbs().maxCoeff(), 1.0);

            // FD 梯度 + FD 海森（对精确 Hessian）
            const double eps = 1e-7;
            fiction::FlatVector gFd(dof);
            Eigen::MatrixXd Hfd = Eigen::MatrixXd::Zero(dof, dof);
            for (int k = 0; k < dof; ++k) {
                fiction::Positions xp = def, xm = def;
                xp(k / fiction::kDim, k % fiction::kDim) += eps;
                xm(k / fiction::kDim, k % fiction::kDim) -= eps;
                gFd(k) = (nh->val(xp) - nh->val(xm)) / (2.0 * eps);
                fiction::Positions gp = fiction::Positions::Zero(4, 3);
                fiction::Positions gm = fiction::Positions::Zero(4, 3);
                nh->addGrad(xp, gp);
                nh->addGrad(xm, gm);
                Hfd.col(k) = (fiction::flatten(gp) - fiction::flatten(gm)) / (2.0 * eps);
            }
            const double gScale = std::max(gA.cwiseAbs().maxCoeff(), 1.0);
            const double gErr = (gA - gFd).cwiseAbs().maxCoeff();
            const double hErr = (H - Hfd).cwiseAbs().maxCoeff();
            std::cout << std::setprecision(3) << std::fixed;
            std::cout << "  " << label << ": E=" << nh->val(def)
                      << " |grad|=" << gA.cwiseAbs().maxCoeff()
                      << " FD err: grad=" << gErr << " hess=" << hErr << "\n";
            check(gErr / gScale < 1e-6, label + " 梯度 FD 相对误差 < 1e-6");
            check(hErr / hScale < 1e-4, label + " 海森(精确) FD 相对误差 < 1e-4");
            check((H - H.transpose()).cwiseAbs().maxCoeff() / hScale < 1e-10, label + " 海森对称");

            // 投影版：拉伸位形（精确 Hessian 已 PSD）应与精确版一致；
            // 压缩/剪切至少应保持 PSD
            Eigen::MatrixXd Hp = Eigen::MatrixXd::Zero(dof, dof);
            {
                std::vector<fiction::Triplet> ptrips;
                nh->addHess(def, ptrips);
                for (const auto& t : ptrips) Hp(t.row(), t.col()) += t.value();
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Hp);
                const double minEig = es.info() == Eigen::Success ? es.eigenvalues().minCoeff() : 1.0;
                check(minEig > -1e-9 * std::max(1.0, hScale), label + " 投影后海森保持 PSD");
                if (label == "轴向拉伸 x1.15") {
                    check((Hp - H).cwiseAbs().maxCoeff() / hScale < 1e-8,
                          "拉伸位形（PSD 域）投影 == 精确海森");
                }
            }
        }
    }

    // (d) 泊松比单调性：体积压缩位形下，ν 越大能量越大（体积刚度）
    {
        Eigen::Matrix3d A = Eigen::Matrix3d::Identity();
        A(0, 0) = A(1, 1) = 0.85; A(2, 2) = 0.9; // 压缩
        fiction::Positions def = rest * A.transpose();
        double E_lo = -1, E_hi = -1;
        for (double nuC : {0.05, 0.45}) {
            fiction::NeoHookeanTetEnergy nh(tets, rest, E_young, nuC); // 参考 = rest
            const double E = nh.val(def);
            if (nuC < 0.3) E_lo = E; else E_hi = E;
            std::cout << "  压缩位形 ν=" << nuC << " → E=" << E << "\n";
        }
        check(E_hi > E_lo && E_hi > 0.0, "体积压缩能量随 ν 单调增大");
    }
}

// -----------------------------------------------------------------------------
// 6) pin-box + 区域外力 / 重力 仿真
// -----------------------------------------------------------------------------
void checkPinBoxSimulation() {
    std::cout << "\n===== 检验 6: pin-box + force-box / 重力 仿真（bar-6） =====\n";
    const std::string barPath = dataPath("bar-6.msh");

    // 场景 A：左端固定，右端 +x 方向拉伸（关闭重力）
    {
        fiction::IpcSceneOptions opt;
        opt.meshPath = barPath;
        opt.enableGravity = false;
        opt.dt = 0.01;
        opt.tol = 1e-8;
        opt.maxIter = 200;
        opt.pinBoxes = {{-0.1, -0.1, -0.1, 0.6, 0.6, 1.1}};
        opt.forceBoxes = {{{9.4, -0.1, -0.1, 10.1, 0.6, 1.1}, fiction::Vec3(30.0, 0.0, 0.0)}};

        const fiction::IpcScene scene = fiction::buildIpcScene(opt);
        fiction::TimeIntegrator::Config cfg;
        cfg.h = opt.dt;
        cfg.tol = opt.tol;
        cfg.maxIter = opt.maxIter;
        cfg.verbose = false;
        fiction::TimeIntegrator integrator(cfg, scene.model.masses);
        for (auto& e : fiction::makeIpcSceneEnergies(scene, opt)) integrator.addEnergy(std::move(e));

        fiction::Positions x = scene.model.x0;
        fiction::Velocities v = fiction::Velocities::Zero(scene.model.numNodes(), fiction::kDim);
        for (int s = 0; s < 50; ++s) {
            const auto state = integrator.stepForward(x, v, scene.dbc);
            x = state.x;
            v = state.v;
        }

        // 固定端位移应为 0（DBC 自由度消去）
        double maxPinDisp = 0.0;
        for (int i : scene.pinnedNodes) maxPinDisp = std::max(maxPinDisp, (x.row(i) - scene.model.x0.row(i)).norm());
        // 自由端（x 接近 10）应被拉向 +x
        double maxTipDx = -1e9;
        for (int i = 0; i < x.rows(); ++i) {
            if (scene.model.x0(i, 0) > 9.9)
                maxTipDx = std::max(maxTipDx, x(i, 0) - scene.model.x0(i, 0));
        }
        std::cout << std::setprecision(6) << std::fixed;
        std::cout << "  A: pin 节点最大位移=" << maxPinDisp << ", 自由端 max dx=" << maxTipDx
                  << ", max|v|=" << v.cwiseAbs().maxCoeff() << "\n";
        check(maxPinDisp < 1e-12, "pin-box 节点保持固定（位移 = 0）");
        check(maxTipDx > 0.0, "外力方向 +x 拉伸自由端");
        check(x.allFinite() && v.allFinite(), "结果有限（无 NaN/Inf）");
    }

    // 场景 B：左端固定，重力 -z 下垂
    {
        fiction::IpcSceneOptions opt;
        opt.meshPath = barPath;
        opt.enableGravity = true;
        opt.gravity = fiction::Vec3(0.0, 0.0, -9.81);
        opt.dt = 0.01;
        opt.tol = 1e-8;
        opt.maxIter = 200;
        opt.pinBoxes = {{-0.1, -0.1, -0.1, 0.6, 0.6, 1.1}};

        const fiction::IpcScene scene = fiction::buildIpcScene(opt);
        fiction::TimeIntegrator::Config cfg;
        cfg.h = opt.dt;
        cfg.tol = opt.tol;
        cfg.maxIter = opt.maxIter;
        cfg.verbose = false;
        fiction::TimeIntegrator integrator(cfg, scene.model.masses);
        for (auto& e : fiction::makeIpcSceneEnergies(scene, opt)) integrator.addEnergy(std::move(e));

        fiction::Positions x = scene.model.x0;
        fiction::Velocities v = fiction::Velocities::Zero(scene.model.numNodes(), fiction::kDim);
        for (int s = 0; s < 50; ++s) {
            const auto state = integrator.stepForward(x, v, scene.dbc);
            x = state.x;
            v = state.v;
        }
        // 自由端（x 接近 10）应下垂（z 减小）
        double maxTipDz = 1e9;
        for (int i = 0; i < x.rows(); ++i) {
            if (scene.model.x0(i, 0) > 9.9)
                maxTipDz = std::min(maxTipDz, x(i, 2) - scene.model.x0(i, 2));
        }
        std::cout << "  B: 自由端 min dz=" << maxTipDz << ", max|v|=" << v.cwiseAbs().maxCoeff() << "\n";
        check(maxTipDz < 0.0, "重力 -z 下垂自由端");
        check(x.allFinite() && v.allFinite(), "结果有限（无 NaN/Inf）");
    }
}

} // namespace

int main() {
    std::cout << "fiction IPC .msh 加载 + mass-spring 数值验证 (纯 3D)\n";
    checkLoaderBasics();
    checkLoaderFallbackAndErrors();
    checkFreeFall();
    checkSpringFiniteDifference();
    checkRegionForceFiniteDifference();
    checkNeoHookean();
    checkPinBoxSimulation();
    std::cout << "\n===== mesh_check " << (g_failCount == 0 ? "全部 PASS" : "存在 FAIL")
              << " =====\n";
    return g_failCount == 0 ? 0 : 1;
}
