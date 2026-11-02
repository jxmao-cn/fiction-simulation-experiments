// S1 前置检查：静态平衡与逐步求解报告语义。
#include "time_integrator.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

namespace {
class QuadraticEnergy final : public fiction::IEnergy {
public:
    explicit QuadraticEnergy(double stiffness) : k_(stiffness) {}
    const char* name() const override { return "quadratic-check"; }
    double hessianWeight(double h) const override { return h * h; }
    double val(const fiction::Positions& x) const override { return 0.5 * k_ * x.squaredNorm(); }
    void addGrad(const fiction::Positions& x, fiction::Positions& g) const override { g += k_ * x; }
    void addHess(const fiction::Positions& x, std::vector<fiction::Triplet>& t) const override {
        (void)x;
        for (int i = 0; i < 3; ++i) t.emplace_back(i, i, k_);
    }
private:
    double k_;
};

class NonfiniteGradient final : public fiction::IEnergy {
public:
    const char* name() const override { return "nonfinite-check"; }
    double hessianWeight(double) const override { return 1.0; }
    double val(const fiction::Positions&) const override { return 0.0; }
    void addGrad(const fiction::Positions&, fiction::Positions& g) const override {
        g(0, 0) = std::numeric_limits<double>::quiet_NaN();
    }
    void addHess(const fiction::Positions&, std::vector<fiction::Triplet>& t) const override {
        for (int i = 0; i < 3; ++i) t.emplace_back(i, i, 1.0);
    }
};

bool near(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol; }
}

int main() {
    using namespace fiction;
    Eigen::VectorXd mass(1); mass << 2.0;
    Positions x(1, 3); x << 1.0, -2.0, 0.5;
    Velocities v(1, 3); v << 0.4, 0.2, -0.1;
    const DBCInfo noDbc;

    auto makeStatic = [&](double h) {
        TimeIntegrator::Config c; c.staticSolve = true; c.h = h; c.tol = 1e-10;
        c.maxIter = 10; c.verbose = false;
        TimeIntegrator solver(c, mass); solver.addEnergy(std::make_shared<QuadraticEnergy>(4.0));
        return solver.stepForward(x, v, noDbc);
    };
    const auto sA = makeStatic(0.001), sB = makeStatic(7.0);
    if (!sA.report.converged || !sB.report.converged || !sA.x.isZero(1e-10) ||
        !sB.x.isZero(1e-10) || !sA.v.isZero(0.0) || !sB.v.isZero(0.0) ||
        !near(sA.report.searchDirectionRate, sB.report.searchDirectionRate) ||
        sA.report.solver.empty()) {
        std::cerr << "FAIL static equilibrium / h independence / solver report\n"; return 1;
    }

    // 默认动态路径与单自由度线性势的隐式欧拉离散更新对照。
    TimeIntegrator::Config dc; dc.h = 0.02; dc.tol = 1e-12; dc.maxIter = 10; dc.verbose = false;
    TimeIntegrator dynamic(dc, mass); dynamic.addEnergy(std::make_shared<QuadraticEnergy>(4.0));
    const auto d = dynamic.stepForward(x, v, noDbc);
    const double den = 1.0 + dc.h * dc.h * 4.0 / mass[0];
    for (int j = 0; j < 3; ++j) {
        const double vRef = (v(0,j) - dc.h * (4.0 / mass[0]) * x(0,j)) / den;
        const double xRef = x(0,j) + dc.h * vRef;
        if (!near(d.x(0,j), xRef, 1e-10) || !near(d.v(0,j), vRef, 1e-10)) {
            std::cerr << "FAIL dynamic implicit-Euler reference\n"; return 1;
        }
    }
    if (!d.report.converged) { std::cerr << "FAIL dynamic convergence report\n"; return 1; }

    TimeIntegrator::Config capped = dc; capped.staticSolve = true; capped.maxIter = 0;
    TimeIntegrator limited(capped, mass); limited.addEnergy(std::make_shared<QuadraticEnergy>(4.0));
    const auto m = limited.stepForward(x, v, noDbc);
    if (m.report.converged || !m.report.maxIterHit || m.report.iterations != 0) {
        std::cerr << "FAIL maxIter=0 non-equilibrium report\n"; return 1;
    }

    TimeIntegrator::Config nc = capped; nc.maxIter = 2;
    TimeIntegrator broken(nc, mass); broken.addEnergy(std::make_shared<NonfiniteGradient>());
    const auto n = broken.stepForward(x, v, noDbc);
    if (n.report.converged || n.report.solver != "failed") {
        std::cerr << "FAIL nonfinite linear solve report\n"; return 1;
    }
    if (!std::isfinite(sA.report.freeGradientInf) || !std::isfinite(sA.report.dbcMaxError) ||
        sA.report.dbcTotal != 0 || sA.report.dbcSatisfied != 0 || sA.report.iterations < 1) {
        std::cerr << "FAIL report diagnostics\n"; return 1;
    }
    std::cout << "PASS static equilibrium; dynamic reference; failure states; StepReport\n";
    return 0;
}
