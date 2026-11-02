// 回归：DBC 刚度翻倍后必须用新罚刚度重算 Newton 方向。
#include "time_integrator.h"

#include <cmath>
#include <iostream>
#include <memory>

namespace {
class TinyQuadratic final : public fiction::IEnergy {
public:
    const char* name() const override { return "tiny-quadratic-check"; }
    double hessianWeight(double h) const override { return h * h; }
    double val(const fiction::Positions& x) const override { return 0.5 * x.squaredNorm(); }
    void addGrad(const fiction::Positions& x, fiction::Positions& grad) const override { grad += x; }
    void addHess(const fiction::Positions&, std::vector<fiction::Triplet>& trips) const override {
        for (int d = 0; d < 3; ++d) trips.emplace_back(d, d, 1.0);
    }
};
}

int main() {
    using namespace fiction;
    Eigen::VectorXd mass(1); mass << 1.0;
    Positions x(1, 3); x << 3e-10, 0.0, 0.0;
    Velocities v = Velocities::Zero(1, 3);
    DBCInfo dbc;
    dbc.stiffness = 1.0;
    const Vec3 zero = Vec3::Zero();
    dbc.addDBC(0, zero, zero);
    dbc.targetPositions.emplace_back(5e-10, 0.0, 0.0);

    TimeIntegrator::Config cfg;
    cfg.h = 1.0;
    cfg.tol = 1e-10;
    cfg.maxIter = 1;
    cfg.staticSolve = true;
    cfg.verbose = false;
    TimeIntegrator solver(cfg, mass);
    solver.addEnergy(std::make_shared<TinyQuadratic>());

    const double initialX = x(0, 0);
    // Compare both states under the refreshed objective (k=2), not across two objectives.
    const double updatedObjectiveInitialEnergy = 0.5 * initialX * initialX +
        0.5 * 2.0 * std::pow(initialX - 5e-10, 2);
    const auto state = solver.stepForward(x, v, dbc);
    const double finalX = state.x(0, 0);
    const double finalEnergy = 0.5 * finalX * finalX +
        0.5 * state.report.dbcStiffness * std::pow(finalX - 5e-10, 2);

    if (!state.report.maxIterHit || state.report.converged || state.report.iterations != 1 ||
        state.report.dbcStiffness != 2.0 || !(finalX > initialX) ||
        !(finalEnergy < updatedObjectiveInitialEnergy)) {
        std::cerr << "FAIL DBC refresh counterexample: x0=" << initialX << " x1=" << finalX
                  << " k=" << state.report.dbcStiffness << " E0=" << updatedObjectiveInitialEnergy
                  << " E1=" << finalEnergy << " iterations=" << state.report.iterations
                  << " maxIterHit=" << state.report.maxIterHit
                  << " converged=" << state.report.converged << '\n';
        return 1;
    }
    std::cout << "PASS DBC penalty refresh reverses old direction and lowers updated objective; "
                 "maxIter report is preserved\n";
    return 0;
}
