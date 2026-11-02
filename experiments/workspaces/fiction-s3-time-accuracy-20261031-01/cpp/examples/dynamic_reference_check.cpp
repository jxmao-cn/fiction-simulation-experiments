// T0：用测试专用二次能量核对隐式 Euler 的离散动态参考。
#include "time_integrator.h"
#include "energies/Gravity.h"

#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace {
constexpr double kStiffness = 20.0;
constexpr double kSolveTol = 1e-10;
constexpr double kDiscreteTol = 1e-10;

class QuadraticEnergy final : public fiction::IEnergy {
public:
    explicit QuadraticEnergy(double k) : k_(k) {}
    const char* name() const override { return "TestQuadratic"; }
    double val(const fiction::Positions& x) const override { return 0.5 * k_ * x.squaredNorm(); }
    void addGrad(const fiction::Positions& x, fiction::Positions& grad) const override { grad += k_ * x; }
    void addHess(const fiction::Positions& x, std::vector<fiction::Triplet>& trips) const override {
        (void)x;
        for (int i = 0; i < 3 * x.rows(); ++i) trips.emplace_back(i, i, k_);
    }
private:
    double k_;
};

bool finiteState(const fiction::TimeIntegrator::State& s) {
    return s.x.allFinite() && s.v.allFinite();
}

void writeRow(std::ofstream* out, const std::string& scenario, double mass, double h,
              int step, int component, double x, double xr, double v, double vr) {
    if (!out || !*out) return;
    *out << scenario << ',' << mass << ',' << h << ',' << step << ',' << step*h << ','
         << component << ',' << x << ',' << xr << ',' << (x-xr) << ','
         << v << ',' << vr << ',' << (v-vr) << '\n';
}

bool oscillatorDiscrete(double mass, double h, int steps, std::ofstream* csv) {
    fiction::Positions x(1,3), v(1,3);
    x << 0.7, -0.25, 0.4;
    v << -0.3, 0.8, 0.15;
    Eigen::VectorXd masses(1); masses << mass;
    fiction::TimeIntegrator::Config cfg;
    cfg.h = h; cfg.tol = kSolveTol; cfg.maxIter = 100; cfg.verbose = false;
    fiction::TimeIntegrator integrator(cfg, masses);
    integrator.addEnergy(std::make_shared<QuadraticEnergy>(kStiffness));
    fiction::DBCInfo dbc;
    fiction::Positions xr = x;
    fiction::Velocities vr = v;
    double maxXErr = 0.0, maxVErr = 0.0;
    for (int n = 1; n <= steps; ++n) {
        const double denom = 1.0 + h*h*kStiffness/mass;
        vr = (vr - h*(kStiffness/mass)*xr) / denom;
        xr = xr + h*vr;
        auto s = integrator.stepForward(x, v, dbc);
        if (!finiteState(s)) return false;
        for (int d = 0; d < 3; ++d) {
            maxXErr = std::max(maxXErr, std::abs(s.x(0,d)-xr(0,d)));
            maxVErr = std::max(maxVErr, std::abs(s.v(0,d)-vr(0,d)));
            writeRow(csv, "oscillator_discrete_cumulative_reference", mass, h, n, d,
                     s.x(0,d), xr(0,d), s.v(0,d), vr(0,d));
        }
        x = s.x; v = s.v;
    }
    const bool pass = maxXErr <= kDiscreteTol && maxVErr <= kDiscreteTol;
    std::cout << "离散振子 m=" << mass << " h=" << h << " steps=" << steps
              << " max_x_error=" << maxXErr << " max_v_error=" << maxVErr
              << " tolerance=" << kDiscreteTol
              << (pass ? " PASS\n" : " FAIL\n");
    return pass;
}

struct ContinuousResult { double error; double energyRatio; };
ContinuousResult oscillatorContinuous(double mass, double h, int steps, std::ofstream* csv) {
    fiction::Positions x(1,3), v(1,3);
    x << 0.7, -0.25, 0.4;
    v << -0.3, 0.8, 0.15;
    const fiction::Positions x0 = x;
    const fiction::Velocities v0 = v;
    Eigen::VectorXd masses(1); masses << mass;
    fiction::TimeIntegrator::Config cfg;
    cfg.h = h; cfg.tol = kSolveTol; cfg.maxIter = 100; cfg.verbose = false;
    fiction::TimeIntegrator integrator(cfg, masses);
    integrator.addEnergy(std::make_shared<QuadraticEnergy>(kStiffness));
    fiction::DBCInfo dbc;
    double maxXError = 0.0, maxVError = 0.0;
    const double omega = std::sqrt(kStiffness/mass);
    for (int n = 1; n <= steps; ++n) {
        auto s = integrator.stepForward(x, v, dbc);
        if (!finiteState(s)) return {INFINITY, INFINITY};
        const double t = h*n;
        for (int d = 0; d < 3; ++d) {
            const double xr = x0(0,d)*std::cos(omega*t) + v0(0,d)/omega*std::sin(omega*t);
            const double vr = -omega*x0(0,d)*std::sin(omega*t) + v0(0,d)*std::cos(omega*t);
            maxXError = std::max(maxXError, std::abs(s.x(0,d)-xr));
            maxVError = std::max(maxVError, std::abs(s.v(0,d)-vr));
            writeRow(csv, "oscillator_continuous", mass, h, n, d, s.x(0,d), xr, s.v(0,d), vr);
        }
        x = s.x; v = s.v;
    }
    const double e0 = 0.5*mass*v0.squaredNorm() + 0.5*kStiffness*x0.squaredNorm();
    const double ef = 0.5*mass*v.squaredNorm() + 0.5*kStiffness*x.squaredNorm();
    // 用无量纲状态误差比较不同物理量纲：x/L0 与 v/(omega*L0)。
    const double L0 = std::max(x0.cwiseAbs().maxCoeff(), v0.cwiseAbs().maxCoeff()/omega);
    if (!(std::isfinite(e0) && std::isfinite(ef) && e0 > 0.0 && ef > 0.0))
        return {INFINITY, INFINITY};
    return {std::max(maxXError/L0, maxVError/(omega*L0)), ef/e0};
}

bool simpleExactCases(double h, int steps, std::ofstream* csv) {
    bool ok = true;
    for (const std::string scenario : {"inertia", "gravity"}) {
        const double mass = 2.0;
        fiction::Positions x(1,3), v(1,3);
        x << 0.2, -0.4, 1.1; v << 0.3, -0.2, 0.5;
        const fiction::Positions x0=x;
        const fiction::Velocities v0=v;
        Eigen::VectorXd masses(1); masses << mass;
        fiction::TimeIntegrator::Config cfg;
        cfg.h=h; cfg.tol=kSolveTol; cfg.maxIter=100; cfg.verbose=false;
        fiction::TimeIntegrator integrator(cfg,masses);
        const fiction::Vec3 g(0.0,0.0,-9.81);
        if (scenario == "gravity") integrator.addEnergy(std::make_shared<fiction::ExampleGravityEnergy>(masses,g));
        fiction::DBCInfo dbc;
        for (int n=1; n<=steps; ++n) {
            fiction::Velocities vr = v0;
            if (scenario == "gravity") vr.row(0) += (n*h*g).transpose();
            fiction::Positions xr = x0 + (n*h)*v0;
            if (scenario == "gravity") xr.row(0) += (0.5*n*(n+1)*h*h*g).transpose();
            auto s=integrator.stepForward(x,v,dbc);
            if (!finiteState(s)) return false;
            for (int d=0; d<3; ++d) {
                const double err=std::max(std::abs(s.x(0,d)-xr(0,d)),std::abs(s.v(0,d)-vr(0,d)));
                ok = ok && err <= kDiscreteTol;
                writeRow(csv,scenario,mass,h,n,d,s.x(0,d),xr(0,d),s.v(0,d),vr(0,d));
            }
            x=s.x; v=s.v;
        }
    }
    std::cout << "纯惯性与恒定重力离散更新 tolerance=" << kDiscreteTol << (ok?" PASS\n":" FAIL\n");
    return ok;
}
} // namespace

int main(int argc, char** argv) {
    std::string csvPath;
    if (argc == 3 && std::string(argv[1]) == "--csv") csvPath = argv[2];
    else if (argc != 1) {
        std::cerr << "用法: dynamic_reference_check [--csv output.csv]\n";
        return 2;
    }
    std::ofstream csv;
    std::ofstream* csvPtr = nullptr;
    if (!csvPath.empty()) {
        csv.open(csvPath);
        if (!csv) { std::cerr << "无法写入 CSV: " << csvPath << '\n'; return 2; }
        csv << "scenario,mass,h,step,time,component,x_numeric,x_reference,x_error,v_numeric,v_reference,v_error\n";
        csv << std::setprecision(17);
        csvPtr = &csv;
    }
    std::cout << std::setprecision(12) << "T0 config: stiffness=" << kStiffness
              << " solve_tol=" << kSolveTol << " discrete_state_tol=" << kDiscreteTol
              << "\n";
    bool ok = simpleExactCases(0.02, 120, csvPtr);
    for (double m : {1.0, 2.0, 5.0}) {
        for (double h : {0.01, 0.02, 0.04}) {
            ok = oscillatorDiscrete(m,h,120,csvPtr) && ok;
        }
    }
    const double m=2.0, finalTime=2.0;
    ContinuousResult coarse=oscillatorContinuous(m,0.04,50,csvPtr);
    ContinuousResult medium=oscillatorContinuous(m,0.02,100,csvPtr);
    ContinuousResult fine=oscillatorContinuous(m,0.01,200,csvPtr);
    const bool trend = std::isfinite(coarse.error) && coarse.error > medium.error && medium.error > fine.error;
    const bool damping = std::isfinite(coarse.energyRatio) && std::isfinite(medium.energyRatio) &&
                         std::isfinite(fine.energyRatio) && coarse.energyRatio > 0.0 &&
                         medium.energyRatio > 0.0 && fine.energyRatio > 0.0 &&
                         coarse.energyRatio < medium.energyRatio && medium.energyRatio < fine.energyRatio &&
                         fine.energyRatio < 1.0;
    std::cout << "连续谐振子 m=" << m << " t=" << finalTime
              << " error(h=.04,.02,.01)=" << coarse.error << ',' << medium.error << ',' << fine.error
              << " convergence=" << (trend?"PASS":"FAIL")
              << " energy_ratio(h=.04,.02,.01)=" << coarse.energyRatio << ',' << medium.energyRatio << ',' << fine.energyRatio
              << " implicit_euler_damping=" << (damping?"PASS":"FAIL") << '\n';
    ok = ok && trend && damping;
    if (csvPtr) { csv.flush(); if (!csv) { std::cerr << "CSV 写入失败\n"; return 2; } }
    return ok ? 0 : 1;
}
