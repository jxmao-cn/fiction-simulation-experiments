// =============================================================================
// utils.cpp — utils.h 的非模板实现（纯 3D 版）
// =============================================================================

#include "utils.h"
#include <Eigen/Eigenvalues>

namespace fiction {

Eigen::Matrix3d makePSD(const Eigen::Matrix3d& H) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(H);
    if (solver.info() != Eigen::Success) {
        return 0.5 * (H + H.transpose());
    }
    Eigen::Vector3d lam = solver.eigenvalues();
    Eigen::Matrix3d V = solver.eigenvectors();
    for (int i = 0; i < 3; ++i) if (lam(i) < 0.0) lam(i) = 0.0;
    return V * lam.asDiagonal() * V.transpose();
}

Eigen::Matrix<double, 6, 6> makePSD(const Eigen::Matrix<double, 6, 6>& H) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> solver(H);
    if (solver.info() != Eigen::Success) {
        return 0.5 * (H + H.transpose());
    }
    Eigen::Matrix<double, 6, 1> lam = solver.eigenvalues();
    Eigen::Matrix<double, 6, 6> V = solver.eigenvectors();
    for (int i = 0; i < 6; ++i) if (lam(i) < 0.0) lam(i) = 0.0;
    return V * lam.asDiagonal() * V.transpose();
}

Eigen::Matrix<double, 9, 9> makePSD(const Eigen::Matrix<double, 9, 9>& H) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 9, 9>> solver(H);
    if (solver.info() != Eigen::Success) {
        return 0.5 * (H + H.transpose());
    }
    Eigen::Matrix<double, 9, 1> lam = solver.eigenvalues();
    Eigen::Matrix<double, 9, 9> V = solver.eigenvectors();
    for (int i = 0; i < 9; ++i) if (lam(i) < 0.0) lam(i) = 0.0;
    return V * lam.asDiagonal() * V.transpose();
}

Eigen::Matrix<double, 12, 12> makePSD(const Eigen::Matrix<double, 12, 12>& H) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 12, 12>> solver(H);
    if (solver.info() != Eigen::Success) {
        return 0.5 * (H + H.transpose());
    }
    Eigen::Matrix<double, 12, 1> lam = solver.eigenvalues();
    Eigen::Matrix<double, 12, 12> V = solver.eigenvectors();
    for (int i = 0; i < 12; ++i) if (lam(i) < 0.0) lam(i) = 0.0;
    return V * lam.asDiagonal() * V.transpose();
}

// -----------------------------------------------------------------------------
// 范数辅助
// -----------------------------------------------------------------------------
double infNorm(const Eigen::VectorXd& v) {
    return v.lpNorm<Eigen::Infinity>();
}

double infNormPerH(const Eigen::VectorXd& p, double h) {
    if (h == 0.0) return p.lpNorm<Eigen::Infinity>();
    return (p / h).lpNorm<Eigen::Infinity>();
}

}  // namespace fiction
