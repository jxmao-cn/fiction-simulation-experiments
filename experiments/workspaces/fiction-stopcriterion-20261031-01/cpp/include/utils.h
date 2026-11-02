#pragma once
// =============================================================================
// utils.h — 数值工具集（纯 3D 版，对应 Python_slovement/utils.py）
// =============================================================================
// 变更：
//   - 移除 2D 特化 2x2/4x4，仅保留 3D 相关的 3x3/6x6/9x9/12x12 + 模板/动态版
//   - 2D 算例经 upgrade2DTo3D 补零后，边单元局部块为 6x6，非 4x4
// =============================================================================

#include <Eigen/Dense>

namespace fiction {

// -----------------------------------------------------------------------------
// Hessian 半正定投影 (PSD projection) — 纯 3D
// -----------------------------------------------------------------------------

/// 对 3x3 对称矩阵做 PSD 投影（单节点 3D 块）
Eigen::Matrix3d makePSD(const Eigen::Matrix3d& H);

/// 对 6x6 对称矩阵做 PSD 投影（边单元：2 节点 ×3 维）
Eigen::Matrix<double, 6, 6> makePSD(const Eigen::Matrix<double, 6, 6>& H);

/// 对 9x9 对称矩阵做 PSD 投影（三角面单元：3 节点 ×3 维）
Eigen::Matrix<double, 9, 9> makePSD(const Eigen::Matrix<double, 9, 9>& H);

/// 对 12x12 对称矩阵做 PSD 投影（四面体单元：4 节点 ×3 维）
Eigen::Matrix<double, 12, 12> makePSD(const Eigen::Matrix<double, 12, 12>& H);

// 遗留兼容：2x2 / 4x4 已废弃，保留内联转发以便历史能量文件仍能编译（内部升维到 3D 再截回）
inline Eigen::Matrix2d makePSD(const Eigen::Matrix2d& H) {
    Eigen::Matrix3d H3 = Eigen::Matrix3d::Zero();
    H3.topLeftCorner<2,2>() = H;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(0.5*(H3+H3.transpose()));
    if (solver.info() != Eigen::Success) return 0.5*(H+H.transpose());
    Eigen::Vector3d lam = solver.eigenvalues();
    Eigen::Matrix3d V = solver.eigenvectors();
    for (int i=0;i<3;++i) if(lam(i)<0) lam(i)=0;
    Eigen::Matrix3d R = V*lam.asDiagonal()*V.transpose();
    return R.topLeftCorner<2,2>();
}
inline Eigen::Matrix4d makePSD(const Eigen::Matrix4d& H) {
    Eigen::Matrix<double,6,6> H6 = Eigen::Matrix<double,6,6>::Zero();
    // 2D 边 4x4 映射到 3D 边 6x6 的对角（x,y 保留，z 零）
    for(int nI=0;nI<2;++nI) for(int nJ=0;nJ<2;++nJ) for(int r=0;r<2;++r) for(int c=0;c<2;++c)
        H6(nI*3+r, nJ*3+c) = H(nI*2+r, nJ*2+c);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6>> solver(0.5*(H6+H6.transpose()));
    if (solver.info() != Eigen::Success) return 0.5*(H+H.transpose());
    Eigen::Matrix<double,6,1> lam = solver.eigenvalues();
    Eigen::Matrix<double,6,6> V = solver.eigenvectors();
    for(int i=0;i<6;++i) if(lam(i)<0) lam(i)=0;
    Eigen::Matrix<double,6,6> R = V*lam.asDiagonal()*V.transpose();
    Eigen::Matrix4d Rout;
    for(int nI=0;nI<2;++nI) for(int nJ=0;nJ<2;++nJ) for(int r=0;r<2;++r) for(int c=0;c<2;++c)
        Rout(nI*2+r, nJ*2+c) = R(nI*3+r, nJ*3+c);
    return Rout;
}

/// 通用 N x N 对称矩阵 PSD 投影模板
template <int N>
Eigen::Matrix<double, N, N> makePSD(const Eigen::Matrix<double, N, N>& H) {
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, N, N>> solver(H);
    if (solver.info() != Eigen::Success) {
        return 0.5 * (H + H.transpose());
    }
    Eigen::Matrix<double, N, 1> lam = solver.eigenvalues();
    Eigen::Matrix<double, N, N> V = solver.eigenvectors();
    for (int i = 0; i < N; ++i) {
        if (lam(i) < 0.0) lam(i) = 0.0;
    }
    return V * lam.asDiagonal() * V.transpose();
}

/// 动态尺寸 PSD 投影（运行时 N，3D 下 N%3==0）
inline Eigen::MatrixXd makePSD_dynamic(const Eigen::MatrixXd& H) {
    if (H.rows() != H.cols()) return 0.5 * (H + H.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
    if (solver.info() != Eigen::Success) {
        return 0.5 * (H + H.transpose());
    }
    Eigen::VectorXd lam = solver.eigenvalues();
    Eigen::MatrixXd V = solver.eigenvectors();
    for (int i = 0; i < lam.size(); ++i) if (lam(i) < 0) lam(i) = 0;
    return V * lam.asDiagonal() * V.transpose();
}

// -----------------------------------------------------------------------------
// 范数辅助
// -----------------------------------------------------------------------------

double infNorm(const Eigen::VectorXd& v);
double infNormPerH(const Eigen::VectorXd& p, double h);

}  // namespace fiction
