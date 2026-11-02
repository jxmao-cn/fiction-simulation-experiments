#pragma once
// =============================================================================
// camera.h — 视图/投影矩阵数学（GL 固定管线与 CPU 软光栅化共用）
// =============================================================================
// 说明：
//   - 列主序矩阵（可直接 glLoadMatrixd / 与 Eigen 列主序乘向量）
//   - 与 viewer.h 的轨道相机约定一致：世界系 z 轴向上，yaw/pitch/radius
//   - 无任何图形 API 依赖，纯 Eigen，可供无显示环境（SSH）复用
// =============================================================================

#include "types.h"

#include <Eigen/Dense>

#include <cmath>

namespace fiction {

constexpr double kDeg2Rad = 3.14159265358979323846 / 180.0;

/// 标准透视投影矩阵（列主序，等价 gluPerspective）
inline Eigen::Matrix4d perspective(double fovDeg, double aspect,
                                   double zNear, double zFar) {
    const double f = 1.0 / std::tan(fovDeg * kDeg2Rad * 0.5);
    Eigen::Matrix4d P = Eigen::Matrix4d::Zero();
    P(0, 0) = f / aspect;
    P(1, 1) = f;
    P(2, 2) = (zFar + zNear) / (zNear - zFar);
    P(3, 2) = -1.0;
    P(2, 3) = 2.0 * zFar * zNear / (zNear - zFar);
    return P;
}

/// lookAt 视图矩阵（列主序，等价 gluLookAt）
inline Eigen::Matrix4d lookAt(const Vec3& eye, const Vec3& target, const Vec3& up) {
    const Vec3 forward = (target - eye).normalized();
    const Vec3 s = forward.cross(up).normalized();
    const Vec3 u = s.cross(forward);
    Eigen::Matrix4d V = Eigen::Matrix4d::Identity();
    V(0, 0) = s(0); V(0, 1) = s(1); V(0, 2) = s(2); V(0, 3) = -s.dot(eye);
    V(1, 0) = u(0); V(1, 1) = u(1); V(1, 2) = u(2); V(1, 3) = -u.dot(eye);
    V(2, 0) = -forward(0); V(2, 1) = -forward(1); V(2, 2) = -forward(2);
    V(2, 3) = forward.dot(eye);
    return V;
}

} // namespace fiction
