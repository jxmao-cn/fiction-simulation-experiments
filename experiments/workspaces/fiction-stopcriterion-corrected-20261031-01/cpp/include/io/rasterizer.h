#pragma once
// =============================================================================
// io/rasterizer.h — CPU 软件光栅化（无 OpenGL 依赖，SSH/无显示环境可用）
// =============================================================================
// 用途：把 Positions + 表面三角面渲染成图像缓冲，供导出视频帧。
//   - z-buffer 深度测试；逐面定向光明暗（|dot| 双面可见）
//   - 超采样（默认 2x）后盒式降采样抗锯齿
//   - 支持三角面、地面网格线、节点高亮
// 坐标/相机约定与 viewer 一致：世界系 z 向上，NDC 深度越小越近。
// =============================================================================

#include "camera.h"
#include "types.h"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace fiction {

/// 轨道相机参数（yaw/pitch/radius 约定与 viewer 相同）
struct SoftCamera {
    Vec3 target = Vec3(0.0, 0.0, 0.0);
    double radius = 5.0;
    double yawDeg = 45.0;
    double pitchDeg = 25.0;
    double fovDeg = 45.0;

    Vec3 eye() const {
        const double cy = std::cos(yawDeg * kDeg2Rad), sy = std::sin(yawDeg * kDeg2Rad);
        const double cp = std::cos(pitchDeg * kDeg2Rad), sp = std::sin(pitchDeg * kDeg2Rad);
        return target + radius * Vec3(cp * cy, cp * sy, sp);
    }
};

class SoftRenderer {
public:
    SoftRenderer(int outWidth, int outHeight, int supersample = 2)
        : outW_(std::max(1, outWidth)),
          outH_(std::max(1, outHeight)),
          ss_(std::max(1, supersample)),
          w_(outW_ * ss_),
          h_(outH_ * ss_),
          depth_(static_cast<size_t>(w_) * static_cast<size_t>(h_),
                 std::numeric_limits<double>::max()),
          color_(static_cast<size_t>(w_) * static_cast<size_t>(h_) * 3, 0) {}

    int width() const { return outW_; }
    int height() const { return outH_; }

    void clear(float r, float g, float b) {
        std::fill(depth_.begin(), depth_.end(), std::numeric_limits<double>::max());
        const unsigned char cr = toByte(r), cg = toByte(g), cb = toByte(b);
        for (size_t i = 0; i < color_.size(); i += 3) {
            color_[i] = cr;
            color_[i + 1] = cg;
            color_[i + 2] = cb;
        }
    }

    void setCamera(const SoftCamera& cam) {
        const Vec3 eye = cam.eye();
        view_ = lookAt(eye, cam.target, Vec3(0.0, 0.0, 1.0));
        proj_ = perspective(cam.fovDeg, static_cast<double>(w_) / h_, kZNear, kZFar);
    }

    /// 三角面（逐面法线定向光；面积为零的退化面跳过）
    void drawTriangles(const Positions& x,
                       const std::vector<std::array<int, 3>>& faces,
                       float r, float g, float b) {
        if (x.rows() == 0 || faces.empty()) return;
        const Vec3 light = Vec3(0.4, -0.6, 0.8).normalized();
        const Eigen::Matrix4d mvp = proj_ * view_;
        for (const auto& f : faces) {
            if (f[0] < 0 || f[0] >= x.rows() || f[1] < 0 || f[1] >= x.rows() ||
                f[2] < 0 || f[2] >= x.rows())
                continue;
            const Vec3 p0 = x.row(f[0]).transpose();
            const Vec3 p1 = x.row(f[1]).transpose();
            const Vec3 p2 = x.row(f[2]).transpose();
            Vec3 n = (p1 - p0).cross(p2 - p0);
            const double area2 = n.norm();
            if (area2 < 1e-24) continue;
            n /= area2;
            const float shade = static_cast<float>(0.35 + 0.65 * std::abs(n.dot(light)));

            Eigen::Vector3d s0, s1, s2;
            if (!project(mvp, p0, s0) || !project(mvp, p1, s1) || !project(mvp, p2, s2))
                continue;
            rasterizeTriangle(s0, s1, s2, shade, r, g, b);
        }
    }

    /// 地面网格（origin/normal 与 viewer::drawGround 相同）
    void drawGrid(const Vec3& origin, const Vec3& normal, double extent, double spacing,
                  float r, float g, float b) {
        if (spacing <= 0.0 || extent <= 0.0) return;
        Vec3 n = normal;
        if (n.squaredNorm() < 1e-24) n = Vec3(0.0, 0.0, 1.0);
        n.normalize();
        const Vec3 ref = std::abs(n(2)) < 0.9 ? Vec3(0.0, 0.0, 1.0) : Vec3(1.0, 0.0, 0.0);
        const Vec3 t1 = ref.cross(n).normalized();
        const Vec3 t2 = n.cross(t1).normalized();

        const Eigen::Matrix4d mvp = proj_ * view_;
        auto line = [&](const Vec3& a, const Vec3& c) {
            Eigen::Vector3d sa, sc;
            if (!project(mvp, a, sa) || !project(mvp, c, sc)) return;
            drawLine(sa, sc, r, g, b);
        };
        for (double s = -extent; s <= extent + 1e-12; s += spacing) {
            line(origin + t1 * s - t2 * extent, origin + t1 * s + t2 * extent);
            line(origin + t2 * s - t1 * extent, origin + t2 * s + t1 * extent);
        }
    }

    /// 节点高亮（如 pin/外力节点）
    void drawNodes(const Positions& x, const std::vector<int>& nodes,
                   int radiusPx, float r, float g, float b) {
        if (x.rows() == 0 || nodes.empty()) return;
        const Eigen::Matrix4d mvp = proj_ * view_;
        const unsigned char cr = toByte(r), cg = toByte(g), cb = toByte(b);
        const int rad = std::max(1, radiusPx * ss_);
        for (int i : nodes) {
            if (i < 0 || i >= x.rows()) continue;
            Eigen::Vector3d s;
            if (!project(mvp, x.row(i).transpose(), s)) continue;
            const int cx = static_cast<int>(std::lround(s.x()));
            const int cy = static_cast<int>(std::lround(s.y()));
            for (int dy = -rad; dy <= rad; ++dy) {
                for (int dx = -rad; dx <= rad; ++dx) {
                    if (dx * dx + dy * dy > rad * rad) continue;
                    const int px = cx + dx, py = cy + dy;
                    if (px < 0 || px >= w_ || py < 0 || py >= h_) continue;
                    const size_t idx = static_cast<size_t>(py) * w_ + px;
                    if (s.z() >= depth_[idx]) continue;
                    depth_[idx] = s.z();
                    color_[3 * idx] = cr;
                    color_[3 * idx + 1] = cg;
                    color_[3 * idx + 2] = cb;
                }
            }
        }
    }

    /// 超采样缓冲 → 输出分辨率 RGB8（盒式降采样）
    void downsample(std::vector<unsigned char>& out) const {
        out.assign(static_cast<size_t>(outW_) * outH_ * 3, 0);
        const int samples = ss_ * ss_;
        for (int oy = 0; oy < outH_; ++oy) {
            for (int ox = 0; ox < outW_; ++ox) {
                unsigned sum[3] = {0, 0, 0};
                for (int sy = 0; sy < ss_; ++sy) {
                    for (int sx = 0; sx < ss_; ++sx) {
                        const size_t idx =
                            (static_cast<size_t>(oy * ss_ + sy) * w_ + (ox * ss_ + sx)) * 3;
                        sum[0] += color_[idx];
                        sum[1] += color_[idx + 1];
                        sum[2] += color_[idx + 2];
                    }
                }
                const size_t o = (static_cast<size_t>(oy) * outW_ + ox) * 3;
                out[o] = static_cast<unsigned char>(sum[0] / samples);
                out[o + 1] = static_cast<unsigned char>(sum[1] / samples);
                out[o + 2] = static_cast<unsigned char>(sum[2] / samples);
            }
        }
    }

private:
    static constexpr double kZNear = 0.01;
    static constexpr double kZFar = 1.0e4;

    static unsigned char toByte(float v) {
        const float c = std::clamp(v, 0.0f, 1.0f);
        return static_cast<unsigned char>(c * 255.0f + 0.5f);
    }

    /// 世界坐标 → 屏幕坐标（x/y 像素，z 为 NDC 深度；w<=0 视为在相机后）
    bool project(const Eigen::Matrix4d& mvp, const Vec3& p, Eigen::Vector3d& screen) const {
        const Eigen::Vector4d c = mvp * Eigen::Vector4d(p.x(), p.y(), p.z(), 1.0);
        if (c.w() <= 1e-9) return false;
        const Eigen::Vector3d ndc = c.head<3>() / c.w();
        screen.x() = (ndc.x() + 1.0) * 0.5 * (w_ - 1);
        screen.y() = (1.0 - ndc.y()) * 0.5 * (h_ - 1);
        screen.z() = ndc.z();
        return true;
    }

    /// 有向边函数（>0 表示点在 a→b 左侧）
    static double edge(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                       double px, double py) {
        return (b.x() - a.x()) * (py - a.y()) - (b.y() - a.y()) * (px - a.x());
    }

    void plot(int x, int y, double z, unsigned char r, unsigned char g, unsigned char b) {
        if (x < 0 || x >= w_ || y < 0 || y >= h_) return;
        const size_t idx = static_cast<size_t>(y) * w_ + x;
        if (z >= depth_[idx]) return;
        depth_[idx] = z;
        color_[3 * idx] = r;
        color_[3 * idx + 1] = g;
        color_[3 * idx + 2] = b;
    }

    void rasterizeTriangle(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
                           const Eigen::Vector3d& p2, float shade,
                           float r, float g, float b) {
        const double area = edge(p0, p1, p2.x(), p2.y());
        if (std::abs(area) < 1e-12) return;
        const bool ccw = area > 0.0;
        const double invArea = 1.0 / area;

        int minX = static_cast<int>(std::floor(std::min({p0.x(), p1.x(), p2.x()})));
        int maxX = static_cast<int>(std::ceil(std::max({p0.x(), p1.x(), p2.x()})));
        int minY = static_cast<int>(std::floor(std::min({p0.y(), p1.y(), p2.y()})));
        int maxY = static_cast<int>(std::ceil(std::max({p0.y(), p1.y(), p2.y()})));
        minX = std::max(minX, 0);
        maxX = std::min(maxX, w_ - 1);
        minY = std::max(minY, 0);
        maxY = std::min(maxY, h_ - 1);

        const unsigned char cr = toByte(r * shade);
        const unsigned char cg = toByte(g * shade);
        const unsigned char cb = toByte(b * shade);

        for (int y = minY; y <= maxY; ++y) {
            for (int x = minX; x <= maxX; ++x) {
                const double px = x + 0.5, py = y + 0.5;
                const double w0 = edge(p1, p2, px, py);
                const double w1 = edge(p2, p0, px, py);
                const double w2 = edge(p0, p1, px, py);
                if (ccw) {
                    if (w0 < 0.0 || w1 < 0.0 || w2 < 0.0) continue;
                } else {
                    if (w0 > 0.0 || w1 > 0.0 || w2 > 0.0) continue;
                }
                const double z = (w0 * p0.z() + w1 * p1.z() + w2 * p2.z()) * invArea;
                const size_t idx = static_cast<size_t>(y) * w_ + x;
                if (z >= depth_[idx]) continue;
                depth_[idx] = z;
                color_[3 * idx] = cr;
                color_[3 * idx + 1] = cg;
                color_[3 * idx + 2] = cb;
            }
        }
    }

    void drawLine(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1,
                  float r, float g, float b) {
        const double dx = p1.x() - p0.x();
        const double dy = p1.y() - p0.y();
        const int steps = static_cast<int>(std::ceil(std::max(std::abs(dx), std::abs(dy))));
        if (steps <= 0) return;
        const unsigned char cr = toByte(r), cg = toByte(g), cb = toByte(b);
        for (int i = 0; i <= steps; ++i) {
            const double t = static_cast<double>(i) / steps;
            const double x = p0.x() + t * dx;
            const double y = p0.y() + t * dy;
            const double z = p0.z() + t * (p1.z() - p0.z());
            // 2x2 覆盖，降采样后约 1px 实线
            const int px = static_cast<int>(std::floor(x));
            const int py = static_cast<int>(std::floor(y));
            plot(px, py, z, cr, cg, cb);
            plot(px + 1, py, z, cr, cg, cb);
            plot(px, py + 1, z, cr, cg, cb);
            plot(px + 1, py + 1, z, cr, cg, cb);
        }
    }

    int outW_, outH_, ss_;
    int w_, h_;
    std::vector<double> depth_;
    std::vector<unsigned char> color_;
    Eigen::Matrix4d view_ = Eigen::Matrix4d::Identity();
    Eigen::Matrix4d proj_ = Eigen::Matrix4d::Identity();
};

} // namespace fiction
