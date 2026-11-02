// =============================================================================
// viewer.cpp — GLFW + OpenGL 可视化模块实现（固定管线版）
// =============================================================================

#include "viewer.h"
#include "camera.h"

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <GL/gl.h>

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace fiction {

namespace {

constexpr double kMinPitch = -89.0;
constexpr double kMaxPitch = 89.0;
constexpr double kMinRadius = 0.05;
constexpr double kMaxRadius = 500.0;

/// 任意单位法向 n 的平面内两个正交基向量
void planeBasis(const Vec3& n, Vec3& t1, Vec3& t2) {
    Vec3 ref = std::abs(n(2)) < 0.9 ? Vec3(0.0, 0.0, 1.0) : Vec3(1.0, 0.0, 0.0);
    t1 = ref.cross(n).normalized();
    t2 = n.cross(t1).normalized();
}

} // anonymous namespace

// =============================================================================
// 构造 / 析构
// =============================================================================
Viewer::Viewer(const Config& cfg, const Vec3& initialTarget, double initialRadius)
    : cfg_(cfg), target_(initialTarget), radius_(initialRadius),
      initTarget_(initialTarget), initRadius_(initialRadius) {
    if (!glfwInit()) {
        std::cerr << "[Viewer] glfwInit failed\n";
        return;
    }
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_VISIBLE, GLFW_TRUE);
    window_ = glfwCreateWindow(cfg_.width, cfg_.height, cfg_.title.c_str(), nullptr, nullptr);
    if (!window_) {
        std::cerr << "[Viewer] glfwCreateWindow failed (无显示环境？)\n";
        glfwTerminate();
        return;
    }
    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1); // vsync
    glfwSetWindowUserPointer(window_, this);

    glfwSetFramebufferSizeCallback(window_, framebufferSizeCb);
    glfwSetCursorPosCallback(window_, cursorPosCb);
    glfwSetMouseButtonCallback(window_, mouseBtnCb);
    glfwSetScrollCallback(window_, scrollCb);
    glfwSetKeyCallback(window_, keyCb);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glClearColor(cfg_.bg[0], cfg_.bg[1], cfg_.bg[2], 1.0f);

    int w = 0, h = 0;
    glfwGetFramebufferSize(window_, &w, &h);
    updateMatrices(w, h);
    initOk_ = true;
}

Viewer::~Viewer() {
    if (window_) {
        glfwDestroyWindow(window_);
        glfwTerminate();
    }
}

// =============================================================================
// 生命周期
// =============================================================================
bool Viewer::ok() const { return initOk_; }

bool Viewer::shouldClose() const {
    return window_ ? glfwWindowShouldClose(window_) != 0 : true;
}

void Viewer::pollEvents() {
    if (window_) glfwPollEvents();
    const auto now = std::chrono::steady_clock::now();
    frameDt_ = std::chrono::duration<double>(now - lastFrameTime_).count();
    lastFrameTime_ = now;
}

void Viewer::beginFrame() {
    if (!window_) return;
    int w = 0, h = 0;
    glfwGetFramebufferSize(window_, &w, &h);
    if (w > 0 && h > 0) glViewport(0, 0, w, h);
    updateMatrices(w, h);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(proj_.data());
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd(view_.data());
}

void Viewer::endFrame() {
    if (window_) glfwSwapBuffers(window_);
}

void Viewer::updateMatrices(int viewportW, int viewportH) const {
    const double aspect = viewportH > 0
        ? static_cast<double>(viewportW) / static_cast<double>(viewportH) : 1.0;
    proj_ = perspective(cfg_.fovDeg, aspect, cfg_.zNear, cfg_.zFar);

    // 球坐标 → 相机位置（世界系 z 轴向上）
    const double cy = std::cos(yaw_ * kDeg2Rad), sy = std::sin(yaw_ * kDeg2Rad);
    const double cp = std::cos(pitch_ * kDeg2Rad), sp = std::sin(pitch_ * kDeg2Rad);
    const Vec3 eye = target_ + radius_ * Vec3(cp * cy, cp * sy, sp);
    view_ = lookAt(eye, target_, Vec3(0.0, 0.0, 1.0));
}

// =============================================================================
// 绘图原语
// =============================================================================
void Viewer::drawPoints(const Positions& x, float size,
                        float r, float g, float b) const {
    if (!window_ || x.rows() == 0) return;
    glPointSize(size);
    glColor3f(r, g, b);
    glBegin(GL_POINTS);
    for (int i = 0; i < x.rows(); ++i) {
        glVertex3d(x(i, 0), x(i, 1), x(i, 2));
    }
    glEnd();
}

void Viewer::drawPoints(const Positions& x, const std::vector<char>& mask,
                        float size, float r, float g, float b) const {
    if (!window_ || x.rows() == 0) return;
    glPointSize(size);
    glColor3f(r, g, b);
    glBegin(GL_POINTS);
    const int M = static_cast<int>(mask.size());
    for (int i = 0; i < x.rows(); ++i) {
        if (i < M && !mask[i]) continue;
        glVertex3d(x(i, 0), x(i, 1), x(i, 2));
    }
    glEnd();
}

void Viewer::drawEdges(const Positions& x,
                       const std::vector<std::pair<int, int>>& edges,
                       float width, float r, float g, float b) const {
    if (!window_ || x.rows() == 0) return;
    glLineWidth(width);
    glColor3f(r, g, b);
    glBegin(GL_LINES);
    for (const auto& [i0, i1] : edges) {
        if (i0 < 0 || i0 >= x.rows() || i1 < 0 || i1 >= x.rows()) continue;
        glVertex3d(x(i0, 0), x(i0, 1), x(i0, 2));
        glVertex3d(x(i1, 0), x(i1, 1), x(i1, 2));
    }
    glEnd();
}

void Viewer::drawTriangles(const Positions& x,
                           const std::vector<std::array<int, 3>>& faces,
                           float r, float g, float b) const {
    if (!window_ || x.rows() == 0) return;
    // 定向光（世界系，斜上方）；用 |dot| 实现双面可见
    const Vec3 light = Vec3(0.4, -0.6, 0.8).normalized();
    glBegin(GL_TRIANGLES);
    for (const auto& f : faces) {
        if (f[0] < 0 || f[0] >= x.rows() || f[1] < 0 || f[1] >= x.rows() ||
            f[2] < 0 || f[2] >= x.rows())
            continue;
        const Vec3 va = x.row(f[0]).transpose();
        const Vec3 vc = x.row(f[2]).transpose();
        const Vec3 vb = x.row(f[1]).transpose();
        Vec3 n = (vb - va).cross(vc - va);
        const double area2 = n.norm();
        if (area2 < 1e-24) continue; // 退化面（面积为零）跳过
        n /= area2;
        const float shade = static_cast<float>(0.35 + 0.65 * std::abs(n.dot(light)));
        glColor3f(r * shade, g * shade, b * shade);
        glVertex3d(va(0), va(1), va(2));
        glVertex3d(vb(0), vb(1), vb(2));
        glVertex3d(vc(0), vc(1), vc(2));
    }
    glEnd();
}

void Viewer::drawGround(const Vec3& origin, const Vec3& normal,
                        double extent, double spacing,
                        float r, float g, float b) const {
    if (!window_ || spacing <= 0.0) return;
    Vec3 n = normal;
    if (n.squaredNorm() < 1e-24) n = Vec3(0.0, 0.0, 1.0);
    n.normalize();
    Vec3 t1, t2;
    planeBasis(n, t1, t2);

    glLineWidth(1.0f);
    glColor3f(r, g, b);
    glBegin(GL_LINES);
    for (double s = -extent; s <= extent + 1e-12; s += spacing) {
        const Vec3 a = origin + t1 * s;
        const Vec3 b1 = origin - t2 * extent, b2 = origin + t2 * extent;
        glVertex3d(a(0), a(1), a(2));
        glVertex3d(b1(0), b1(1), b1(2));
        glVertex3d(a(0), a(1), a(2));
        glVertex3d(b2(0), b2(1), b2(2));
        const Vec3 c = origin + t2 * s;
        const Vec3 d1 = origin - t1 * extent, d2 = origin + t1 * extent;
        glVertex3d(c(0), c(1), c(2));
        glVertex3d(d1(0), d1(1), d1(2));
        glVertex3d(c(0), c(1), c(2));
        glVertex3d(d2(0), d2(1), d2(2));
    }
    glEnd();
}

// =============================================================================
// 交互状态
// =============================================================================
bool Viewer::resetRequested() {
    const bool req = resetReq_;
    resetReq_ = false;
    return req;
}

// =============================================================================
// GLFW 回调
// =============================================================================
void Viewer::framebufferSizeCb(GLFWwindow* w, int width, int height) {
    (void)w;
    if (width > 0 && height > 0) glViewport(0, 0, width, height);
}

void Viewer::cursorPosCb(GLFWwindow* w, double xpos, double ypos) {
    Viewer* self = static_cast<Viewer*>(glfwGetWindowUserPointer(w));
    if (!self) return;
    const double dx = xpos - self->lastX_;
    const double dy = ypos - self->lastY_;
    self->lastX_ = xpos;
    self->lastY_ = ypos;
    if (!self->leftDrag_ && !self->rightDrag_) return;

    if (self->leftDrag_) {
        self->yaw_ += dx * 0.3;
        self->pitch_ += dy * 0.3;
        self->pitch_ = std::clamp(self->pitch_, kMinPitch, kMaxPitch);
    } else { // 右键平移：沿相机右/上方向移动注视点
        const double scale = self->radius_ * 0.0015;
        const double cy = std::cos(self->yaw_ * kDeg2Rad);
        const double sy = std::sin(self->yaw_ * kDeg2Rad);
        const double cp = std::cos(self->pitch_ * kDeg2Rad);
        const double sp = std::sin(self->pitch_ * kDeg2Rad);
        const Vec3 right(-sy, cy, 0.0);
        const Vec3 forward = -Vec3(cp * cy, cp * sy, sp);
        const Vec3 up = right.cross(forward);
        self->target_ += right * (-dx * scale) + up * (dy * scale);
    }
}

void Viewer::mouseBtnCb(GLFWwindow* w, int button, int action, int mods) {
    Viewer* self = static_cast<Viewer*>(glfwGetWindowUserPointer(w));
    if (!self) return;
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        self->leftDrag_ = (action == GLFW_PRESS) && !(mods & GLFW_MOD_SHIFT);
        self->rightDrag_ = (action == GLFW_PRESS) && (mods & GLFW_MOD_SHIFT);
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        self->rightDrag_ = (action == GLFW_PRESS);
    }
    glfwGetCursorPos(w, &self->lastX_, &self->lastY_);
}

void Viewer::scrollCb(GLFWwindow* w, double, double yoff) {
    Viewer* self = static_cast<Viewer*>(glfwGetWindowUserPointer(w));
    if (!self) return;
    self->radius_ *= std::pow(0.9, yoff);
    self->radius_ = std::clamp(self->radius_, kMinRadius, kMaxRadius);
}

void Viewer::keyCb(GLFWwindow* w, int key, int, int action, int) {
    Viewer* self = static_cast<Viewer*>(glfwGetWindowUserPointer(w));
    if (!self || action != GLFW_PRESS) return;
    switch (key) {
        case GLFW_KEY_ESCAPE:
            glfwSetWindowShouldClose(w, GLFW_TRUE);
            break;
        case GLFW_KEY_SPACE:
            self->paused_ = !self->paused_;
            break;
        case GLFW_KEY_R:
            self->target_ = self->initTarget_;
            self->radius_ = self->initRadius_;
            self->yaw_ = self->initYaw_;
            self->pitch_ = self->initPitch_;
            self->resetReq_ = true;
            break;
        case GLFW_KEY_UP:
            self->stepsPerFrame_ = std::min(self->stepsPerFrame_ + 1, 64);
            break;
        case GLFW_KEY_DOWN:
            self->stepsPerFrame_ = std::max(self->stepsPerFrame_ - 1, 0);
            break;
        default:
            break;
    }
}

} // namespace fiction
