#pragma once
// =============================================================================
// viewer.h — GLFW + OpenGL 可视化模块（固定管线版）
// =============================================================================
// 职责：
//   - 独立于 TimeIntegrator / IEnergy，只消费 Positions [N,3] 等只读数据
//   - 提供窗口生命周期：beginFrame / endFrame / pollEvents / shouldClose
//   - 提供绘图原语：drawPoints / drawEdges / drawGround
//   - 轨道相机：左键拖拽旋转，滚轮缩放，右键拖拽平移
//   - 交互状态供外层仿真循环读取：paused() / stepsPerFrame() / resetRequested()
//     键位：Space 暂停，R 复位相机，Up/Down 每帧步数 ±1，ESC 关闭窗口
//
// 渲染方式：OpenGL 兼容管线（固定管线），矩阵用 Eigen 自算后 glLoadMatrixd，
// 无 GLM / glad 依赖。教学型项目，网格规模小，立即模式足够。
//
// 使用约定：单窗口单实例（构造时 glfwInit，析构时 glfwTerminate）。
// =============================================================================

#include "types.h"

#include <array>
#include <chrono>
#include <string>
#include <utility>
#include <vector>

struct GLFWwindow;

namespace fiction {

class Viewer {
public:
    struct Config {
        int width = 1280;                ///< 初始窗口宽
        int height = 720;                ///< 初始窗口高
        std::string title = "fiction viewer";
        float bg[3] = {0.08f, 0.09f, 0.11f};  ///< 背景色 RGB
        double fovDeg = 45.0;            ///< 视场角（度）
        double zNear = 0.01;
        double zFar = 1000.0;
    };

    /// 默认配置构造
    Viewer()
        : Viewer(Config(), Vec3(0.0, 0.0, 0.5), 6.0) {}
    /// @param initialTarget 相机初始注视点
    /// @param initialRadius 初始轨道半径（相机到注视点距离）
    explicit Viewer(const Config& cfg,
                    const Vec3& initialTarget,
                    double initialRadius);
    ~Viewer();

    Viewer(const Viewer&) = delete;
    Viewer& operator=(const Viewer&) = delete;

    // --- 生命周期 ---
    bool ok() const;             ///< 窗口/OpenGL 初始化是否成功
    bool shouldClose() const;    ///< 用户请求关闭（叉掉窗口或 ESC）
    void pollEvents();           ///< 处理输入事件（每帧一次）
    void beginFrame();           ///< 清屏 + 应用相机矩阵（每帧绘制前调用）
    void endFrame();             ///< 交换缓冲（每帧绘制后调用）

    // --- 绘图原语（单位与 Positions 一致，颜色为 0..1 RGB） ---
    void drawPoints(const Positions& x, float size,
                    float r, float g, float b) const;
    /// 仅绘制 mask 中非零的节点（如 DBC 节点高亮）
    void drawPoints(const Positions& x, const std::vector<char>& mask, float size,
                    float r, float g, float b) const;
    void drawEdges(const Positions& x,
                   const std::vector<std::pair<int, int>>& edges,
                   float width, float r, float g, float b) const;
    /// 三角面绘制（逐面法线做定向光 lambert 明暗，双面可见；用于 IPC 模型表面）
    void drawTriangles(const Positions& x,
                       const std::vector<std::array<int, 3>>& faces,
                       float r, float g, float b) const;
    /// 以 origin 为中心、normal 为法向的地面网格（如 z=0 平面网格）
    void drawGround(const Vec3& origin, const Vec3& normal,
                    double extent, double spacing,
                    float r, float g, float b) const;

    // --- 外层仿真循环读取的交互状态 ---
    bool paused() const { return paused_; }
    int stepsPerFrame() const { return stepsPerFrame_; }
    bool resetRequested();       ///< 读取并清除"复位相机"请求
    double frameDt() const { return frameDt_; }  ///< 上一帧耗时（秒），可控制步进节奏

private:
    // --- GLFW 回调（静态转发到成员） ---
    static void framebufferSizeCb(GLFWwindow* w, int width, int height);
    static void cursorPosCb(GLFWwindow* w, double xpos, double ypos);
    static void mouseBtnCb(GLFWwindow* w, int button, int action, int mods);
    static void scrollCb(GLFWwindow* w, double xoff, double yoff);
    static void keyCb(GLFWwindow* w, int key, int scancode, int action, int mods);

    void updateMatrices(int viewportW, int viewportH) const;  ///< 重算投影/视图矩阵

    Config cfg_;
    GLFWwindow* window_ = nullptr;
    bool initOk_ = false;

    // 轨道相机状态
    Vec3 target_ = Vec3(0.0, 0.0, 0.5);  ///< 注视点
    double radius_ = 6.0;                ///< 轨道半径
    double yaw_ = 45.0;                  ///< 方位角（度）
    double pitch_ = 25.0;                ///< 俯仰角（度）
    Vec3 initTarget_;                    ///< R 键复位用初始值
    double initRadius_ = 6.0;
    double initYaw_ = 45.0;
    double initPitch_ = 25.0;
    mutable Eigen::Matrix4d proj_ = Eigen::Matrix4d::Identity();
    mutable Eigen::Matrix4d view_ = Eigen::Matrix4d::Identity();

    // 鼠标/键盘交互状态
    bool leftDrag_ = false;
    bool rightDrag_ = false;
    double lastX_ = 0.0;
    double lastY_ = 0.0;
    bool paused_ = false;
    int stepsPerFrame_ = 1;
    bool resetReq_ = false;
    double frameDt_ = 0.0;
    std::chrono::steady_clock::time_point lastFrameTime_ =
        std::chrono::steady_clock::now();
};

} // namespace fiction
