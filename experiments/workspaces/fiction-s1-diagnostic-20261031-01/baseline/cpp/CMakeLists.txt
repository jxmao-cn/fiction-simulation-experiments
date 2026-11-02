cmake_minimum_required(VERSION 3.20)
project(fiction_sim CXX)

# -----------------------------------------------------------------------------
#  C++ 标准与构建类型（用户要求可升级至 C++20，便于 Concept）
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "[fiction] C++ standard = ${CMAKE_CXX_STANDARD}, buildType=${CMAKE_BUILD_TYPE}")

# -----------------------------------------------------------------------------
#  选项：是否启用 CHOLMOD（Eigen 封装 + CHOLMOD 计算主线，可选回退）
#  - macOS: 通过 brew --prefix suite-sparse 自动探测
#  - Ubuntu: 通过 find_library / find_path 探测 /usr/include/suitesparse
#  - 保持 Eigen 封装主线：无论是否启用 CHOLMOD，代码均走 Eigen::SparseMatrix 组装
# -----------------------------------------------------------------------------
option(FICTION_USE_CHOLMOD "Enable Eigen::CholmodSupernodalLLT (CHOLMOD) for Newton solve" ON)

# =============================================================================
#  Eigen3 探测
#  - 优先 find_package(CONFIG)，失败则回退到 brew --prefix eigen / 系统路径
#  - Eigen 为 header-only，无需链接，仅需 include
# =============================================================================
set(EIGEN_FOUND FALSE)
find_package(Eigen3 3.4 QUIET CONFIG)
if(Eigen3_FOUND)
    message(STATUS "[fiction] Eigen3 found via find_package CONFIG: ${Eigen3_VERSION} @ ${Eigen3_DIR}")
    set(EIGEN_FOUND TRUE)
else()
    find_package(Eigen3 3.4 QUIET)
    if(Eigen3_FOUND)
        message(STATUS "[fiction] Eigen3 found via find_package: ${Eigen3_VERSION}")
        set(EIGEN_FOUND TRUE)
    endif()
endif()

if(NOT EIGEN_FOUND)
    # 回退：brew / 系统路径探测
    execute_process(
        COMMAND brew --prefix eigen
        OUTPUT_VARIABLE EIGEN_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(EXISTS "${EIGEN_PREFIX}/include/eigen3")
        set(EIGEN3_INCLUDE_DIR "${EIGEN_PREFIX}/include/eigen3")
        message(STATUS "[fiction] Eigen3 found via brew: ${EIGEN3_INCLUDE_DIR}")
        set(EIGEN_FOUND TRUE)
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIR})
    elseif(EXISTS "/usr/include/eigen3")
        set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3")
        message(STATUS "[fiction] Eigen3 found at /usr/include/eigen3 (Ubuntu apt)")
        set(EIGEN_FOUND TRUE)
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIR})
    elseif(EXISTS "/opt/homebrew/include/eigen3")
        set(EIGEN3_INCLUDE_DIR "/opt/homebrew/include/eigen3")
        message(STATUS "[fiction] Eigen3 found at /opt/homebrew/include/eigen3")
        set(EIGEN_FOUND TRUE)
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        target_include_directories(Eigen3::Eigen INTERFACE ${EIGEN3_INCLUDE_DIR})
    else()
        message(FATAL_ERROR "[fiction] Eigen3 not found. macOS: brew install eigen | Ubuntu: sudo apt install libeigen3-dev")
    endif()
endif()

# =============================================================================
#  zlib 探测（可选）：视频帧 PNG 导出用；缺失时 image_writer 回退 PPM
# =============================================================================
find_package(ZLIB QUIET)
if(ZLIB_FOUND)
    message(STATUS "[fiction] zlib found, PNG frame export enabled")
else()
    message(STATUS "[fiction] zlib not found, PNG frame export falls back to PPM")
endif()

# =============================================================================
#  SuiteSparse / CHOLMOD 探测（可选，跨平台 macOS + Ubuntu）
# =============================================================================
set(FICTION_HAS_CHOLMOD FALSE)
set(SUITESPARSE_PREFIX "")
set(SUITESPARSE_INCLUDE_DIR "")
set(SUITESPARSE_LIB_DIR "")

if(FICTION_USE_CHOLMOD)
    # 尝试 brew 路径（macOS）
    execute_process(
        COMMAND brew --prefix suite-sparse
        OUTPUT_VARIABLE SUITESPARSE_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(SUITESPARSE_PREFIX AND EXISTS "${SUITESPARSE_PREFIX}/include/suitesparse/cholmod.h")
        set(SUITESPARSE_INCLUDE_DIR "${SUITESPARSE_PREFIX}/include")
        set(SUITESPARSE_LIB_DIR "${SUITESPARSE_PREFIX}/lib")
        message(STATUS "[fiction] SuiteSparse via brew: ${SUITESPARSE_PREFIX}")
        set(FICTION_HAS_CHOLMOD TRUE)
    elseif(EXISTS "/opt/homebrew/include/suitesparse/cholmod.h")
        set(SUITESPARSE_INCLUDE_DIR "/opt/homebrew/include")
        set(SUITESPARSE_LIB_DIR "/opt/homebrew/lib")
        message(STATUS "[fiction] SuiteSparse at /opt/homebrew (macOS)")
        set(FICTION_HAS_CHOLMOD TRUE)
    elseif(EXISTS "/usr/include/suitesparse/cholmod.h")
        # Ubuntu apt: libsuitesparse-dev
        set(SUITESPARSE_INCLUDE_DIR "/usr/include")
        set(SUITESPARSE_LIB_DIR "/usr/lib/x86_64-linux-gnu")
        if(NOT EXISTS "${SUITESPARSE_LIB_DIR}/libcholmod.so")
            set(SUITESPARSE_LIB_DIR "/usr/lib")
        endif()
        message(STATUS "[fiction] SuiteSparse at /usr (Ubuntu apt)")
        set(FICTION_HAS_CHOLMOD TRUE)
    else()
        # 最后尝试 find_library
        find_path(CHOLMOD_INCLUDE_DIR cholmod.h PATH_SUFFIXES suitesparse)
        find_library(CHOLMOD_LIB cholmod)
        if(CHOLMOD_INCLUDE_DIR AND CHOLMOD_LIB)
            get_filename_component(SUITESPARSE_INCLUDE_DIR "${CHOLMOD_INCLUDE_DIR}" DIRECTORY)
            get_filename_component(SUITESPARSE_LIB_DIR "${CHOLMOD_LIB}" DIRECTORY)
            message(STATUS "[fiction] SuiteSparse via find_library: inc=${CHOLMOD_INCLUDE_DIR} lib=${CHOLMOD_LIB}")
            set(FICTION_HAS_CHOLMOD TRUE)
        else()
            message(WARNING "[fiction] FICTION_USE_CHOLMOD=ON but SuiteSparse not found. Fallback to Eigen SimplicialLDLT (pure Eigen). "
                            "Install: macOS `brew install suite-sparse` | Ubuntu `sudo apt install libsuitesparse-dev`")
            set(FICTION_HAS_CHOLMOD FALSE)
        endif()
    endif()
endif()

if(FICTION_HAS_CHOLMOD)
    add_library(suitesparse INTERFACE)
    # 同时暴露 .../include 与 .../include/suitesparse，以便 <cholmod.h> 与 <suitesparse/cholmod.h> 均可
    target_include_directories(suitesparse INTERFACE ${SUITESPARSE_INCLUDE_DIR} ${SUITESPARSE_INCLUDE_DIR}/suitesparse)
    if(SUITESPARSE_LIB_DIR)
        target_link_directories(suitesparse INTERFACE ${SUITESPARSE_LIB_DIR})
    endif()
    # 常见 CHOLMOD 依赖链：cholmod + AMD/CAMD/CCOLAMD/COLAMD/SuiteSparseConfig
    # 不同系统库名略有差异，缺失的链接项会由 find 阶段报错，这里保持完整
    target_link_libraries(suitesparse INTERFACE cholmod amd camd ccolamd colamd suitesparseconfig)
    # 额外常见依赖（Ubuntu 下 cholmod 依赖 blas/lapack，若链接失败可取消注释）
    # find_package(BLAS QUIET)
    # if(BLAS_FOUND)
    #     target_link_libraries(suitesparse INTERFACE ${BLAS_LIBRARIES})
    # endif()
    message(STATUS "[fiction] CHOLMOD enabled: inc=${SUITESPARSE_INCLUDE_DIR} lib=${SUITESPARSE_LIB_DIR}")
else()
    message(STATUS "[fiction] CHOLMOD disabled -> using Eigen::SimplicialLDLT (header-only fallback)")
endif()

# =============================================================================
#  核心库 fiction_core：time_integrator + utils
#  - 对外暴露 include/ 目录
#  - 依赖 Eigen（必选） + SuiteSparse（可选）
#  - 编译定义 FICTION_HAS_CHOLMOD / FICTION_USE_CHOLMOD 供 src/time_integrator.cpp 切换求解器
# =============================================================================
add_library(fiction_core STATIC
    src/utils.cpp
    src/time_integrator.cpp
)
target_include_directories(fiction_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(fiction_core PUBLIC Eigen3::Eigen)
if(FICTION_HAS_CHOLMOD)
    target_link_libraries(fiction_core PUBLIC suitesparse)
    target_compile_definitions(fiction_core PUBLIC FICTION_HAS_CHOLMOD=1)
    if(FICTION_USE_CHOLMOD)
        target_compile_definitions(fiction_core PUBLIC FICTION_USE_CHOLMOD=1)
    endif()
endif()
# 警告与优化
target_compile_options(fiction_core PRIVATE -Wall -Wextra -Wpedantic)

# =============================================================================
#  可执行：cholmod_test（保留原始裸 CHOLMOD 验证，仅当 CHOLMOD 可用时编译）
# =============================================================================
if(FICTION_HAS_CHOLMOD AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cholmod_test.cpp")
    add_executable(cholmod_test cholmod_test.cpp)
    target_link_libraries(cholmod_test PRIVATE Eigen3::Eigen suitesparse)
endif()

# =============================================================================
#  可执行：simulator（控制台验证，对标 Python_slovement/simulator.py，无可视化）
#  - 演示 TimeIntegrator + Inertia 内置 + IEnergy 扩展
#  - 仅打印 residual / step size / 位置变化
# =============================================================================
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/simulator.cpp")
    add_executable(simulator examples/simulator.cpp)
    target_link_libraries(simulator PRIVATE fiction_core)
endif()

# =============================================================================
#  可执行：barrier_check（BarrierEnergy 数值验证：FD 梯度/海森 + CCD + 落地场景）
# =============================================================================
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/barrier_check.cpp")
    add_executable(barrier_check examples/barrier_check.cpp)
    target_link_libraries(barrier_check PRIVATE fiction_core)
endif()

# =============================================================================
#  可执行：ipc_mesh_demo（IPC .msh 模型加载 + mass-spring 仿真 + OBJ/视频导出）
# =============================================================================
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/ipc_mesh_demo.cpp")
    add_executable(ipc_mesh_demo examples/ipc_mesh_demo.cpp)
    target_link_libraries(ipc_mesh_demo PRIVATE fiction_core)
    if(ZLIB_FOUND)
        target_link_libraries(ipc_mesh_demo PRIVATE ZLIB::ZLIB)
        target_compile_definitions(ipc_mesh_demo PRIVATE FICTION_HAS_ZLIB=1)
    endif()
endif()

# =============================================================================
#  可执行：mesh_check（msh 解析 + mass-spring/外力数值验证，供 CI）
# =============================================================================
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/mesh_check.cpp")
    add_executable(mesh_check examples/mesh_check.cpp)
    target_link_libraries(mesh_check PRIVATE fiction_core)
    target_compile_definitions(mesh_check PRIVATE
        FICTION_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
endif()

# =============================================================================
#  可视化模块（GLFW + OpenGL 固定管线，可选）
#  - FICTION_USE_VIEWER=ON 时探测 glfw3 + OpenGL，找到则构建 fiction_viewer 库
#    与 simulator_visual 可执行；找不到则警告并跳过，不影响主线编译
# =============================================================================
option(FICTION_USE_VIEWER "Build GLFW+OpenGL real-time viewer module" ON)
set(FICTION_HAS_VIEWER FALSE)
set(FICTION_GLFW_TARGET "")
if(FICTION_USE_VIEWER)
    find_package(OpenGL QUIET)
    find_package(glfw3 3.3 QUIET CONFIG)
    if(TARGET glfw)
        set(FICTION_GLFW_TARGET glfw)
    else()
        find_library(GLFW3_LIB glfw3)
        find_path(GLFW3_INCLUDE_DIR GLFW/glfw3.h)
        if(GLFW3_LIB AND GLFW3_INCLUDE_DIR)
            add_library(glfw3::glfw UNKNOWN IMPORTED)
            set_target_properties(glfw3::glfw PROPERTIES
                IMPORTED_LOCATION "${GLFW3_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES "${GLFW3_INCLUDE_DIR}")
            set(FICTION_GLFW_TARGET glfw3::glfw)
        endif()
    endif()
    if(OPENGL_FOUND AND FICTION_GLFW_TARGET)
        set(FICTION_HAS_VIEWER TRUE)
        message(STATUS "[fiction] Viewer enabled: OpenGL + GLFW3 found")
    else()
        message(WARNING "[fiction] FICTION_USE_VIEWER=ON but glfw3/OpenGL not found, viewer skipped. "
                        "Install: Ubuntu `sudo apt install libglfw3-dev libgl1-mesa-dev`")
    endif()
endif()

if(FICTION_HAS_VIEWER)
    add_library(fiction_viewer STATIC
        src/viewer.cpp
    )
    target_include_directories(fiction_viewer PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(fiction_viewer PUBLIC fiction_core OpenGL::GL ${FICTION_GLFW_TARGET})
    target_compile_options(fiction_viewer PRIVATE -Wall -Wextra -Wpedantic)

    # 可执行：simulator_visual（交互式实时可视化 demo）
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/simulator_visual.cpp")
        add_executable(simulator_visual examples/simulator_visual.cpp)
        target_link_libraries(simulator_visual PRIVATE fiction_viewer)
    endif()

    # 可执行：ipc_mesh_visual（IPC .msh 模型交互式可视化）
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/ipc_mesh_visual.cpp")
        add_executable(ipc_mesh_visual examples/ipc_mesh_visual.cpp)
        target_link_libraries(ipc_mesh_visual PRIVATE fiction_viewer)
    endif()
endif()

# =============================================================================
#  可执行：time_integrator_demo（最小 Newton 验证，用于 CI）
# =============================================================================
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/time_integrator_demo.cpp")
    add_executable(time_integrator_demo examples/time_integrator_demo.cpp)
    target_link_libraries(time_integrator_demo PRIVATE fiction_core)
endif()

# =============================================================================
#  安装与测试（预留）
# =============================================================================
enable_testing()
# 阶段零：材料实现与动态离散参考的独立验收入口。
add_executable(material_validation examples/material_validation.cpp)
target_link_libraries(material_validation PRIVATE fiction_core)
add_test(NAME stage0_material_validation COMMAND material_validation)
add_executable(dynamic_reference_check examples/dynamic_reference_check.cpp)
target_link_libraries(dynamic_reference_check PRIVATE fiction_core)
add_test(NAME stage0_dynamic_reference COMMAND dynamic_reference_check)
add_executable(integrator_report_check examples/integrator_report_check.cpp)
target_link_libraries(integrator_report_check PRIVATE fiction_core)
add_test(NAME integrator_report_check COMMAND integrator_report_check)
add_executable(tension_check examples/tension_check.cpp)
target_link_libraries(tension_check PRIVATE fiction_core)
if(TARGET simulator)
    add_test(NAME simulator_smoke COMMAND simulator)
    set_tests_properties(simulator_smoke PROPERTIES FAIL_REGULAR_EXPRESSION "\\[FAIL\\]")
endif()
if(TARGET barrier_check)
    add_test(NAME barrier_check COMMAND barrier_check)
endif()
if(TARGET mesh_check)
    add_test(NAME mesh_check COMMAND mesh_check)
endif()
if(TARGET ipc_mesh_demo)
    add_test(NAME ipc_mesh_demo_smoke
             COMMAND ipc_mesh_demo --mesh ${CMAKE_CURRENT_SOURCE_DIR}/tests/data/cube.msh
                     --no-gravity --steps 2 --tol 1e-6)
endif()
if(TARGET simulator_visual)
    add_test(NAME simulator_visual_build_only
             COMMAND ${CMAKE_COMMAND} -E echo "GUI demo, no automated test")
endif()
if(TARGET ipc_mesh_visual)
    add_test(NAME ipc_mesh_visual_build_only
             COMMAND ${CMAKE_COMMAND} -E echo "GUI demo, no automated test")
endif()
