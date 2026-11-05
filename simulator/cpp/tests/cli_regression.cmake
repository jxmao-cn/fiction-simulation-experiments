# 使用真实 CLI 检查返回码与诊断，防止错误输入或未收敛被误报为成功。
if(NOT DEFINED TOOL OR NOT DEFINED MESH)
    message(FATAL_ERROR "TOOL and MESH are required")
endif()

function(expect_failure diagnostic)
    execute_process(
        COMMAND "${TOOL}" --mesh "${MESH}" ${ARGN}
        RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error
        TIMEOUT 10)
    if("${result}" STREQUAL "0" OR NOT "${output}${error}" MATCHES "${diagnostic}")
        message(FATAL_ERROR "Expected failure matching '${diagnostic}' for ${ARGN}; result=${result}\n${output}${error}")
    endif()
endfunction()

# 用同名目录阻止文件写入；不依赖文件系统权限（root 下 chmod 负例无效）。
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef suffix)
set(scratch "${CMAKE_CURRENT_BINARY_DIR}/cli-output-${suffix}")
file(MAKE_DIRECTORY "${scratch}/obj-initial/0.obj" "${scratch}/obj-later/1.obj")
expect_failure("OBJ 写出失败" --steps 0 --export "${scratch}/obj-initial")
expect_failure("OBJ 写出失败" --steps 1 --no-gravity --export-every 1 --export "${scratch}/obj-later")

# 同时覆盖 PNG/zlib 和无 zlib 的 PPM 构建，失败前不应进入 ffmpeg 编码。
file(MAKE_DIRECTORY "${scratch}/video.mp4.frames/frame_000000.png"
                    "${scratch}/video.mp4.frames/frame_000000.ppm")
expect_failure("帧写出失败" --steps 0 --video-size 16 16 --video "${scratch}/video.mp4")
file(REMOVE_RECURSE "${scratch}")

# 超出 int 范围的步数过去会窄化为 0；NaN/Inf 过去会绕过大小比较。
expect_failure("参数错误" --steps 4294967296)
expect_failure("参数错误" --steps 2147483648)
expect_failure("参数错误" --steps -1)
expect_failure("参数错误" --dt nan)
expect_failure("参数错误" --poisson nan)
expect_failure("参数错误" --video-fps inf)
expect_failure("参数错误" --init-vel nan 0 0)
expect_failure("参数错误" --gravity 0 0 -inf)
expect_failure("参数错误" --pin-box nan 0 0 1 1 1)
expect_failure("参数错误" --tol 0)
expect_failure("参数错误" --youngs -1)

# 给真实的非线性场景仅一次 Newton 迭代；应在第一步失败并立即停止。
expect_failure("step=1.*stop_reason=maxIter.*completed_steps=0"
    --steps 2 --max-iter 1 --tol 1e-12
    --pin-box -0.1 -0.1 -0.1 0.1 1.1 1.1
    --force-box 0.9 -0.1 -0.1 1.1 1.1 1.1 100 0 0)

# 零步是合法的场景检查；NH 材料标题不能继续标成 mass-spring。
execute_process(
    COMMAND "${TOOL}" --mesh "${MESH}" --material nh --steps 0 --no-gravity
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 10)
if(NOT "${result}" STREQUAL "0" OR NOT "${output}" MATCHES "Neo-Hookean FEM")
    message(FATAL_ERROR "Valid NH scene failed or wrong material label: ${result}\n${output}${error}")
endif()
