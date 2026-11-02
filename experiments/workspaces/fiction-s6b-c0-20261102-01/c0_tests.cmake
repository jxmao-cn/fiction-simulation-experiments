add_executable(fiction_s6b_c0_unit "${CMAKE_CURRENT_SOURCE_DIR}/../c0_unit.cpp")
target_link_libraries(fiction_s6b_c0_unit PRIVATE ipc::toolkit)
add_test(NAME fiction_s6b_c0_unit COMMAND fiction_s6b_c0_unit)

# C0只编译上游barrier/distance/基础CCD测试源，不包含物理或基准大测试。
get_target_property(C0_ALL_TEST_SOURCES ipc_toolkit_tests SOURCES)
set(C0_SELECTED_TEST_SOURCES)
foreach(SOURCE IN LISTS C0_ALL_TEST_SOURCES)
  if(SOURCE MATCHES "/(main|utils)\\.(cpp|hpp)$" OR
     SOURCE MATCHES "/config\\.hpp$" OR
     SOURCE MATCHES "/barrier/test_barrier\\.cpp$" OR
     SOURCE MATCHES "/distance/test_[^/]*\\.cpp$" OR
     SOURCE MATCHES "/ccd/(test_point_point_ccd|test_point_edge_ccd|test_point_triangle_ccd|test_edge_edge_ccd|collision_generator)\\.(cpp|hpp)$")
    list(APPEND C0_SELECTED_TEST_SOURCES "${SOURCE}")
  endif()
endforeach()
set_property(TARGET ipc_toolkit_tests PROPERTY SOURCES "${C0_SELECTED_TEST_SOURCES}")
file(WRITE "${CMAKE_BINARY_DIR}/c0-selected-upstream-sources.txt" "${C0_SELECTED_TEST_SOURCES}\n")

set(C0_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")
while(NOT EXISTS "${C0_ROOT}/fiction.project.json")
  get_filename_component(C0_ROOT "${C0_ROOT}" DIRECTORY)
endwhile()
file(READ "${C0_ROOT}/fiction.project.json" C0_LAYOUT)
string(JSON C0_SIMULATOR GET "${C0_LAYOUT}" simulator)
set(FICTION_USE_VIEWER OFF CACHE BOOL "C0 has no viewer" FORCE)
add_subdirectory("${C0_ROOT}/${C0_SIMULATOR}/cpp" "${CMAKE_BINARY_DIR}/production")
