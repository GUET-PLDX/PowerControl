cmake_minimum_required(VERSION 3.22)
project(power_control_host_tests LANGUAGES C CXX)

get_filename_component(MODULES_OWNER_ROOT "${CMAKE_CURRENT_LIST_DIR}/../.." ABSOLUTE)
get_filename_component(WORKSPACE_ROOT "${MODULES_OWNER_ROOT}/../.." ABSOLUTE)
if(NOT LIBXR_ROOT)
  set(LIBXR_ROOT "${WORKSPACE_ROOT}/Middlewares/Third_Party/LibXR")
endif()

set(LIBXR_SYSTEM linux)
set(LIBXR_DEFAULT_SCALAR float)
add_subdirectory("${LIBXR_ROOT}" "${CMAKE_CURRENT_BINARY_DIR}/libxr")
target_include_directories(xr PUBLIC "${LIBXR_ROOT}/driver/linux")
target_sources(xr PRIVATE "${LIBXR_ROOT}/driver/linux/linux_timebase.cpp")

add_executable(power_control_integration_test power_control_integration_test.cpp)
target_compile_features(power_control_integration_test PRIVATE cxx_std_20)
target_compile_options(power_control_integration_test PRIVATE -UNDEBUG)
target_include_directories(power_control_integration_test PRIVATE
  "${MODULES_OWNER_ROOT}/PowerControl"
  "${MODULES_OWNER_ROOT}/SuperPower"
  "${MODULES_OWNER_ROOT}/Referee"
  "${MODULES_OWNER_ROOT}/CMD"
)
target_link_libraries(power_control_integration_test PRIVATE xr pthread)

enable_testing()
add_test(NAME power_control_integration COMMAND power_control_integration_test)
