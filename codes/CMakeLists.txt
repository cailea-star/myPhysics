# ==================== 项目配置 ====================
# 仅生成 compile_commands.json 供 clangd 使用
cmake_minimum_required(VERSION 3.15)


# ==================== vcpkg 路径配置 ====================
# 根据平台自动设置 vcpkg 路径和 triplet (寻找vcpkg)
if(WIN32)
    set(CMAKE_TOOLCHAIN_FILE "$ENV{USERPROFILE}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "vcpkg toolchain file")
    set(VCPKG_TARGET_TRIPLET "x64-mingw-static" CACHE STRING "vcpkg target triplet")
else()
    set(CMAKE_TOOLCHAIN_FILE "$ENV{HOME}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "vcpkg toolchain file")
    set(VCPKG_TARGET_TRIPLET "x64-linux" CACHE STRING "vcpkg target triplet")
endif()


# ==================== 编译设置 ====================
project(UTILS LANGUAGES CXX)                            # 定义项目名称和使用的语言
set(CMAKE_CXX_STANDARD 17)                              # 设置 C++ 标准为 C++17
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)                   # 生成 compile_commands.json 供 clangd 使用
enable_testing()                                        # 启用 CTest 测试功能


# ==================== 搜索依赖库 ====================
find_package(Eigen3 REQUIRED)
find_package(GSL REQUIRED)
if(MINGW AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")                                   # Clang 沿用 MinGW 的 GNU OpenMP 运行库
    set(OpenMP_CXX_FLAGS "-fopenmp=libgomp" CACHE STRING "OpenMP flags" FORCE)         # 避免 Clang 默认链接缺失的 libomp
endif()
find_package(OpenMP REQUIRED)


# ================ 源代码编译 + PUBLIC 传递依赖 ====================
# 收集 cc 库的源文件, 配置库目标
file(GLOB UTILS_CPPS CONFIGURE_DEPENDS "utils/*.cpp")
add_library(mylib INTERFACE)                                        # 添加库目标, INTERFACE or SHARED 模式
target_include_directories(mylib INTERFACE utils)                   # 传递头文件搜索路径, INTERFACE or PUBLIC 模式
target_link_libraries(mylib INTERFACE Eigen3::Eigen)                # 传递依赖库, INTERFACE or PUBLIC 模式
target_link_libraries(mylib INTERFACE GSL::gsl GSL::gslcblas)       # 传递依赖库, INTERFACE or PUBLIC 模式
target_link_libraries(mylib INTERFACE OpenMP::OpenMP_CXX)           # 传递依赖库, INTERFACE or PUBLIC 模式


# ==================== 复制 json ====================
add_custom_target(copy_compile_commands ALL                                             # 默认构建时同步编译数据库
    COMMAND ${CMAKE_COMMAND} -E copy_if_different                                       # 仅在内容变化时复制
        ${CMAKE_BINARY_DIR}/compile_commands.json                                       # 源文件
        ${CMAKE_SOURCE_DIR}/compile_commands.json                                       # 目标文件
    COMMENT "Copied compile_commands.json to source directory"                          # 打印说明
)


# ==================== 编译可执行文件 ====================
file(GLOB TEST_CPPS CONFIGURE_DEPENDS "test/*_test.cpp")                                # 收集测试源文件
foreach(TEST_FILE ${TEST_CPPS})                                                         # 遍历每个测试源文件
    get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)                              # 提取测试文件名(NAME_WE模式)->TEST_NAME
    add_executable(${TEST_NAME} ${TEST_FILE})                                           # 配置可执行文件目标
    target_link_libraries(${TEST_NAME} PRIVATE mylib)                                   # 链接库(依赖自动传递)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})                                    # 将测试可执行文件注册到 CTest
endforeach()
