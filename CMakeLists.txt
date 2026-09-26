cmake_minimum_version_required(4.0)

project(AOCS VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

target_include_directories(AOCS PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(test)