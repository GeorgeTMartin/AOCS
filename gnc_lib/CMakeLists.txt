add_library(gnc_lib STATIC
    main.cpp
    src/conversions/frames/FK5.cpp
    src/navigation/navigation.cpp
    src/navigation/orbit_determination.cpp
)

target_include_directories(gnc_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/shared
    ${PROJECT_SOURCE_DIR}/sim/include
    ${PROJECT_SOURCE_DIR}/shared/external/eigen-5.0.0/eigen-5.0.0
)

target_compile_definitions(gnc_lib PRIVATE EIGEN_NO_MALLOC)

set(CMAKE_CXX_FLAGS "-O2 -Wall -Wextra -Wpedantic")
