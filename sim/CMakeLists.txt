add_executable(sim
    main.cpp
)

configure_file(config.json config.json COPYONLY)

set(EIGEN_INCLUDE_DIR
    ${PROJECT_SOURCE_DIR}/shared/external/eigen-5.0.0/eigen-5.0.0
)

if(NOT EXISTS "${EIGEN_INCLUDE_DIR}/Eigen/Dense")
    message(FATAL_ERROR "Vendored Eigen headers were not found: ${EIGEN_INCLUDE_DIR}")
endif()

target_include_directories(sim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/shared
    ${EIGEN_INCLUDE_DIR}
)

if(BUILD_TESTING)
    add_executable(sim_math_tests
        tests/math_tests.cpp
    )
    target_include_directories(sim_math_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/shared
        ${EIGEN_INCLUDE_DIR}
    )
    add_test(NAME sim_math_tests COMMAND sim_math_tests)
endif()

find_package(jsoncpp CONFIG QUIET)

if(TARGET JsonCpp::JsonCpp)
    target_link_libraries(sim PRIVATE JsonCpp::JsonCpp)
elseif(TARGET jsoncpp_lib)
    target_link_libraries(sim PRIVATE jsoncpp_lib)
elseif(TARGET jsoncpp_static)
    target_link_libraries(sim PRIVATE jsoncpp_static)
else()
    find_path(JSONCPP_INCLUDE_DIR NAMES json/value.h PATH_SUFFIXES jsoncpp)
    find_library(JSONCPP_LIBRARY NAMES jsoncpp)

    if(NOT JSONCPP_INCLUDE_DIR OR NOT JSONCPP_LIBRARY)
        message(FATAL_ERROR "JsonCpp was not found. Install its development headers and library.")
    endif()

    add_library(JsonCpp::JsonCpp UNKNOWN IMPORTED)
    set_target_properties(JsonCpp::JsonCpp PROPERTIES
        IMPORTED_LOCATION "${JSONCPP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${JSONCPP_INCLUDE_DIR}"
    )
    target_link_libraries(sim PRIVATE JsonCpp::JsonCpp)
endif()