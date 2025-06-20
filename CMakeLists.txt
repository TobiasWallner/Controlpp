cmake_minimum_required(VERSION 3.16)

include(cmake/CPM.cmake)

CPMAddPackage(
    NAME Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen/
    GIT_TAG 3.4.0
    OPTIONS
      "EIGEN_BUILD_TESTING OFF"
      "BUILD_TESTING OFF"
)

project(Controlpp LANGUAGES CXX)

add_library(controlpp INTERFACE)

# target_compile_features(controlpp PRIVATE cxx_std_20)

target_include_directories(controlpp INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/include
)

target_link_libraries(controlpp INTERFACE
    Eigen3::Eigen
)


if (CONTROLPP_COMPILE_TEST)

    CPMAddPackage("gh:google/googletest@1.17.0")

    enable_testing()

    add_executable(controlpp_tests)

    include(tests/sources.cmake)

    target_compile_features(controlpp_tests PRIVATE cxx_std_23)

    target_link_libraries(controlpp_tests
        controlpp
        GTest::gtest_main
    )

    include(GoogleTest)
    add_test(NAME controlpp_tests COMMAND controlpp_tests)

endif()
