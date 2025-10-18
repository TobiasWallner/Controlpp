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

target_include_directories(controlpp SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})

target_link_libraries(controlpp INTERFACE
    Eigen3::Eigen
)

target_compile_features(controlpp INTERFACE cxx_std_23)

##################################################################################
#                               TESTING
##################################################################################

if (CONTROLPP_COMPILE_TEST)

    CPMAddPackage("gh:google/googletest@1.17.0")

    enable_testing()

    add_executable(controlpp_tests
        tests/main.cpp

        #tests/H2Controller_Tests.cpp
        tests/LQRController_Tests.cpp
        
        tests/Pade_Tests.cpp
        tests/TimeVariantControl_Test.cpp
        tests/KalmanFilter_Tests.cpp
        tests/Estimators_Tests.cpp
        tests/Transformation_Tests.cpp
        tests/StateSpace_Tests.cpp
        tests/Polynom_Tests.cpp
        tests/math_Tests.cpp
        tests/analysis_Tests.cpp
        tests/Generator_Tests.cpp

    )

    set_target_properties(controlpp_tests PROPERTIES LINKER_LANGUAGE CXX)

    target_link_libraries(controlpp_tests
        controlpp
        GTest::gtest_main
    )

    if ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR (CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
        target_compile_options(controlpp_tests PUBLIC
            -Wall
            -Wextra
            -Werror
            -Wno-array-bounds # Gets triggered by SIMD
            -Wno-maybe-uninitialized
            #-fmax-errors=1
        )
    endif()

    include(GoogleTest)
    add_test(NAME controlpp_tests COMMAND controlpp_tests)

endif()


##################################################################################
#                               Examples
##################################################################################

