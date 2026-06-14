cmake_minimum_required(VERSION 3.16)

project(controlpp LANGUAGES CXX)

# ---------------------------------------------------------------------------
#                       Add Libraries
# ---------------------------------------------------------------------------

include(cmake/BPM.cmake)

# Add eigen for the main matrix/vector algebraics
BPMAddInstallPackage("https://gitlab.com/libeigen/eigen#^5.0.0" PACKAGES Eigen3)

# Add blas and lapack for Dare solvers
find_package(BLAS)
if(BLAS_FOUND)
    message(STATUS "BLAS found: ${BLAS_LIBRARIES}")
else()
    message(STATUS "BLAS not found")
endif()

find_package(LAPACK)
if(LAPACK_FOUND)
    message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
else()
    message(STATUS "LAPACK not found")
endif()

if(BLAS_FOUND AND LAPACK_FOUND)
    message(STATUS "BALS and LAPACK will be used as a DARE solver")
else()
    message(STATUS "A fallback DARE solver will be used since BLAS and LAPACK have not been found")
endif()

# Add a csv reader to read measurement data
BPMAddSourcePackage("https://github.com/TobiasWallner/CSVd#v4.0.2")

BPMAddInstallPackage("https://github.com/google/googletest#1.17.0" PACKAGES GTest PRIVATE)

BPMMakeAvailable()

# ---------------------------------------------------------------------------
#                       Define this library
# ---------------------------------------------------------------------------

add_library(controlpp STATIC
    src/Bode.cpp
)

# target_compile_features(controlpp PRIVATE cxx_std_20)

target_include_directories(controlpp PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}/include
)

target_include_directories(controlpp SYSTEM PRIVATE ${EIGEN3_INCLUDE_DIR})

target_link_libraries(controlpp PUBLIC
    Eigen3::Eigen
    csvd
)

target_compile_features(controlpp PUBLIC cxx_std_23)

##################################################################################
#                               TESTING
##################################################################################

if (CONTROLPP_COMPILE_TEST)

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
        tests/algorithm_Tests.cpp
        tests/analysis_Tests.cpp
        tests/Generator_Tests.cpp
        tests/Bode_Tests.cpp

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
        )
    endif()

    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        target_compile_options(controlpp_tests PUBLIC
            -Wno-maybe-uninitialized
            -Wno-array-bounds # Gets triggered by SIMD
            -Wno-maybe-uninitialized
            -fmax-errors=3
        )
    endif()

    include(GoogleTest)
    add_test(NAME controlpp_tests COMMAND controlpp_tests)

endif()


##################################################################################
#                               Examples
##################################################################################

