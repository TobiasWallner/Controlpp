Controlpp
=========

[![gcc-build-and-test](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_gcc.yml/badge.svg)](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_gcc.yml)

[![clang-build-and-test](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_clang.yml/badge.svg)](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_clang.yml)

A C++ library for control systems: from classic transfer-function workflows to state-space, estimators, and optimal control. Design continuous or discrete controllers, run time-variant logic, and bring in LQR, Kalman filtering, and even H₂ synthesis building blocks.

- Documentation: https://tobiaswallner.github.io/Controlpp/
- GitHub Repository: https://github.com/TobiasWallner/Controlpp

Features:
---------

- **Performance**
  Supports Vectorisation for: SSE2, SSE3, SSE4, AVX, AVX2, AVX512, AltiVec/VSX, ARM NEON and now S390x SIMD (ZVector) through the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library.
- **Classical Control, Transfer Functions and State Space**  
  Work directly with continuous and discrete transfer functions, numerator/denominator polynomials, zeros/poles, and conversions from and to state space realisations.
- **Time-Variant Controllers**  
  Have high jitter or changing sampling rates? E.g.: Controlling with a sensor (e.g.: camera) that changes its frame rate?
  The time variant controllers offer you two inputs: the value and the sample-time.
- **Modern Controllers & Optimal Control**  
  Includes linear-quadratic controller support and estimation tooling (LQR / Kalman), with room to extend into H2/H-inf synthesis.
- **Identification & Estimation**
  Recursive Least Squares and Kalman filters to estimate parameters and states, ready to slot into your controller pipelines.


Quick start
----

#### CMake and CPM.cmake

Read more about [CMake](https://cmake.org/) (build file generator) [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake) (package manager) here.

```cmake
# include the package manager
include(cmake/CPM.cmake)

# add the package
CPMAddPackage("gh:TobiasWallner/Controlpp#v1.5.0")

# create your project executable
add_executable(my_app src/main.cpp)

# link the library
target_link_libraries(my_app PRIVATE controlpp)
```

Example:
--------

```cpp
#include <controlpp.hpp>

int main(){
    const auto s = controlpp::tf::s<float>;

    // define a PT2 element:
    const float f = 1000;
    const float omega = 2 * 3.1415 * f;
    const float D = 0.3;
    const auto PT2_s = 
                               1 
      / //-------------------------------------------------
      (1 + (2 * D * s / omega) + (s * s) / (omega * omega));

    // convert from transfer function to state space
    const auto PT2_ss = to_state_space(PT2_s);

    // convert from s (continuous) to z (discrete) domain
    const float Ts = 1.0/10'000.0;
    const auto PT2_z = discretise_zoh(PT2_ss, Ts);
}

```
