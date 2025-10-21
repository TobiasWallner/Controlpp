Controlpp
=========

[![gcc-build-and-test](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_gcc.yml/badge.svg)](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_gcc.yml)

[![clang-build-and-test](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_clang.yml/badge.svg)](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_clang.yml)

A C++ library for control systems: from classic transfer-function workflows to state-space, estimators, and optimal control. Design continuous or discrete controllers, run time-variant logic, and bring in LQR, Kalman filtering, and even H₂ synthesis building blocks.

Features:
---------

- **Performance**
  Supports [Vectorisation](https://eigen.tuxfamily.org/index.php?title=FAQ#Vectorization) for: SSE2, SSE3, SSE4, AVX, AVX2, AVX512, AltiVec/VSX, ARM NEON and now S390x SIMD (ZVector) through the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library.
- **Classical Control, Transfer Functions and State Space**  
  Work directly with continuous and discrete transfer functions, numerator/denominator polynomials, zeros/poles, and conversions from and to state space realisations.
  Switch domains as needed with zero order hold or tustin in a type save manner:  
  s -- zero-order hold --> z -- inverser-bilinear (tustin) --> q -- bilinear (tustin) --> z
- **Time-Variant Controllers**  
  Have high jitter or changing sampling rates? E.g.: Controlling with a sensor (e.g.: camera) that changes its and sampling rate (ROI)?
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
CPMAddPackage("gh:TobiasWallner/Controlpp#main")

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
    const auto s controlpp::tf::s<float>;

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

---- 

### Known Issues:
- `-W-array-bounds`
  might get triggered by Eigen's use of SIMD instructions. The operations are all safe and heavily tested by Eigen. It might occur when for example a vector with a size of 3 is being used with a SIMD instruction that uses 4 values. To get rid of the warning, use: `-Wno-array-bounds`
- `-W-maybe-uninitialized`
  Might get triggered by Eigen's solvers if the compiler fails to proove that an uninitialised value is not used. To get rid of the warning, use: `-Wno-maybe-uninitialized`

----

### Set your CMake default compiler

#### On Linux:
Add to `~/.bashrc`, `~/.zshrc` or the config script of your terminal:
```
export CC=/path/to/gcc
export CXX=/path/to/g++
```

likely path: `/usr/bin/`.

#### On Windows:
add the environment variables
```
setx CC "\path\to\gcc.exe"
setx CXX "\path\to\g++.exe"
```

likely path: `C:\mingw64\bin\` or `C:\Program Files\LLVM\bin\`.

----

###  Set your CMake default generator

#### On Linux:
Add to `~/.bashrc`, `~/.zshrc` or the config script of your terminal:
```
export CMAKE_GENERATOR="Ninja Multi-Config"
```

#### On Windows:
Add the environment variable
```
setx CMAKE_GENERATOR "Ninja Multi-Config"
```

----

### Set your default CPM Package Manager library cache

#### On Linux:
Add to `~/.bashrc`, `~/.zshrc` or the config script of your terminal:
```
export CPM_SOURCE_CACHE=$HOME/.cache/CPM
```

#### On Windows:
Add the environment variable
```
setx CPM_SOURCE_CACHE "%USERPROFILE%\.cache\CPM"
```