Controlpp
=========

[![gcc-build-and-test](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_gcc.yml/badge.svg)](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_gcc.yml)

[![clang-build-and-test](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_clang.yml/badge.svg)](https://github.com/TobiasWallner/Controlpp/actions/workflows/test_with_clang.yml)

A systems control library for C++.

Usage:
------

```
make config
```

```
make build
```

```
make test
```



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