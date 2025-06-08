cmake_minimum_required(VERSION 3.16)

include(cmake/CPM.cmake)

CPMAddPackage(
    NAME eigen
    https://gitlab.com/libeigen/eigen/
)

project(Controlpp LANGUAGES CXX)

add_library(controlpp)

add_executable(controlpp_tests

)
