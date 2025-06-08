config = Release

CC = gcc
CXX = g++

.PHONY: build
test:
	cmake -S . -B build -G "Ninja Multi-Config" -DCONTROLPP_COMPILE_TEST=ON -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX)
	cmake --build build --config Release
	ctest --test-dir build -C $(config) -V