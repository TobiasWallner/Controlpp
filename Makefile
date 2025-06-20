config = Release
CC = gcc
CXX = g++

build/config.timestamp: CMakeLists.txt
	cmake -S . -B build -G "Ninja Multi-Config" -DCONTROLPP_COMPILE_TEST=ON -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX)
	$(shell echo > $@)
config: build/config.timestamp


.PHONY: build
build: config
	cmake --build build --config $(config)

.PHONY: test
test: build
	ctest --test-dir build --build-config $(config) -V