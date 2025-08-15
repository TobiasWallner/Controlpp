config = Release
CC = gcc
CXX = g++

build/CMakeChache.txt: CMakeLists.txt
	cmake -S . -B build -G "Ninja Multi-Config" -DCONTROLPP_COMPILE_TEST=ON -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX)

config: build/CMakeChache.txt


.PHONY: build
build: config
	cmake --build build --config $(config)

.PHONY: test
test: build
	ctest --test-dir build --build-config $(config) -V

.PHONY: clean
clean:
	cmake --build build --target clean --config Release
	cmake --build build --target clean --config Debug