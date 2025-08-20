config = Release
CC = gcc
CXX = g++
G = "Ninja"

.PHONY: config
config:
	cmake -S . -B build/$(config)/ -G $(G) -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_BUILD_TYPE=$(config) -DCONTROLPP_COMPILE_TEST=ON


.PHONY: build
build: config
	cmake --build build/$(config)/ --config $(config) --verbose

.PHONY: test
test: build
	ctest --test-dir build/$(config)/ --build-config $(config) -V

.PHONY: clean
clean:
	cmake --build build/Release/ --target clean --config Release
	cmake --build build/Debug/ --target clean --config Debug