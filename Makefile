config = Release


build/$(config)/config.timestampfile: CMakeLists.txt
	cmake -S . -B build/$(config)/ -G $(G) -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) -DCMAKE_BUILD_TYPE=$(config) -DCONTROLPP_COMPILE_TEST=ON
	echo "" > $@


.PHONY: build
build: build/$(config)/config.timestampfile
	cmake --build build/$(config)/ --config $(config)

.PHONY: test
test: build
	ctest --test-dir build/$(config)/ --build-config $(config) -V

.PHONY: clean
clean:
	cmake --build build/Release/ --target clean --config Release
	cmake --build build/Debug/ --target clean --config Debug