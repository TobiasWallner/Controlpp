config = Release

build/$(config)/config.timestampfile: CMakeLists.txt
	cmake -S . -B build/$(config)/ -DCMAKE_BUILD_TYPE=$(config) -DCONTROLPP_COMPILE_TEST=ON
	echo "" > $@


.PHONY: build
build: build/$(config)/config.timestampfile
	cmake --build build/$(config)/ --config $(config)

.PHONY: test
test: build
	ctest --test-dir build/$(config)/ --build-config $(config) --output-on-failure

.PHONY: clean
clean:
	cmake --build build/Release/ --target clean --config Release
	cmake --build build/Debug/ --target clean --config Debug