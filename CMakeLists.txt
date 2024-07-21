cmake_minimum_required(VERSION 3.18.0)

set(PROJECT_NAME test)
project(${PROJECT_NAME})

# Basic compiler settings
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)

# project building process: 
add_executable(${PROJECT_NAME}
	test.cpp
)

target_compile_definitions(${PROJECT_NAME} PUBLIC
	
)

target_include_directories(${PROJECT_NAME} PUBLIC
	.
	Matrix
	Vector
)

target_link_libraries(${PROJECT_NAME} PUBLIC

)

# Compiler Options for Clang:
# ===========================
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	set(COMPILER_FLAGS 
		
		# set compiler flags for debugging
		#---------------------------------
		$<$<CONFIG:Debug>:
			-O0	# disable optimisations
			-g	# enable debug information
		>
		
		# set compiler flags for release
		#-------------------------------
		$<$<CONFIG:Release>:
			-O3 # Enable all Optimisations
		>
	)
# Compiler Options for GNU:
# =========================
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(COMPILER_FLAGS
		
		# Enable link time optimisations
		-flto
		
		# set compiler flags for debugging
		#---------------------------------
		$<$<CONFIG:Debug>:
			-O0	# disable optimisations
			-g	# enable debug information
			-ggdb # extra flag for the debugger
		>
		
		# set compiler flags for release
		#-------------------------------
		$<$<CONFIG:Release>:
			-O3 # Enable all Optimisations
		>
	)
# Compiler Options for MSVC:
# ==========================
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	set(COMPILER_FLAGS
		/W3
		/Zc:__cplusplus # enable __cplusplus
		# set compiler flags for debugging
		#---------------------------------
		$<$<CONFIG:Debug>:
		
			# disable optimisations for debugging
			/Od 
			
			# produces a PDB file in a format that supports the Edit and Continue feature. 
			# To use Edit and Continue debugging features, you must use this option. 
			/ZI 
		>
		
		# set compiler flags for release versions
		#----------------------------------------
		$<$<CONFIG:Release>:
			
			# Creates fast code
			/O2 
		>
	)
else()
	message( SEND_ERROR "The compiler you are using is not handeled. Add your compiler and compiler flags to the if(CMAKE_CXX_COMPILER_ID STREQUAL <your compiler>) in the CMakeLists.txt of the root folder")
endif()


target_compile_options(${PROJECT_NAME} PUBLIC
	${COMPILER_FLAGS}
)

