cmake_minimum_required(VERSION 3.8)

project(sample_statistics LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

##
## MAIN_PROJECT CHECK
## determine if sample_statistics is built as a subproject (using add_subdirectory) or if it is the main project
##
set(MAIN_PROJECT OFF)
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  set(MAIN_PROJECT ON)
endif()


##
## TARGET
## create target and add include path
##
add_library(sample_statistics INTERFACE)
if (${CMAKE_VERSION} VERSION_LESS "3.8.0")
    target_compile_features(sample_statistics INTERFACE cxx_range_for)
else()
    target_compile_features(sample_statistics INTERFACE cxx_std_11)
endif()

target_include_directories(
    sample_statistics
    INTERFACE
    include
)

add_subdirectory(examples)


# Testing with GoogleTest (optional, but recommended)
option (INCLUDE_TESTS 
        "Unit testing with GoogleTest" ON) 

if(INCLUDE_TESTS)

	# Download and unpack googletest at configure time
	configure_file(CMakeLists.txt.gtest ${PROJECT_BINARY_DIR}/googletest-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "Build step for googletest failed: ${result}")
	endif()

	# Prevent overriding the parent project's compiler/linker
	# settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

	# Add googletest directly to our build. This defines
	# the gtest and gtest_main targets.
	add_subdirectory(${PROJECT_BINARY_DIR}/googletest-src
					 ${PROJECT_BINARY_DIR}/googletest-build
					 EXCLUDE_FROM_ALL)

	# The gtest/gtest_main targets carry header search path
	# dependencies automatically when using CMake 2.8.11 or
	# later. Otherwise we have to add them here ourselves.
	if (CMAKE_VERSION VERSION_LESS 2.8.11)
	  include_directories("${gtest_SOURCE_DIR}/include")
	endif()

	enable_testing()

	add_subdirectory(test)

endif()