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

	include(FetchContent)
	FetchContent_Declare(
			googletest
			URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip # v1.12.1 (last version to support C++11)
	)
	# For Windows: Prevent overriding the parent project's compiler/linker settings
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	FetchContent_MakeAvailable(googletest)

	enable_testing()

	add_subdirectory(test)

endif()