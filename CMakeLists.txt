cmake_minimum_required (VERSION 3.0)
project(geohash)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_SHARED_MODULE_CREATE_C_FLAGS
        "${CMAKE_SHARED_MODULE_CREATE_C_FLAGS} -undefined dynamic_lookup")
    cmake_policy(SET CMP0042 NEW)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -Wall -pedantic -std=c99")

add_library(geohash OBJECT
    geohash.c
    geohashHelper.c)

add_library(geohash-shared  MODULE $<TARGET_OBJECTS:geohash>)
add_library(geohash-static  STATIC $<TARGET_OBJECTS:geohash>)
add_library(geohash-library SHARED $<TARGET_OBJECTS:geohash>)

# Tell CMake we want a .so and not .dylib
set_target_properties(geohash-shared PROPERTIES SUFFIX .so)

# Create a geohashTest file including main
add_executable(geohashTest geohashMain.c)
target_link_libraries(geohashTest geohash-static)

# Compile with cmake -DBMI2=ON ..
# Check for BMI2 so we can use assembly intrinsics for geohash
include(CheckCSourceRuns)
set(CMAKE_REQUIRED_FLAGS "-mbmi2")
check_c_source_runs("
#include <immintrin.h>
int main() {
    long long a = _pdep_u64(0xFFFF, 0x3333);
    _pext_u64(a, 0x5555);
    return 0;
}"
CC_HAS_BMI2_INTRINSICS)

option(BMI2 "Enable compile-time BMI2 support." OFF)
if(CC_HAS_BMI2_INTRINSICS)
    set(BMI2 ON)
else()
    # If you request BMI2, but can't run BMI2 code, we disable BMI2.
    set(BMI2 OFF)
endif()


if(BMI2)
    set_target_properties(geohash PROPERTIES COMPILE_FLAGS "-mbmi2")
endif()
