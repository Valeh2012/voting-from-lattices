#cmake version
cmake_minimum_required(VERSION 2.8)

#project name
project(CumulativeVoting-BP)

#flags
if(APPLE)
  cmake_policy(SET CMP0042 NEW)
endif()

option(NFL_OPTIMIZED "Build an optimized version of NFL" OFF)
if(NFL_OPTIMIZED)
    message(STATUS "Tring to use optimized version of NFL")
    add_definitions(-DNFL_OPTIMIZED)

    include(CheckCXXSourceCompiles)
    set(CMAKE_CXX_FLAGS_OLD ${CMAKE_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

    check_cxx_source_compiles("#include <immintrin.h>\nint main() { return 0;}" NFLLIB_USE_AVX)
    check_cxx_source_compiles("#include <pmmintrin.h>\nint main() { return 0;}" NFLLIB_USE_SSE)
    if(NFLLIB_USE_AVX)
        message(STATUS "Using AVX vector engine")
        add_definitions(-DNTT_AVX2)
    elseif(NFLLIB_USE_SSE)
        message(STATUS "Using SSE vector engine")
        add_definitions(-DNTT_SSE)
    else()
        message(ERROR "Asking for optimized mode but neither SSE not AVX support detected")
    endif()
endif()


#files
file(GLOB SRC_FILES ./lib/*.c ./lib/*.cpp ./src/vericrypt.cpp ./main.cpp)

#include
include(CMakeCompilers.txt)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
#find_package(NFL REQUIRED)
include_directories(include)
include_directories(lib)

#target
add_executable(vericrypt ${SRC_FILES} )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl, stack-size=10000000")

set(NFLLIB "./NFLlib") # change it to NFLlib download path
add_library(nfllib SHARED IMPORTED)
set_target_properties(nfllib PROPERTIES IMPORTED_LOCATION "../deps/libnfllib.so") #change it to NFLlib installation path
target_link_libraries(vericrypt nfllib ${GMP_LIBRARY} ${MPFR_LIBRARY})
target_include_directories(vericrypt PUBLIC "${NFLLIB}/include/")

