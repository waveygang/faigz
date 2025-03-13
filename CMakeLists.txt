cmake_minimum_required(VERSION 3.10)
project(faigz LANGUAGES C CXX)

# Find htslib
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTSLIB REQUIRED htslib)
include_directories(${HTSLIB_INCLUDE_DIRS})
link_directories(${HTSLIB_LIBRARY_DIRS})

# Main library target - header only
add_library(faigz INTERFACE)
target_include_directories(faigz INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# C benchmark executable
add_executable(bench_faigz bench_faigz.c)
target_link_libraries(bench_faigz ${HTSLIB_LIBRARIES} pthread)

# C++ test target
add_executable(test_faigz_cpp test_faigz.cpp)
target_link_libraries(test_faigz_cpp ${HTSLIB_LIBRARIES} pthread)
# Force C++ compilation for this target
set_target_properties(test_faigz_cpp PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES)
