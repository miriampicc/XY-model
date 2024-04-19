
cmake_minimum_required(VERSION 3.24)

### Setup project name and language
project(CMT VERSION 1.1
        DESCRIPTION "A CMake project template"
        HOMEPAGE_URL "https://github.com/DavidAce/CMakeTemplate"
        LANGUAGES CXX
        )

set(CMAKE_OSX_SYSROOT /Library/Developer/CommandLineTools/SDKs/MacOSX14.sdk CACHE PATH "The path to the macOS SDK to use")


### Set options
option(CMT_ENABLE_OPENMP "Enable OpenMP flags such as -fopenmp"   OFF)
option(CMT_ENABLE_MPI    "Enables use of MPI (work in progress)"  OFF)
option(CMT_ENABLE_TESTS  "Enables CTest testing"                  OFF)
option(CMT_CMAKE_DEBUG   "Print info during CMake configuration"  OFF)


################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git identifiers     ###
################################################################
include(cmake/gitversion.cmake)


if(CMT_CMAKE_DEBUG)
    ### Print operating system details
    include(cmake/PrintHostInfo.cmake)
endif()

### Add all source files
add_executable(CMT
        source/main.cpp
        source/measures.cpp
        source/rng.cpp
        source/memory_check.h
        source/o2.h
        source/robust_filesystem.h
        source/initialization.cpp
        source/initialization.h
        source/02.h
        source/montecarlo.h
        source/montecarlo.cpp
        source/parallel_temp.cpp
        source/parallel_temp.h
)
target_include_directories(CMT PRIVATE source)
set_target_properties(CMT PROPERTIES LINK_WHAT_YOU_USE TRUE)


### Create a helper targets that can be reused to build tests
add_library(cmt-libs INTERFACE)
add_library(cmt-flags INTERFACE)
#include(cmake/CompilerFlags.cmake)



### Find dependencies
#find_package(h5pp   1.11.0 REQUIRED)
#find_package(Eigen3 3.4.0  REQUIRED)
#target_link_libraries(cmt-libs INTERFACE h5pp::h5pp Eigen3::Eigen)


if(CMT_ENABLE_OPENMP)
    # Try to find the package
    find_package(OpenMP REQUIRED COMPONENTS CXX)
    target_link_libraries(cmt-flags INTERFACE OpenMP::OpenMP_CXX)

#    if(OpenMP_CXX_FOUND)
#        message(STATUS "OpenMP found: ${OpenMP_CXX_FLAGS}")
#        target_compile_options(cmt-flags INTERFACE ${OpenMP_CXX_FLAGS})
#        target_link_libraries(cmt-flags INTERFACE OpenMP::OpenMP_CXX)
#    else()
#        message(STATUS "OpenMP not found, trying manual specification...")
#        # Set the paths manually
#        set(OpenMP_BASE_DIR "/opt/homebrew/opt/libomp")
#
#        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${OpenMP_BASE_DIR}/include")
#        set(OPENMP_LIBRARY "${OpenMP_BASE_DIR}/lib/libomp.dylib")
#
#        target_compile_options(cmt-flags INTERFACE ${OpenMP_CXX_FLAGS})
#        target_link_libraries(cmt-flags INTERFACE ${OPENMP_LIBRARY})
#
#        include_directories("${OpenMP_BASE_DIR}/include")
#
#    endif()
endif()

if(CMT_ENABLE_MPI)
    find_package(MPI COMPONENTS CXX REQUIRED)
    target_link_libraries(cmt-flags INTERFACE MPI::MPI_CXX)
    #include(cmake/PrintTargetProperties.cmake)
    #print_target_properties(MPI::MPI_CXX)
endif()

if(APPLE)
    # Assume default Homebrew installation path
    #include_directories(/usr/local/opt/libomp/include)
    #link_directories(/usr/local/opt/libomp/lib)
endif()



### Link to executable
target_link_libraries(CMT PUBLIC cmt-libs cmt-flags)  # Link everything to main project


### CTest
if(CMT_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()


if(CMT_CMAKE_DEBUG)
    # Print summary of CMake configuration
    include(cmake/PrintTargetInfo.cmake)
    print_and_write_project_summary(CMT)
endif()
