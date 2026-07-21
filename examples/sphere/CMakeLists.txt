cmake_minimum_required(VERSION 3.19)

project("Surface Energy" VERSION 1.0.0 LANGUAGES CXX)
message(STATUS "Version: ${PROJECT_VERSION}")

set(CMAKE_CXX_STANDARD 14)

# Add MPI
option(USE_MPI "Use MPI library" ON)

# Add OpenMP
option(USE_OPENMP "Use OpenMP library" ON)
# Add OpenMP
#option(USE_OPENMP "Use OpenMP library" OFF)

# Add SLEPc (Scalable Library for Eigenvalue Problem Computations)
option(USE_SLEPc "Use SLEPc library" OFF)

# Include coverage
option(ENABLE_COVERAGE "Use Check library" OFF)
if(ENABLE_COVERAGE)
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/CMake-codecov/cmake" ${CMAKE_MODULE_PATH})
endif()

SET(CMAKE_COLOR_MAKEFILE ON)
SET(CMAKE_VERBOSE_MAKEFILE OFF)
set(EXECUTABLE_OUTPUT_PATH "../")

# GNU specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

# Intel specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

# Clang specific settings
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
endif()

# Check if Ninja
set(USED_CMAKE_GENERATOR "${CMAKE_GENERATOR}" CACHE STRING "Expose CMAKE_GENERATOR" FORCE)
if (USED_CMAKE_GENERATOR MATCHES "Ninja")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
endif()

# Check OS
IF (APPLE)
    message("MACOSX system")
    add_definitions("-DSYS_MACOSX")
    SET(BLA_VENDOR Apple)
ELSEIF (UNIX)
    message("UNIX system")
    add_definitions("-DSYS_LINUX")
    SET(BLA_VENDOR OpenBLAS)
ELSE ()
    message(FATAL_ERROR "Your system has not been recongnised." )
ENDIF ()

# Compile settings
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE MATCHES "Debug")
message(STATUS "Debug mode")
add_definitions("-DDEBUG_MODE")
set(CMAKE_CXX_FLAGS_DEBUG "--coverage -O0 -g -Wall -Wpedantic -Wextra -Wunused-variable")
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
message(STATUS "Release mode")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

endif()


###############################################################################
## dependencies ###############################################################
###############################################################################

# Add math library
link_libraries(m)

# Include directories
include_directories(BEFORE ${EIGEN_DIR})

# Add PkgConfig
find_package(PkgConfig REQUIRED)

# Add PETSc
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
if (PETSC_FOUND)
 link_libraries(PkgConfig::PETSC)
 message ("PETSC version linked: ${PETSC_VERSION}")
 if(${PETSC_VERSION} VERSION_GREATER_EQUAL "3.21.1")
    add_definitions("-DPETSC_NEW")
  else()
    add_definitions("-DPETSC_OLD")
  endif()
else (PETSC_FOUND)
    message (FATAL_ERROR "${PROJECT_NAME}: Unable to find the PETSC folder.")
endif(PETSC_FOUND)

# Add SLEPc
if (USE_SLEPc)
  pkg_search_module(SLEPC REQUIRED IMPORTED_TARGET SLEPc)
  if (SLEPC_FOUND)
    link_libraries(PkgConfig::SLEPC)
    message ("SLEPc version linked: ${SLEPC_VERSION}")
    add_definitions("-DUSE_SLEPC")
  else (SLEPC_FOUND)
    message (FATAL_ERROR "${PROJECT_NAME}: Unable to find the SLEPc folder.")
  endif(SLEPC_FOUND)  
endif()

# Add MPI
if(USE_MPI)
  find_package(MPI REQUIRED)
  if (MPI_FOUND)
    add_definitions("-DUSE_MPI")
    include_directories(${MPI_CXX_INCLUDE_DIRS})
    link_libraries(${MPI_CXX_LIBRARIES})
  endif()
endif()

# Add BLOB3D
find_package(PkgConfig REQUIRED)
pkg_check_modules(BLOB3D REQUIRED IMPORTED_TARGET blob3D)

# If blob3D.pc was built elsewhere, Libs: -L... may not match this machine (e.g. lib vs lib64).
if(DEFINED ENV{BLOB3D_DIR})
  foreach(_blob3D_lib IN ITEMS "$ENV{BLOB3D_DIR}/lib64" "$ENV{BLOB3D_DIR}/lib")
    if(EXISTS "${_blob3D_lib}/libblob3D.a")
      link_directories("${_blob3D_lib}")
      break()
    endif()
  endforeach()
endif()

link_libraries(PkgConfig::BLOB3D)

# Add VTK (required by static libblob3D VTK I/O).
find_package(VTK REQUIRED)
if(VTK_FOUND)
  link_libraries(${VTK_LIBRARIES})
  message("VTK version linked: ${VTK_VERSION}")
endif()

# Add OpenMP (required when linking static BLOB3D builds compiled with OpenMP).
if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OpenMP_CXX_FOUND)
    link_libraries(OpenMP::OpenMP_CXX)
    add_definitions("-DUSE_OPENMP")
    message("OpenMP linked")
  endif()
endif()

###############################################################################
## Create executable ##########################################################
###############################################################################

add_executable("exe" main.cpp)
target_compile_definitions("exe" PUBLIC -DMESH_UPDATE=0)
target_compile_definitions("exe" PUBLIC -DCOMPUTE_ANGULAR_TERMS=1)
target_compile_definitions("exe" PUBLIC -DMONOLITHIC_SCHEME=0)
target_compile_definitions("exe" PUBLIC -DSTAGGERED_SCHEME=1)
target_compile_definitions("exe" PUBLIC -DGH3TH_INTEGRAL)
target_compile_definitions("exe" PUBLIC -DNUMERICAL_DERIVATIVES=0)
target_compile_definitions("exe" PUBLIC -DACTIVATE_MIGRATION=0)
