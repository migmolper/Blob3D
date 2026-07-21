#  .$$$$$$\...$$$$$$\..$$\.......$$$$$$$$\.$$$$$$$\...$$$$$$\..
#  $$..__$$\.$$..__$$\.$$.|......$$.._____|$$..__$$\.$$..__$$\.
#  $$./..\__|$$./..$$.|$$.|......$$.|......$$.|..$$.|$$./..$$.|
#  \$$$$$$\..$$.|..$$.|$$.|......$$$$$\....$$$$$$$..|$$$$$$$$.|
#  .\____$$\.$$.|..$$.|$$.|......$$..__|...$$..__$$<.$$..__$$.|
#  $$\...$$.|$$.|..$$.|$$.|......$$.|......$$.|..$$.|$$.|..$$.|
#  \$$$$$$..|.$$$$$$..|$$$$$$$$\.$$$$$$$$\.$$.|..$$.|$$.|..$$.|
#  .\______/..\______/.\________|\________|\__|..\__|\__|..\__|
#  ............................................................
#  ............................................................
#  ............................................................

cmake_minimum_required(VERSION 3.19)

project(BLOB3D VERSION 24.12 LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add MPI (Message Passing Interface)
option(USE_MPI "Use MPI library" ON)

# Add VTK (Visualization Toolkit)
option(USE_VTK "Use VTK library" ON)

# Add HDF5
option(USE_HDF5 "Use HDF5 library" ON)

# Add OpenMP (Open Multi-Processing)
option(USE_OPENMP "Use OpenMP library" OFF)

# Add SLEPc (Scalable Library for Eigenvalue Problem Computations)
option(USE_SLEPc "Use SLEPc library" OFF)

# Build tests
option(BLOB3D_BUILD_TESTING "Build BLOB3D tests" ON)

# Include coverage
option(ENABLE_COVERAGE "Enable coverage utilities" OFF)
if(ENABLE_COVERAGE)
  set(CMAKE_MODULE_PATH
      "${CMAKE_SOURCE_DIR}/external/CMake-codecov/cmake" ${CMAKE_MODULE_PATH})
endif()

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)

set(EIGEN_DIR "${CMAKE_SOURCE_DIR}/external/eigen3" CACHE PATH
    "Path to Eigen headers")
set(BLOB3D_SPLINES_DIR_DEFAULT "${CMAKE_SOURCE_DIR}/src/Potentials/splines"
    CACHE PATH "Default folder containing BLOB3D spline potential sets")

message(STATUS "****************************")
message(STATUS "* ░█▀▄░█░░░█▀█░█▀▄░▀▀█░█▀▄ *")
message(STATUS "* ░█▀▄░█░░░█░█░█▀▄░░▀█░█░█ *")
message(STATUS "* ░▀▀░░▀▀▀░▀▀▀░▀▀░░▀▀░░▀▀░ *")
message(STATUS "****************************")
message(STATUS "Version: ${PROJECT_VERSION}")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wno-undefined-var-template")
endif()

set(USED_CMAKE_GENERATOR "${CMAKE_GENERATOR}" CACHE STRING
    "Expose CMAKE_GENERATOR" FORCE)
if(USED_CMAKE_GENERATOR MATCHES "Ninja")
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
endif()

if(APPLE)
  message("MACOSX system")
  add_compile_definitions(SYS_MACOSX)
  set(BLA_VENDOR Apple)
elseif(UNIX)
  message("UNIX system")
  add_compile_definitions(SYS_LINUX)
  set(BLA_VENDOR OpenBLAS)
else()
  message(FATAL_ERROR "Your system has not been recongnised.")
endif()

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_BUILD_TYPE MATCHES "Debug")
  message(STATUS "Debug mode")
  add_compile_definitions(DEBUG_MODE)
  set(CMAKE_CXX_FLAGS_DEBUG
      "--coverage -O0 -g -Wall -Wpedantic -Wextra -Wunused-variable")
elseif(CMAKE_BUILD_TYPE MATCHES "Release")
  message(STATUS "Release mode")
endif()

find_package(PkgConfig REQUIRED)
pkg_search_module(PETSC REQUIRED IMPORTED_TARGET PETSc)
if(PETSC_FOUND)
  message(STATUS "PETSC version linked: ${PETSC_VERSION}")
endif()

set(BLOB3D_LINK_LIBS m PkgConfig::PETSC)

if(USE_SLEPc)
  pkg_search_module(SLEPC REQUIRED IMPORTED_TARGET SLEPc)
  if(SLEPC_FOUND)
    message(STATUS "SLEPc version linked: ${SLEPC_VERSION}")
    add_compile_definitions(USE_SLEPC)
    list(APPEND BLOB3D_LINK_LIBS PkgConfig::SLEPC)
  endif()
endif()

if(USE_MPI)
  find_package(MPI REQUIRED)
  if(MPI_FOUND)
    add_compile_definitions(USE_MPI)
    list(APPEND BLOB3D_LINK_LIBS MPI::MPI_CXX)
  endif()
endif()

if(USE_OPENMP)
  find_package(OpenMP REQUIRED)
  if(OPENMP_FOUND)
    add_compile_definitions(USE_OPENMP)
    list(APPEND BLOB3D_LINK_LIBS OpenMP::OpenMP_CXX)
  endif()
endif()

if(USE_HDF5)
  add_compile_definitions(USE_HDF5)
endif()

if(USE_VTK)
  find_package(VTK REQUIRED)
  if(VTK_FOUND)
    message(STATUS "VTK version linked: ${VTK_VERSION}")
    add_compile_definitions(USE_VTK)
    list(APPEND BLOB3D_LINK_LIBS ${VTK_LIBRARIES})
  endif()
endif()

add_subdirectory(src)

if(BLOB3D_BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

# pkg-config file for downstream: pkg_check_modules(BLOB3D IMPORTED_TARGET blob3D)
set(BLOB3D_PC_REQUIRES "PETSc")
if(USE_SLEPc)
  set(BLOB3D_PC_REQUIRES "${BLOB3D_PC_REQUIRES} SLEPc")
endif()
set(BLOB3D_PC_CFLAGS "")
if(EIGEN_DIR)
  set(BLOB3D_PC_CFLAGS " -I${EIGEN_DIR}")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/blob3D.pc.in
  ${CMAKE_CURRENT_BINARY_DIR}/blob3D.pc
  @ONLY
)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BLOB3DConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/BLOB3DConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BLOB3D
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/BLOB3DConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(EXPORT BLOB3DTargets
  FILE BLOB3DTargets.cmake
  NAMESPACE BLOB3D::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BLOB3D
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/BLOB3DConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/BLOB3DConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BLOB3D
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/blob3D.pc
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
