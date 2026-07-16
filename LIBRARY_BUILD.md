# Building and using the SOLERA library

A copy of this guide also lives in the documentation tree: [docs/LIBRARY_BUILD.md](docs/LIBRARY_BUILD.md).

This document describes how to configure, build, and install **SOLERA** as a **static library** (`libsolera`), and how to link it from a separate C++ project.

## Install with Spack

For a reproducible install across systems using [Spack](https://spack.readthedocs.io/), see **[SPACK.md](SPACK.md)** (also [docs/SPACK.md](docs/SPACK.md)).

## Prerequisites

- **CMake** 3.19 or newer  
- **C++14** compiler  
- **PETSc** (must be discoverable via **pkg-config**, e.g. `PETSc.pc` on `PKG_CONFIG_PATH`)  
- **Eigen** headers (`EIGEN_DIR` pointing to the directory that contains the `Eigen/` folder)  
- Optional, depending on CMake options: **MPI**, **VTK**, **HDF5**, **OpenMP**, **SLEPc**, **CGAL** (must match how you built SOLERA)

Ensure `PETSC_DIR`, `PETSC_ARCH`, and `PKG_CONFIG_PATH` include your PETSc (and SLEPc, if used) `pkgconfig` directories before running CMake.

## Build and install (typical)

From the **repository root** (the directory that contains the top-level `CMakeLists.txt`):

```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$HOME/opt/solera-24.12" \
  -DSOLERA_BUILD_TESTING=OFF \
  -DEIGEN_DIR="/path/to/eigen" \
  -DCMAKE_CXX_COMPILER=/path/to/mpicxx

cmake --build build -j
cmake --install build
```

Common CMake options (see top-level `CMakeLists.txt`):

| Option | Default | Meaning |
|--------|---------|---------|
| `USE_MPI` | ON | MPI support |
| `USE_VTK` | ON | VTK linking |
| `USE_HDF5` | ON | Compile definitions for HDF5 |
| `USE_OPENMP` | OFF | OpenMP |
| `USE_SLEPc` | OFF | SLEPc |
| `USE_CGAL` | OFF | CGAL |
| `SOLERA_BUILD_TESTING` | ON | Build unit tests |

After installation you should have:

- **Library:** `<prefix>/lib/libsolera.a` (static)
- **Headers:** `<prefix>/include/SOLERA/**/*.hpp`
- **CMake package:** `<prefix>/lib/cmake/SOLERA/` (`SOLERAConfig.cmake`, `SOLERATargets.cmake`, …)
- **pkg-config:** `<prefix>/lib/pkgconfig/solera.pc`

## Environment for downstream projects

### Using CMake `find_package`

Point CMake at the install prefix:

```bash
export CMAKE_PREFIX_PATH="/path/to/install/prefix:${CMAKE_PREFIX_PATH}"
```

### Using pkg-config

Include both **SOLERA** and **PETSc** (and any other `Requires` from `solera.pc`) on `PKG_CONFIG_PATH`:

```bash
export PKG_CONFIG_PATH="/path/to/install/prefix/lib/pkgconfig:${PETSC_DIR}/${PETSC_ARCH}/lib/pkgconfig:${PKG_CONFIG_PATH}"
```

Verify:

```bash
pkg-config --cflags --libs solera
```

## Using SOLERA from another project (recommended: CMake)

Create a separate directory for your application, e.g. `my_app/` with:

**`my_app/CMakeLists.txt`**

```cmake
cmake_minimum_required(VERSION 3.19)
project(MyApp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SOLERA CONFIG REQUIRED)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE SOLERA::solera)
```

**`my_app/main.cpp`** (example skeleton only—replace with real SOLERA APIs you need)

```cpp
// Include headers from the installed tree: <prefix>/include/SOLERA/...
// Example (adjust to the headers your code actually uses):
#include "Variables.hpp"

int main(int argc, char** argv) {
  (void)argc;
  (void)argv;
  // Initialize PETSc/MPI etc. as required by your SOLERA workflow.
  return 0;
}
```

Configure and build your app:

```bash
cmake -S my_app -B my_app/build \
  -DCMAKE_PREFIX_PATH="/path/to/solera/install/prefix"
cmake --build my_app/build -j
```

`SOLERA::solera` carries **include directories** and **linked libraries** (PETSc, MPI, VTK, …) as configured when SOLERA was built, so you normally only need `target_link_libraries(... SOLERA::solera)`.

## Using SOLERA with pkg-config in CMake

If you prefer pkg-config:

```cmake
cmake_minimum_required(VERSION 3.19)
project(MyApp LANGUAGES CXX)

find_package(PkgConfig REQUIRED)
pkg_check_modules(solera REQUIRED IMPORTED_TARGET solera)

add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE PkgConfig::solera)
```

Ensure `PKG_CONFIG_PATH` is set as described above before running CMake.

## Using `add_subdirectory` (same machine, no install)

From another CMake project that can see this repository:

```cmake
add_subdirectory(/path/to/DMD dmd_build EXCLUDE_FROM_ALL)
add_executable(my_app main.cpp)
target_link_libraries(my_app PRIVATE SOLERA::solera)
```

This is convenient for development but couples both trees; **install + `find_package`** is better for redistribution.

## Troubleshooting

- **PETSc not found during SOLERA configure:** set `PKG_CONFIG_PATH` so `pkg-config --modversion PETSc` works.
- **`INTERFACE_INCLUDE_DIRECTORIES` error when installing:** ensure you use a current `SOLERA` `CMakeLists.txt` that uses `$<BUILD_INTERFACE:...>` / `$<INSTALL_INTERFACE:...>` for includes.
- **Missing Eigen in downstream compile:** the installed `SOLERA::solera` target and `solera.pc` are set up to reference Eigen as configured at SOLERA build time; pass `-DEIGEN_DIR=...` when building SOLERA if Eigen is not the default bundled path.

For a scripted workflow matching a local PETSc layout, you can adapt `configure.sh` in this repository.
