# How to compile SOLERA in Polaris (ALCF)

Developed in collaboration with Hewlett Packard Enterprise (HPE), Polaris is a leading-edge system that will give scientists and application developers a platform to test and optimize codes for Aurora, Argonne's upcoming Intel-HPE exascale supercomputer.

The Polaris software environment is equipped with the HPE Cray programming environment, HPE Performance Cluster Manager (HPCM) system software, and the ability to test programming models, such as OpenMP and SYCL, that will be available on Aurora and the next-generation of DOE’s high performance computing systems. Polaris users will also benefit from NVIDIA’s HPC software development kit, a suite of compilers, libraries, and tools for GPU code development. 

## Requirements

In order to compile SOLERA, first you need to install the following libraries:
- PETSc
- Eigen

## Compilation instructions

Load CMake:
```console
module use /soft/modulefiles
module load spack-pe-base cmake
```

Load GNU compilers:
```console
module swap PrgEnv-nvhpc PrgEnv-gnu
module load nvhpc-mixed
```

Execute configure 
```console
chmod +x configure.sh
./configure.sh
```

## Run SOLERA test 
