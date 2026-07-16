# How to compile SOLERA in Hercules (CICA)

- Date: 27/02/2025
- Author: M.Molinos

The Scientific Computing Center of Andalucía (CICA) provides the academic and research community of Andalucía with access to the Hércules High-Performance Computing (HPC) cluster. This system consists of 232 computing nodes, totaling 47.5 terabytes of RAM, 11,136 cores, 1 petabyte of storage, and a computing power of 855 teraflops. Node connectivity is established through a low-latency Infiniband HDR network operating at 200/100Gbps. The new Hércules supercomputer also includes six servers, each equipped with an Nvidia A100 GPU with 40GB of memory.

Continuously updated and improved, access to the cluster is managed using the open-source Slurm software, which organizes user workloads into different queues. This queue system ensures that short jobs experience reduced waiting times while preventing long jobs from monopolizing cluster resources, thereby optimizing core usage based on workload distribution across the queues.

Currently, the cluster hosts a wide range of computational software frequently used by research groups, and it remains open to the installation of new software upon user request.

The CICA Operations team is responsible for implementing and maintaining the service, as well as providing technical support to users when needed. This service includes activities such as backup management, system monitoring, and the installation and updating of software.

## Requiered modules

The following modules need to be loaded to compile any code with SOLERA

```console
module load make/4.4.1-GCCcore-12.3.0

module load CMake/3.26.3-GCCcore-12.3.0

module load OpenMPI/4.1.5-GCC-12.3.0

module load GCC/12.3.0

module load Eigen/3.4.0-GCCcore-12.3.0

module load PETSc/3.21.6-foss-2023a
```

You will need to define some enviromental variables:
```console
export SOLERA_DIR=$HOME/DMD
export CGAL_DIR=$HOME/cgal
PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig
````

These commands are embebed in the configure.sh file located in this directory.

## Compilation instructions

* 1º Resources allocation:
```console
salloc --mem=16G -c 4 -t 00:30:00 -p standard srun --pty /bin/bash -i
```
* 2º Compilation:
```console
chmod +x configure-Hercules.sh
./configure-Hercules.sh
exit
```

## Execute a code compiled with SOLERA

* Execution:
```console
sbatch run.sh
```



