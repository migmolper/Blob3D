# How to compile SOLERA in MareNostrum 5 (BSC)

- Date: 19/01/2026
- Author: M.Molinos

The Barcelona Supercomputing Center – Centro Nacional de Supercomputación (BSC-CNS) provides the national and international research community with access to MareNostrum, one of Europe’s leading high-performance computing (HPC) infrastructures. MareNostrum is a large-scale supercomputing system composed of thousands of compute nodes, delivering several petaflops of peak performance, supported by a high-capacity memory subsystem and a multi-petabyte parallel storage system designed to efficiently handle data-intensive workloads. Node interconnectivity is ensured by a high-bandwidth, low-latency InfiniBand network, enabling efficient large-scale parallel computations.

The system is continuously upgraded and maintained, and user access is managed through the Slurm workload manager, which organizes computational jobs into multiple partitions and queues. This scheduling strategy ensures fair and efficient resource allocation, minimizing queue times for short and medium jobs while preventing long-running simulations from monopolizing computational resources, thereby maximizing overall system throughput.
MareNostrum provides users with access to a comprehensive software stack, including optimized compilers, mathematical libraries, and widely used scientific applications for large-scale simulations and data analysis. Additional software can be installed upon request, ensuring flexibility to accommodate the specific needs of different research projects.

The operation, maintenance, and user support of MareNostrum are carried out by the BSC Operations and User Support teams, who ensure reliable system performance and provide technical assistance. These services include system monitoring, software installation and updates, data management support, and user guidance to facilitate efficient and productive use of the infrastructure.

## Requiered modules

The following modules need to be loaded to compile any code with SOLERA

```console
module load oneapi/2023.2.0

module load gcc/13.2.0

module load cmake/3.29.2

module load eigen/3.4.0

module load petsc/3.23.0-mumps
```

You will need to define some enviromental variables:
```console
export SOLERA_DIR=$HOME/DMD
export CGAL_DIR=$HOME/cgal
PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$PETSC_DIR/$PETSC_ARCH/lib/pkgconfig
````

These commands are embebed in the configure.sh file located in this directory.

## Compilation instructions

* Compilation:
```console
chmod +x configure-Marenostrum.sh
./configure-Marenostrum.sh
exit
```

## Execute a code compiled with SOLERA

* Execution (from the **DMD repository root**):
```console
sbatch configs/marenostrum/run.sh
```



