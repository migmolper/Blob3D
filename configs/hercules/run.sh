#!/bin/bash
#SBATCH --job-name=tasting-SOLERA               # Nombre del trabajo
#SBATCH --output=%x-%j.out                      # Archivo de salida
#SBATCH --error=%x-%j.err                       # Archivo de errores
#SBATCH --partition=standard                    # Partición a utilizar
#SBATCH --time=24:00:00                         # Tiempo máximo de ejecución (horas:minutos:segundos)
#SBATCH --ntasks=8                              # Número total de tareas MPI
#SBATCH --ntasks-per-node=3                     # Número de tareas en cada nodo
#SBATCH --cpus-per-task=16                      # N?mero máximo de CPUs por tarea
#SBATCH --mem=16G                               # Memoria total para el trabajo

## Configure environments
clear
module purge
module load GCC/12.3.0
module load CMake/3.26.3-GCCcore-12.3.0
module load PETSc/3.21.6-foss-2023a

# Establecer variables de entorno específicas si es necesario.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
MPI_RUN=mpirun
MPI_P=$SLURM_NTASKS
SOLERA_EXE=exe-tasting-SOLERA

##./${SOLERA_EXE} "[test-ADP-MgHx][init-adp]"

##./${SOLERA_EXE} "[boundary-conditions][apply-symmetry-kinematic-restrictions]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[boundary-conditions][enforce_periodic_boundary_conditions]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[init-DMD-simulation-DMDA]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[init-DMD-simulation][domain-decomposition-Mg-cube]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-init-DMD-simulation][domain-decomposition-MgHx-cube]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[init-DMD-simulation][domain-decomposition-hex-nanowire-Mg-R25A]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[Numerical][Search-neighbors]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[Numerical][Search-diffusive-sites]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-compute-energy-density][cube-Mg-x17-x10-x10-periodic]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][cube-Mg-x20-x15-x15-periodic]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][cube-Mg-x17-x10-x10-finite]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][cube-Mg-x5-x5-x5]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][Mg-hcp-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][nanowire-MgHx-60A]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][hex-nanowire-Mg-R-25A]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][nanocube-17-x-10-x-10]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][nanoparticle-Mg-20A]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-initialize-stdv-q][nanoparticle-MgHx-20A]"


##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[mechanical-relaxation][Mg-hcp-box-x20-x15-x15-finite]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[mechanical-relaxation][Mg-R-20A-nanoparticle]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[mechanical-relaxation][Mg-hcp-nanowire]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[mechanical-relaxation-bulk][Mg-hcp-bulk-x20-x15-x15-periodic]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[mechanical-relaxation-bulk][MgHx-rutile-bulk-x5-x5-x5]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[mechanical-relaxation-bulk][Mg-hcp-nanowire]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-dF-SNES-PETSc][Mg-hcp-bulk-x20-x15-x15]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-dF-SNES-PETSc][Mg-hcp-bulk-x12-x9-x9]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-dF-SNES-PETSc][Mg-hcp-bulk-x12-x9-x9]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-dF-SNES-PETSc][Mg-hcp-nanowire]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-dF-TAO-PETSc][Mg-hcp-bulk-x20-x15-x15]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-stdv-q-SNES-PETSc][Mg-hcp-periodic-cell-x12-x9-x9]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-stdv-q-TAO-PETSc][Mg-hcp-periodic-cell-x9-x4-x4]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-stdv-q-TAO-PETSc][Mg-hcp-periodic-cell-x12-x9-x9]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE}  "[min-F0-d-stdv-q-SNES-PETSc][MgHx-hcp-periodic-cell]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE}  "[min-F0-d-stdv-q-TAO-PETSc][MgHx-hcp-periodic-cell]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE}  "[min-F0-d-stdv-q-SNES-PETSc][Mg-hcp-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE}  "[min-F0-d-stdv-q-SNES-PETSc][MgH2-rutile-periodic-cell]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE}  "[min-F0-d-stdv-q-TAO-PETSc][Mg-hcp-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE}  "[min-F0-d-stdv-q-TAO-PETSc][MgH2-rutile-periodic-cell]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE}  "[min-F0-d-mean-q-SNES-PETSc][Mg-hcp-circ-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-mean-q-SNES-PETSc][Mg-hcp-hex-nanowire-R-25]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-mean-q-stdv-q-monolithic-SNES-PETSc][Mg-hcp-cube-x9-x4-x4-free]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-mean-q-stdv-q-monolithic-SNES-PETSc][Mg-system-nanoparticle]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-mean-q-stdv-q-SNES-PETSc][MgHx-system-nanoparticle]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-mean-q-stdv-q-SNES-PETSc][MgH2-system-nanoparticle]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-mean-q-stdv-q-SNES-PETSc][Mg-hcp-circ-section-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-mean-q-stdv-q-SNES-PETSc][Mg-hcp-hex-section-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-d-mean-q-stdv-q-SNES-PETSc][MgHx-hcp-circ-section-nanowire]"
${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-SNES-PETSc][Mg-hcp-nanowire]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-min-F0-dF-d-stdv-q-SNES-PETSc][MgH2-rutile-periodic-cell-x3-x3-x3]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-dF-d-stdv-q-SNES-PETSc][Mg-hcp-periodic-cell-x12-x9-x9]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[min-F0-dF-d-stdv-q-TAO-PETSc][Mg-hcp-periodic-cell-x9-x4-x4]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-Thermostat-TAO-PETSc][MgH2-rutile-periodic-cell-x3-x3-x3]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-IO-dump][Mg-MgHx-system]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-DMSwarmWriteHDF5][Mg-MgHx-system]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-DMSwarmReadHDF5][Mg-MgHx-system]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[enforce-Neumann-boundary-condition][MgHx-hcp-x4-x3-x5-XY-periodic]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[enforce-Neumann-boundary-condition][MgHx-hcp-x6-x4-x2-Z-periodic]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[enforce-Neumann-boundary-condition][MgHx-hcp-nanowire-Z-periodic]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-chemical-multiplier-MgHx][chemical-potential-sphere-bcc]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-d-mf-gamma][MgHx-hcp-cube-x9-x6-x6]"


##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[mf-critical-state-energy-barrier][MgHx-hcp-x7-x5-x3]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-Compute-Arrhenius-TM][MgHx-hcp-cube-x12-x6-x6]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-Mass-Transport-Master-Equation-PETSc-Forward-Euler][MgHx-hcp-cube-x12-x6-x6]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-Mass-Transport-Master-Equation-PETSc-Backward-Euler][MgHx-hcp-cube-x12-x6-x6]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-Mass-Transport-Master-Equation-PETSc-Runge-Kutta][MgHx-hcp-cube-x12-x6-x6]"

##./${SOLERA_EXE} "[Numerical][Measure-evaluation-Vij]"
##./${SOLERA_EXE} "[Numerical][Measure-evaluation-Vijk]"
##./${SOLERA_EXE} "[Numerical][Integral-evaluation-6d-3th]"
##./${SOLERA_EXE} "[Numerical][Integral-evaluation-9d-3th]"

##./${SOLERA_EXE} "[Numerical][Integral-evaluation-6d-mp]"
##./${SOLERA_EXE} "[Numerical][Integral-evaluation-9d-mp]"

##./${SOLERA_EXE} "[mf-MgHx-adp][energy-contributions-dV0-dstdvq]"
##./${SOLERA_EXE} "[mf-MgHx-adp][energy-contributions]"

##./${SOLERA_EXE} "[mf-MgHx-adp][gradient-energy-contributions]"
##./${SOLERA_EXE} "[MgHx-adp][energy-contributions-MgHx-adp]"
##./${SOLERA_EXE} "[test-ADP-MgHx][init-adp]"
##./${SOLERA_EXE} "[test-ADP-MgHx][rho]"
##./${SOLERA_EXE} "[test-ADP-MgHx][V-phi]"
##./${SOLERA_EXE} "[test-ADP-MgHx][V-mu]"
##./${SOLERA_EXE} "[test-ADP-MgHx][V-lambda]"
##./${SOLERA_EXE} "[test-mf-MgHx-adp][V-phi-vs-mf-V-phi]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[MgHx-mf-V-bulk][evaluate-V0]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[MgHx-mf-V-bulk][evaluate-DV0-Dxi]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[MgHx-mf-V-bulk][evaluate-DL0-DU]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[MgHx-mf-V-bulk][evaluate-DL0-Dstdvq]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[MgHx-mf-V-bulk][comparison-MP-vs-GH3th]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[mechanical-relaxation][Mg-system]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-Thermostat-bulk-SNES-PETSc][Mechanical-eqs]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-Thermal-Expansion][Mechanical-eqs]"
##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE}  "[test-Thermal-Expansion][aniso-bulk-TAO-PETSc]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[boundary-conditions][periodic-box-expansion]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[Numerical][Search-neighbors-periodic]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[test-Equilibrium-U-q-TAO-PETSc][Mechanical-eqs]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[compute-1-Hx-diffusion][Chemical-MgHx]"

##${MPI_RUN} -np ${MPI_P} ./${SOLERA_EXE} "[compute-rate-molar-fraction-single-Hx][Chemical-eqs]"
