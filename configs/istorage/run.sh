#$ -S /bin/bash
#$ -N tasting-SOLERA
#$ -wd /home/mmolinos/DMD
#$ -m a
#$ -M mmolinos@us.es
#$ -o tasting-SOLERA.out
#$ -e tasting-SOLERA.err
#$ -q all.q@istorage-03
#$ -pe mpi 8

#PETSC_VERSION="Debug-3.21"
PETSC_VERSION="Release-3.21"


if [[ "$PETSC_VERSION" == "Release-3.21" ]]; then
    echo Release version
    module load gcc-10.2.0
    module load cmake-3.24.0
    module load petsc-3.21.0-openmpi-slepc-nodebug
    MPI_RUN=/home/software/petsc-3.21.0/installation/bin/mpiexec

elif [[ "$PETSC_VERSION" == "Debug-3.21" ]]; then
    echo Debug version
    module load gcc-10.2.0
    module load cmake-3.24.0
    module load petsc-3.21.0-openmpi-slepc-debug
    MPI_RUN=/home/software/petsc-3.21.0/installation-debug/bin/mpiexec
else
    echo "Unrecognised option" $PETSC_VERSION
    exit
fi


MPI_P=8
export OMP_NUM_THREADS=10

##./exe-tasting-SOLERA "[test-ADP-MgHx][init-adp]"

##./exe-tasting-SOLERA "[boundary-conditions][apply-symmetry-kinematic-restrictions]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[boundary-conditions][enforce_periodic_boundary_conditions]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[init-DMD-simulation-DMDA]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[init-DMD-simulation][domain-decomposition-Mg-cube]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-init-DMD-simulation][domain-decomposition-MgHx-cube]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[init-DMD-simulation][domain-decomposition-hex-nanowire-Mg-R25A]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[Numerical][Search-neighbors]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[Numerical][Search-diffusive-sites]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-compute-energy-density][cube-Mg-x17-x10-x10-periodic]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][cube-Mg-x20-x15-x15-periodic]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][cube-Mg-x17-x10-x10-finite]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][cube-Mg-x5-x5-x5]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][Mg-hcp-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][nanowire-MgHx-60A]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][hex-nanowire-Mg-R-25A]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][nanocube-17-x-10-x-10]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][nanoparticle-Mg-20A]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-initialize-stdv-q][nanoparticle-MgHx-20A]"


##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation][Mg-hcp-box-x20-x15-x15-finite]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation][Mg-R-20A-nanoparticle]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation][Mg-hcp-nanowire]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation-bulk][Mg-hcp-bulk-x20-x15-x15-periodic]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation-bulk][MgHx-rutile-bulk-x5-x5-x5]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation-bulk][Mg-hcp-nanowire]"
${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation-bulk][Al-fcc-cube-x10-x10-x10-periodic]"
${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation-bulk][Cu-fcc-cube-x10-x10-x10-periodic]"
${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation-bulk][Al2Cu-fcc-x5-x5-x5-periodic]"

${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-dF-SNES-PETSc][Mg-hcp-bulk-x20-x15-x15]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-dF-SNES-PETSc][Mg-hcp-bulk-x12-x9-x9]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-dF-SNES-PETSc][Mg-hcp-bulk-x12-x9-x9]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-dF-SNES-PETSc][Mg-hcp-nanowire]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-dF-TAO-PETSc][Mg-hcp-bulk-x20-x15-x15]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-stdv-q-SNES-PETSc][Mg-hcp-periodic-cell-x12-x9-x9]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-stdv-q-TAO-PETSc][Mg-hcp-periodic-cell-x9-x4-x4]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA  "[min-F0-d-stdv-q-SNES-PETSc][MgHx-hcp-periodic-cell]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA  "[min-F0-d-stdv-q-TAO-PETSc][MgHx-hcp-periodic-cell]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA  "[min-F0-d-stdv-q-SNES-PETSc][Mg-hcp-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA  "[min-F0-d-stdv-q-SNES-PETSc][MgH2-rutile-periodic-cell]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA  "[min-F0-d-stdv-q-TAO-PETSc][Mg-hcp-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA  "[min-F0-d-stdv-q-TAO-PETSc][MgH2-rutile-periodic-cell]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA  "[min-F0-d-mean-q-SNES-PETSc][Mg-hcp-circ-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-SNES-PETSc][Mg-hcp-hex-nanowire-R-25]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-monolithic-SNES-PETSc][Mg-hcp-cube-x9-x4-x4-free]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-monolithic-SNES-PETSc][Mg-system-nanoparticle]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-SNES-PETSc][MgHx-system-nanoparticle]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-SNES-PETSc][MgH2-system-nanoparticle]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-SNES-PETSc][Mg-hcp-circ-section-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-SNES-PETSc][Mg-hcp-hex-section-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-SNES-PETSc][MgHx-hcp-circ-section-nanowire]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-d-mean-q-stdv-q-SNES-PETSc][Mg-hcp-hex-nanowire-R-25]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-min-F0-dF-d-stdv-q-SNES-PETSc][MgH2-rutile-periodic-cell-x3-x3-x3]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-dF-d-stdv-q-SNES-PETSc][Mg-hcp-periodic-cell-x12-x9-x9]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[min-F0-dF-d-stdv-q-TAO-PETSc][Mg-hcp-periodic-cell-x9-x4-x4]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-Thermostat-TAO-PETSc][MgH2-rutile-periodic-cell-x3-x3-x3]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-IO-dump][Mg-MgHx-system]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-DMSwarmWriteHDF5][Mg-MgHx-system]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-DMSwarmReadHDF5][Mg-MgHx-system]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[enforce-Neumann-boundary-condition][MgHx-hcp-x4-x3-x5-XY-periodic]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[enforce-Neumann-boundary-condition][MgHx-hcp-x6-x4-x2-Z-periodic]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[enforce-Neumann-boundary-condition][MgHx-hcp-nanowire-Z-periodic]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-chemical-multiplier-MgHx][chemical-potential-sphere-bcc]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-d-mf-gamma][MgHx-hcp-cube-x9-x6-x6]"


##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mf-critical-state-energy-barrier][MgHx-hcp-x7-x5-x3]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-Compute-Arrhenius-TM][MgHx-hcp-cube-x9-x6-x6-1H-middle]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-Mass-Transport-Master-Equation-PETSc-Forward-Euler][MgHx-hcp-cube-x12-x6-x6]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-Mass-Transport-Master-Equation-PETSc-Backward-Euler][MgHx-hcp-cube-x12-x6-x6]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-Mass-Transport-Master-Equation-PETSc-Runge-Kutta][MgHx-hcp-cube-x12-x6-x6]"

##./exe-tasting-SOLERA "[Numerical][Measure-evaluation-Vij]"
##./exe-tasting-SOLERA "[Numerical][Measure-evaluation-Vijk]"
##./exe-tasting-SOLERA "[Numerical][Integral-evaluation-6d-3th]"
##./exe-tasting-SOLERA "[Numerical][Integral-evaluation-9d-3th]"

##./exe-tasting-SOLERA "[Numerical][Integral-evaluation-6d-mp]"
##./exe-tasting-SOLERA "[Numerical][Integral-evaluation-9d-mp]"

##./exe-tasting-SOLERA "[mf-MgHx-adp][energy-contributions-dV0-dstdvq]"
##./exe-tasting-SOLERA "[mf-MgHx-adp][energy-contributions]"

##./exe-tasting-SOLERA "[mf-MgHx-adp][gradient-energy-contributions]"
##./exe-tasting-SOLERA "[MgHx-adp][energy-contributions-MgHx-adp]"
##./exe-tasting-SOLERA "[test-ADP-MgHx][init-adp]"
##./exe-tasting-SOLERA "[test-ADP-MgHx][rho]"
##./exe-tasting-SOLERA "[test-ADP-MgHx][V-phi]"
##./exe-tasting-SOLERA "[test-ADP-MgHx][V-mu]"
##./exe-tasting-SOLERA "[test-ADP-MgHx][V-lambda]"
##./exe-tasting-SOLERA "[test-mf-MgHx-adp][V-phi-vs-mf-V-phi]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[MgHx-mf-V-bulk][evaluate-V0]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[MgHx-mf-V-bulk][evaluate-DV0-Dxi]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[MgHx-mf-V-bulk][evaluate-DL0-DU]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[MgHx-mf-V-bulk][evaluate-DL0-Dstdvq]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[MgHx-mf-V-bulk][comparison-MP-vs-GH3th]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[mechanical-relaxation][Mg-system]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-Thermostat-bulk-SNES-PETSc][Mechanical-eqs]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-Thermal-Expansion][Mechanical-eqs]"
##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA  "[test-Thermal-Expansion][aniso-bulk-TAO-PETSc]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[boundary-conditions][periodic-box-expansion]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[Numerical][Search-neighbors-periodic]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[test-Equilibrium-U-q-TAO-PETSc][Mechanical-eqs]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[compute-1-Hx-diffusion][Chemical-MgHx]"

##${MPI_RUN} -np ${MPI_P} ./exe-tasting-SOLERA "[compute-rate-molar-fraction-single-Hx][Chemical-eqs]"
