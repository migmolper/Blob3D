
#!/bin/bash

PETSC_VERSION="Release-3.22.1"
##PETSC_VERSION="Debug-3.22.1"

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    LOCAL_HOST_NAME=$HOSTNAME
elif [[ "$OSTYPE" == "darwin"* ]]; then
    LOCAL_HOST_NAME="$(scutil --get LocalHostName)"
else
    echo "Failure: Unsuported OS" >&2
    exit 1
fi

echo We are in $LOCAL_HOST_NAME

## Configure environments
if [[ "$PETSC_VERSION" == "Release-3.22.1" ]]; then
MPI_RUN=~/petsc/arch-darwin-c-release/bin/mpirun
else
MPI_RUN=~/petsc/arch-darwin-c-debug/bin/mpirun
fi

export OMP_NUM_THREADS=1

MPI_P=8
SIZE_MPI_X=2
SIZE_MPI_Y=2
SIZE_MPI_Z=2

NUMBER_STEPS=60
KAPPA=1
${MPI_RUN} -np ${MPI_P} ./exe \
${SIZE_MPI_X} ${SIZE_MPI_Y} ${SIZE_MPI_Z} ${NUMBER_STEPS} ${KAPPA} \
-minJKO_dx_tao_gatol 1.e-5 \
-minJKO_dx_tao_gttol 1.e-6 \
-minJKO_dx_tao_max_it 100 \
-minJKO_dx_tao_type cg \
-minJKO_dx_tao_cg_type prp \
-minJKO_dx_tao_cg_eta 0.01 \
-minJKO_dx_tao_max_funcs 100 \
-minJKO_dx_tao_monitor_globalization \
-minJKO_dx_tao_converged_reason \
-log_view

echo "------------------------------------"
echo "    THE SIMULATION HAS FINISHED     "
echo "------------------------------------"
