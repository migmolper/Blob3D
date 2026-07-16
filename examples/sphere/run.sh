
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

MPI_P=8
${MPI_RUN} -np ${MPI_P} ./exe     

echo "------------------------------------"
echo "    THE SIMULATION HAS FINISHED     "
echo "------------------------------------"
