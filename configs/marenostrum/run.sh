#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=test-SOLERA

# Directorio de trabajo
#SBATCH --chdir=.

# Archivo de salida
#SBATCH --output=batch-%x-%j.out                  

# Archivo de errores
#SBATCH --error=batch-%x-%j.err                   

# Partición a utilizar
#SBATCH --partition=standard                

# Tiempo máximo de ejecución (horas:minutos:segundos)
#SBATCH --time=24:00:00                    

# Número total de tareas MPI
#SBATCH --ntasks=8             

# Número máximo de CPUs por task
#SBATCH --cpus-per-task=16

## Common modules
module purge
module load gcc/13.2.0
module load cmake/3.29.2
module load openmpi/4.1.6-gcc
module load ucx/1.16.0-gcc
module load lapack/3.12.1-gcc
module load petsc/3.24.1-gcc12.3-ompi4.1.6-double

# Establecer variables de entorno específicas si es necesario.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
MPI_RUN=mpirun
MPI_P=$SLURM_NTASKS
SOLERA_EXE=exe-tasting-SOLERA

# Binario en build/test/ respecto a la raíz del repo. Convención: lanzar sbatch
# desde la raíz de DMD (p. ej. sbatch configs/marenostrum/run.sh); --chdir=. deja
# el cwd en SLURM_SUBMIT_DIR. Opcional: SOLERA_DIR apunta al clon si no coincide.
if [[ -n "${SOLERA_DIR}" ]] && [[ -d "${SOLERA_DIR}/build/test" ]]; then
  cd "${SOLERA_DIR}" || exit 1
elif [[ -n "${SLURM_SUBMIT_DIR}" ]]; then
  cd "${SLURM_SUBMIT_DIR}" || exit 1
else
  _here="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
  cd "$(cd "${_here}/../.." && pwd)" || exit 1
fi

EXE_DIR="build/test"
REL_EXE="${EXE_DIR}/exe-tasting-SOLERA"

if [[ ! -x "${REL_EXE}" ]]; then
  echo "No se encuentra el binario de tests en ${REL_EXE} (cwd: $(pwd))."
  echo "Compila antes con ./configure.sh (modo Lib o Lib&Test)."
  exit 1
fi

# Ruta absoluta: después hacemos "cd test" y una ruta relativa rompería mpirun.
EXE_BIN="$(pwd)/${REL_EXE}"

# Ejecutamos desde el directorio de tests para que las rutas relativas
# de ficheros de entrada (por ejemplo IO/inputs/...) sean válidas.
cd "test" || exit 1

# Si pasas argumentos a este script, se los reenviamos como expresión de filtro de Catch
if [[ $# -gt 0 ]]; then
  ${MPI_RUN} -np ${MPI_P} "${EXE_BIN}" "$@"
else
  # Ejemplo por defecto (puedes editarlo a tu gusto)
  ${MPI_RUN} -np ${MPI_P} "${EXE_BIN}" "[test-Mass-Transport-Master-Equation-PETSc][MgHx-hcp-cube-x12-x6-x6-Forward-Euler]"
fi
