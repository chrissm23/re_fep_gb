#!/bin/bash
#
# Submit REMD to cluster
#
#SBATCH --job-name=%residue_number%_%aminoacid_mutant%_remd
#SBATCH --time=4-23:59:59
#SBATCH --nodes=%n_nodes%
#SBATCH --ntasks=%n_tasks%
#SBATCH --gres=gpu:%n_gpus%
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ch.sustay@physik.uni-muenchen.de
#

module load amber

cd RE

echo "Running H-REMD..."
mpirun -n %n_cores% pmemd.cuda.MPI -ng %n_replicas% -groupfile groupfile

cd ..
echo "Done."
