#!/bin/bash
#
# Submit equilibration simulations to run on the cluster.
#
#SBATCH --job-name=%residue_number%_%aminoacid_mutant%_equil
#SBATCH --time=1-23:59:59
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ch.sustay@physik.uni-muenchen.de
#

module load amber

cd equilibration

for i in $(ls -d * | sort -n); do
       cd $i
       echo "Equilibrating window $i"

       if [ $i == 0 ]; then
	       inpcrd=minimization.rst7
       else
	       inpcrd=../$((i-1))/equilibration.rst7
       fi

       echo "Heating..."
       pmemd.cuda -i heat.in -c $inpcrd -ref $inpcrd -p topology.parm7 \
	       -O -o heat.out -inf heat.info -r heat.rst7 -x heat.nc \
	       -l heat.log

       echo "Equilibrating..."
       pmemd.cuda -i equilibration.in -c heat.rst7 -p topology.parm7 \
	       -O -o equilibration.out -inf equilibration.inf \
	       -r equilibration.rst7 -x equilibration.nc \
	       -l equilibration.log

       if [ $i < 9 ]; then
	       cp equilibration.rst7 ../../RE/equil.rep.00$((i+1))
       else
	       cp equilibration.rst7 ../../RE/equil.rep.0$((i+1))
       fi
       
       cd ..
done

cd ..
echo "Done."