#!/bin/bash
#
# Run minimization of WT for later equilibraion of all topologies
#

pmemd=$AMBERHOME/bin/pmemd.MPI
fe_dir=./FE/minimization/%wt_or_mt%

echo "Minimizing..."

mpirun -n 6 $pmemd -i $fe_dir/minimization.in -p $fe_dir/topology.parm7 \
	-c $fe_dir/coordinates.rst7 -ref $fe_dir/coordinates.rst7 \
	-O -o $fe_dir/minimization.out -inf $fe_dir/min.info \
	-r $fe_dir/minimization.rst7

echo "Minimization finished."
