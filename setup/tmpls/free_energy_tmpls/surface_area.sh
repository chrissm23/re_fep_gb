#!/bin/bash
#
# Run Solvent Accessible Surface Area calculation.
#

cd WT

cp ../../equilibration/WT/0/equilibration.rst7 equilibration.rst7
echo "Calculating SASA of WT..."
pmemd.cuda -i sasa.in -c equilibration.rst7 -p topology.parm7 \
	-O -o sasa.out -inf sasa.info -r sasa.rst7 -x sasa.nc -l sasa.log

cd ../MT

cp ../../equilibration/MT/0/equilibration.rst7 equilibration.rst7
echo "Calculating SASA of MT..."
pmemd.cuda -i sasa.in -c equilibration.rst7 -p topology.parm7 \
       -O -o sasa.out -inf sasa.inf -r sasa.rst7 -x sasa.nc -l sasa.log

cd ..
echo "Done."
