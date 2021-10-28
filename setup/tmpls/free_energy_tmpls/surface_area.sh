#!/bin/bash
#
# Run Solvent Accessible Surface Area calculation.
#

cd WT

echo "Calculating SASA of WT..."
pmemd.cuda -i sasa.in -c equilibration.rst7 -p topology.parm7 \
	-O -o sasa.out -inf sasa.info -r sasa.rst7 -x sasa.nc -l sasa.log

cd ../MT

echo "Calculating SASA of MT..."
pmemd.cuda -i sasa.in -c equilibration.rst7 -p topology.parm7 \
       -O -o sasa.out -inf sasa.inf -r sasa.rst7 -x sasa.nc -l sasa.log

cd ..
echo "Done."
