#!/bin/bash
#
# Run minimization of WT for later equilibraion of all topologies
#

pmemd=$AMBERHOME/bin/pmemd

echo "Minimizing"

$pmemd -i minimization.in -p topology.parm7 -c coordinates.rst7 \
	-ref coordinates.rst7 -O -o minimization.out -inf min.info \
	-r minimization.rst7

echo "Done"
