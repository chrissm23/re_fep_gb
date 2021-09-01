#!/bin/bash
#
# Run genremdinputs.py to generate all necessary files for REMD
#

genremdinputs.py -inputs hamiltonians.dat -groupfile groupfile.ref \
	-i mdin.ref -O
