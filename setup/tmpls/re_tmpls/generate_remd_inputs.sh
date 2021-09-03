#!/bin/bash
#
# Run genremdinputs.py to generate all necessary files for REMD
#

fe_dir=FE/RE/%wt_or_mt

cd $fe_dir

genremdinputs.py -inputs hamiltonians.dat \
	-groupfile groupfile.ref -i mdin.ref -O

cd ../../
