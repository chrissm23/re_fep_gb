#!/bin/bash
#
# Create 10 restart files for REMD from long configuration sample
#

cpptraj -p topology.parm7 -y sample.nc <<_EOF

for i=500;i<5001;i+=500
  trajout restart_\$i.rst7 restartnc onlyframes \$i
done

run

quit
_EOF

pattern="restart_*.rst7"
restart_files=( $pattern )
wt_or_mt=${PWD##*/}

for i in ${!restart_files[@]}; do
  if [ ! -d "../../RE/$wt_or_mt/$i" ]; then
    mkdir ../../RE/$wt_or_mt/$i
  fi
  cp ${restart_files[$i]} ../../RE/$wt_or_mt/$i/restart.rst7
done