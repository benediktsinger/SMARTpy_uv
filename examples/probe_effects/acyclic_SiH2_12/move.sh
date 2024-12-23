#!/bin/bash

for file in $PWD/*_S_SiH2_12_acyclic_Rh0_cplx.sdf; do
name=$(basename $file _S_SiH2_12_acyclic_Rh0_cplx.sdf)
mv ${name}* DONE/
done

#for file in $PWD/*_S_SiH2_12_cyclic_Rh1_cplx.sdf; do
#echo $file
#name=$(basename $file _S_SiH2_12_cyclic_Rh1_cplx.sdf)
#echo $name
#mv ${name}* DONE/
#done
