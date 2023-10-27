#!/bin/bash

#SBATCH --job-name=dvDESI
###SBATCH --output=dvDESIKL-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=qual_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu

module load gsl
module load openmpi3
WORKDIR=/home/u17/jiachuanxu/CosmoLike/KL_WFIRST
cd ${WORKDIR}
### LSST
for (( c=0; c<2; c++ ))
do
	for (( k=0; k<12; k++ ))
	do
		./test_cosmic_shear LSST ${c} ${k}
	done
done

### DESI2
for (( c=0; c<6; c++ ))
do 
	for (( k=0; k<12; k++ ))
	do
		./test_cosmic_shear DESI2 ${c} ${k}
	done
done
