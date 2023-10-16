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
for (( c=0; c<2; c++ ))
do
	for bary in dmo mb2 illustris eagle HzAGN TNG100 cowls_AGN cowls_AGN_T8p5 cowls_AGN_T8p7 BAHAMAS BAHAMAS_T7p6 BAHAMAS_T8p0
	do
		./like_fourier ${c} 0 shear_shear ${bary}
	done
done
