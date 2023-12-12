#!/bin/bash

#SBATCH --job-name=likeDSA
#SBATCH --output=log/likeDSA-%A_%a.out
#SBATCH --error=log/likeDSA-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

module load gsl
module swap openmpi3 mpich/3.3.1
WORKDIR=/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST
cd ${WORKDIR}

./like_fourier 0 0 shear_shear dmo
