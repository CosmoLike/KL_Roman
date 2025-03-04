#!/bin/bash

#SBATCH --job-name=likeSKA
#SBATCH --output=log/likeSKA-%A.out
#SBATCH --error=log/likeSKA-%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

if [ "$#" -ne 1 ]; then
    echo "Error: Need to pass the parameter file name" >&2
    exit 1
fi

module load gsl
WORKDIR=/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST
cd ${WORKDIR}

./like_fourier shear_shear dmo $1
