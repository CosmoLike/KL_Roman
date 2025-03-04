#!/bin/bash

#SBATCH --job-name=covSKA
#SBATCH --output=log/covSKA-%A_%a.out
#SBATCH --error=log/covSKA-%A_%a.err
#SBATCH --array=1-770
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=1:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

if [ "$#" -ne 1 ]; then
    echo "Error: Need to pass the parameter file name" >&2
    exit 1
fi

module load gsl
WORKDIR=/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST
cd ${WORKDIR}

for (( c=0; c<4; c++ ))
do
	hit=$(( ${SLURM_ARRAY_TASK_ID} + c * 770 ))
	./compute_covariances_fourier ${hit} $1
done



