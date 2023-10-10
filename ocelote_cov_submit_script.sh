#!/bin/bash

#SBATCH --job-name=covDESI
###SBATCH --output=covDESIKL-%A_%a.out
#SBATCH --array=1-770
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=0:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu

module load gsl
module load openmpi3
WORKDIR=/home/u17/jiachuanxu/CosmoLike/KL_WFIRST
cd ${WORKDIR}
for (( c=0; c<4; c++ ))
do
	hit=$(( ${SLURM_ARRAY_TASK_ID} + c * 770 ))
	./compute_covariances_fourier ${hit}
done



