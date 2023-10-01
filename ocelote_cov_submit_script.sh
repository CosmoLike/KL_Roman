#!/bin/bash

#SBATCH --job-name=covDESI
###SBATCH --output=covDESIKL-%A_%a.out
#SBATCH --array=1-990
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=240:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu

module load gsl
module load openmpi3
WORKDIR=/home/u17/jiachuanxu/CosmoLike/KL_WFIRST
cd ${WORKDIR}
for (( c=0; c<2; c++ ))
do
	hit=$(( ${SLURM_ARRAY_TASK_ID} + c * 990 ))
	./compute_covariances_fourier ${hit}
done



