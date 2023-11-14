#!/bin/bash

#SBATCH --job-name=covDSA
#SBATCH --output=log/covDSA-%A_%a.out
#SBATCH --array=1-770
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --qos=qual_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=0:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

module load gsl
module swap openmpi3 mpich/3.3.1
WORKDIR=/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST
cd ${WORKDIR}
for (( c=0; c<4; c++ ))
do
	hit=$(( ${SLURM_ARRAY_TASK_ID} + c * 770 ))
	./compute_covariances_fourier ${hit}
done



