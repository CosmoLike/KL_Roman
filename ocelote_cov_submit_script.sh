#!/bin/bash

#SBATCH --job-name=covPIT
#SBATCH --output=/xdisk/timeifler/jiachuanxu/job_logs/covRomanPIT-%A_%a.out
#SBATCH --array=1-999
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=0:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu

module load gsl
module load openmpi3
WORKDIR=/home/u17/jiachuanxu/CosmoLike/KL_WFIRST

### 5 neff x 5 area x  4 ell_max scenarios
### 10 tomography bins, 55 power spectra
### In total 30,800 covariance blocks
### Submit in 880 arrays with 35 runs per job

cd ${WORKDIR}
for (( c=0; c<155; c++ ))
do
	hit=$(( ${SLURM_ARRAY_TASK_ID} + c * 999 ))
	./compute_covariances_fourier ${hit}
done

#154000




#for (( c=0; c<30; c++ ))
#do
#	bad_id=$(( ${SLURM_ARRAY_TASK_ID} + c * 2 - 1 ))
#    hit=${bad_hits[$bad_id]}
#	./compute_covariances_fourier ${hit}
#done

