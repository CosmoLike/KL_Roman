#!/bin/bash

#SBATCH --job-name=covDESIKL
#SBATCH --output=/xdisk/timeifler/jiachuanxu/job_logs/covDESI2KL-%A.out
#SBATCH --array=1-220
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

###SBATCH --partition=high_priority
###SBATCH --qos=user_qos_timeifler
#SBATCH --partition=standard
#SBATCH --qos=qual_qos_timeifler

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

for (( c=0; c<1; c++ ))
do
	hit=$(( ${SLURM_ARRAY_TASK_ID} + c * 220))
	./compute_covariances_fourier ${hit}
done

#154000

#bad_hits=(101830 38984 118343 93084 110429 114959 43773 89056 121472 138312 78672 96869 80994 52477 54776 75436 24626 9226 116719 77877 76207 53868 130695 50416 147767 17127 80722 44557 136281 85873 112924 52423 138053 6760 130830 66587 80082 42746 109405 153900 21925 124588 87316 19900 48464 85005 69245 128897 22782 116690 31140 129592 77005 11778 112485 116448 106961 27541 21461 34562 30061 150596 120287 92082 13542 18770 54900 49519 57891 81581 116282 137128 47556 113548 74940 84598)

#for (( c=0; c<4; c++ ))
#do
#	bad_id=$(( ${SLURM_ARRAY_TASK_ID} + c * 19 - 1 ))
#    hit=${bad_hits[$bad_id]}
#	./compute_covariances_fourier ${hit}
#done

