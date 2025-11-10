#!/bin/bash
#SBATCH --job-name=RomanKL
#SBATCH --output=RomanKL-%A_%a.out

### 1. puma
###SBATCH --nodes=1
###SBATCH --ntasks-per-node=80
###SBATCH --ntasks-per-socket=40
###SBATCH --cpus-per-task=1
###SBATCH --partition=high_priority
###SBATCH --qos=user_qos_timeifler
###SBATCH --account=timeifler
###SBATCH --time=96:00:00

### 2. ocelote
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --array=1-770
#SBATCH --qos=qual_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=96:00:00

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu


module load gsl
module swap openmpi3/3.1.4 mpich/3.3.1

cd /home/u17/jiachuanxu/CosmoLike/KL_WFIRST
for i in $(seq 1 2);
do
IDX=$(( ${SLURM_ARRAY_TASK_ID} + ( i - 1 ) * 770 ))
./compute_covariances_fourier $IDX
done


