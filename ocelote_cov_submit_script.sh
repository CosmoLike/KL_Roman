#!/bin/bash
#SBATCH --job-name=cov3x2
#SBATCH --output=log/cov3x2-%A_%a.out

### puma
##SBATCH --array=1-770
##SBATCH --ntasks=1
##SBATCH --cpus-per-task=1
##SBATCH --partition=high_priority
##SBATCH --qos=user_qos_timeifler
##SBATCH --account=timeifler
##SBATCH --time=2:00:00

### ocelote
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --array=1-770
#SBATCH --qos=qual_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=96:00:00

#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

module load gsl
module load openmpi5
WORKDIR=/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST

cd ${WORKDIR}

for (( c=0; c<3; c++ ))
do
    hit=$(( ${SLURM_ARRAY_TASK_ID} + c * 770 ))
    ./compute_covariances_fourier ${hit} 1
done
