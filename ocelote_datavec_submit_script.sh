#!/bin/bash

#SBATCH --job-name=dv3x2
#SBATCH --output=log/dv3x2-%A.out
#SBATCH --error=log/dv3x2-%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

if [ "$#" -ne ]; then
    echo "Error: Need to pass the parameter file" >&2
    exit 1
fi

module load gsl
module load openmpi5
WORKDIR=/home/u15/yhhuang/cosmology/CosmoLike/KL_WFIRST

cd $WORKDIR

# Use: ./like_fourier [like_flag]
# if like_flag is 0, it computes the data vector and saves it to a file
# if like_flag is 1, it computes the likelihood and returns it
./like_fourieri $1 0
