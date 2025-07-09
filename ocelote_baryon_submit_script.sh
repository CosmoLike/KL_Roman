#!/bin/bash

#SBATCH --job-name=bary
#SBATCH --output=log/bary-%A.out
#SBATCH --error=log/bary-%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --account=timeifler
#SBATCH --time=01:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

if [ "$#" -ne 1 ]; then
    echo "Error: Need to pass the parameter file" >&2
    exit 1
fi

module load gsl
module load openmpi5

echo Running on host `hostname`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOB_ID
cd $SLURM_SUBMIT_DIR

sims=("dmo" "mb2" "illustris" "eagle" "HzAGN" "TNG100" 
    "owls_AGN" "owls_DBLIMFV1618" "owls_NOSN" "owls_NOSN_NOZCOOL" "owls_NOZCOOL"
    "owls_REF" "owls_WDENS" "owls_WML1V848" "owls_WML4"
    "cowls_AGN" "cowls_AGN_T8p5" "cowls_AGN_T8p7"
    "BAHAMAS" "BAHAMAS_T7p6" "BAHAMAS_T8p0")

for sim in "${sims[@]}"; do
    ./like_fourier $1 0 $sim
done

echo Simulations calculated: ${sims[*]}
date
