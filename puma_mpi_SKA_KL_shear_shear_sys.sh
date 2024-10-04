#!/bin/bash 
#SBATCH --job-name=SKA_KL_cn
#SBATCH --output=/xdisk/timeifler/yhhuang/log/cnSKA_KL-%A_%a.out
#SBATCH --error=log/cnSKA_KL-%A_%a.err
#SBATCH --nodes=1
#SBATCH --array=1-12
#SBATCH --ntasks-per-node=80
#SBATCH --ntasks-per-socket=40
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

module load anaconda
module load gsl
module swap openmpi3 mpich/3.3.1
conda init bash
source ~/.bashrc

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
cd $SLURM_SUBMIT_DIR
conda activate forecast

export MPI_DSM_DISTRIBUTE
date
/usr/bin/time mpiexec -n 400 python runSKA_shear_shear_sys_KL.py
date
