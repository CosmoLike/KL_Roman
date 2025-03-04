#!/bin/bash 
#SBATCH --job-name=SKA_KL_cn
#SBATCH --output=/xdisk/timeifler/yhhuang/log/cnSKA_KL-%A_%a.out
#SBATCH --error=log/cnSKA_KL-%A_%a.err
#SBATCH --nodes=1
#SBATCH --array=1-10
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
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
export LD_LIBRARY_PATH=/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`

export MPI_DSM_DISTRIBUTE
date
/usr/bin/time mpiexec -n 20 python runSKA_shear_shear_sys_KL.py
#/usr/bin/time mpiexec -n 10 python test.py
date
