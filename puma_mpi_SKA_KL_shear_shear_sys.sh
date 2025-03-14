#!/bin/bash 
#SBATCH --job-name=SKA_KL_cn
#SBATCH --output=/xdisk/timeifler/yhhuang/log/cnSKA_KL-%A.out
#SBATCH --error=log/cnSKA_KL-%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=80
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
module load openmpi5
conda init bash
source ~/.bashrc

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
cd $SLURM_SUBMIT_DIR
conda activate forecast-puma

export MPI_DSM_DISTRIBUTE
export LD_LIBRARY_PATH=/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`
date
#/usr/bin/time mpiexec -n 80 python runSKA_shear_shear_sys_KL.py
mpirun --mca pml ob1 --mca btl tcp,self -n 80 python runSKA_shear_shear_sys_KL.py
date
