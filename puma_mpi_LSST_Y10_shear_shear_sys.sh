#!/bin/bash
#SBATCH --job-name=LSST10cn
#SBATCH --output=LSST_Y10-%A_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --ntasks-per-socket=40
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"

module load gsl
module swap openmpi3 mpich/3.3.1
module load anaconda
conda init bash
source ~/.bashrc

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOBID
cd $SLURM_SUBMIT_DIR
conda activate python2

export MPI_DSM_DISTRIBUTE
date
#/usr/bin/time mpiexec -n 560 python runWFIRST_shear_shear_sys_opti.py
/usr/bin/time mpiexec -n 400 python runLSST_shear_shear_sys_Y10.py
date
