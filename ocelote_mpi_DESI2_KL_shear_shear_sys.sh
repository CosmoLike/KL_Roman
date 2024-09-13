#!/bin/bash
#SBATCH --job-name=DESI2_KL_cn
#SBATCH --output=DESI2_KL-%A_%a.out
#SBATCH --nodes=10
#SBATCH --array=1-36
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=120:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu

iSELECT=$(( (${SLURM_ARRAY_TASK_ID}-1) / 6 ))
iSN=$(( (${SLURM_ARRAY_TASK_ID}-1) - ( ${iSELECT} * 6 ) ))

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
echo "Target Selection Scenario = ${iSELECT}"
echo "Shape Noise Scenario = ${iSN}"

module purge > /dev/null 2>&1
module load gsl
module load mpich
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
which mpiexec
export MPI_DSM_DISTRIBUTE
date
#/usr/bin/time mpiexec -n 560 python runWFIRST_shear_shear_sys_opti.py
/usr/bin/time mpiexec -n 400 python runWFIRST_shear_shear_sys_KL.py ${iSELECT} ${iSN}
date
