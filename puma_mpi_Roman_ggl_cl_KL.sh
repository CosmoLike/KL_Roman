#!/bin/bash -l
#SBATCH --job-name=2x2pt
#SBATCH --output=/xdisk/timeifler/yhhuang/log/2x2pt-%A_%a.out
#SBATCH --error=log/2x2pt-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=80
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yhhuang@arizona.edu

echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_NTASKS = ${SLURM_NTASKS}"

module load gsl
module load openmpi5
module load anaconda
conda init bash
source ~/.bashrc

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Slurm job NAME is $SLURM_JOB_NAME
echo Slurm job ID is $SLURM_JOB_ID
cd $SLURM_SUBMIT_DIR

conda activate forecast-puma
export MPI_DSM_DISTRIBUTE
export LD_LIBRARY_PATH=/opt/ohpc/pub/libs/gnu8/gsl/2.6/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:`pwd`

date
mpirun --mca pml ob1 --mca btl tcp,self -n 80 python runRoman_ggl_cl_KL.py -nsteps=4000 -nwalkers=400
date
