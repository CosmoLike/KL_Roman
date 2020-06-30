#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q qualified
#PBS -J 1-1540
#PBS -l select=1:ncpus=1:mem=6GB
#PBS -l place=free:shared
#PBS -l walltime=1:00:00
#PBS -N WL_W1st_KL_cov
#PBS -e /home/u17/jiachuanxu/output/
#PBS -o /home/u17/jiachuanxu/output/

module load gsl/2/2.1
module load mpich/ge/gcc/64/3.2.1
module load openmpi

cd $PBS_O_WORKDIR
/home/u17/jiachuanxu/CosmoLike/KL_WFIRST/./compute_covariances_fourier $PBS_ARRAY_INDEX >& /home/u17/jiachuanxu/output/job_output_$PBS_ARRAY_INDEX.log




