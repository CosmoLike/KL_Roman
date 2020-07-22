#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q standard
#PBS -l select=1:ncpus=1:mem=12GB
#PBS -l place=free:shared
#PBS -l walltime=5:00:00
#PBS -N KL_W1st_30_cov
#PBS -e /home/u17/jiachuanxu/output/
#PBS -o /home/u17/jiachuanxu/output/

module load gsl/2/2.1

cd $PBS_O_WORKDIR
python inv_cov.py
#/home/u17/jiachuanxu/CosmoLike/KL_WFIRST/./compute_covariances_fourier $PBS_ARRAY_INDEX




