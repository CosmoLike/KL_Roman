#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q standard
#PBS -J 1-1540
#PBS -l select=1:ncpus=1:mem=6GB
#PBS -l place=free:shared
#PBS -l walltime=6:00:00
#PBS -N KL_W1st_30_cov
#PBS -e /home/u17/jiachuanxu/output/
#PBS -o /home/u17/jiachuanxu/output/

module load gsl/2/2.1

cd $PBS_O_WORKDIR
#for (( c=0; c<233; c++ ))
#do
#	hit=$(( $PBS_ARRAY_INDEX + $c * 465 ))
#	/home/u17/jiachuanxu/CosmoLike/KL_WFIRST/./compute_covariances_fourier $hit
#done
/home/u17/jiachuanxu/CosmoLike/KL_WFIRST/./compute_covariances_fourier $PBS_ARRAY_INDEX




