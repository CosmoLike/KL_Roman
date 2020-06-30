#!/bin/bash
#PBS -S /bin/bash
#PBS -V
#PBS -W group_list=cosmo
#PBS -q qualified
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l place=free:shared
#PBS -l walltime=1:00:00
#PBS -N KL_W1st_data_vec
#PBS -e /home/u17/jiachuanxu/output/
#PBS -o /home/u17/jiachuanxu/output/

module load gsl/2/2.1

cd $PBS_O_WORKDIR
/home/u17/jiachuanxu/CosmoLike/KL_WFIRST/./like_fourier opti WFIRST_KL shear_shear >& /home/u17/jiachuanxu/output/job_output_datavec.log




