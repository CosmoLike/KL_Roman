#!/bin/bash
#PBS -S /bin/bash
#PBS -W group_list=cosmo
#PBS -q high_pri
### Set the number of nodes,cores and memory that will be used for this job
### select=3 is the node count, ncpus=28 are the cores in each node,
### mem=168gb is memory per node, pcmem=6gb is the memory per core - optional
#PBS -l select=20:ncpus=28:mem=168GB
#PBS -l place=free:shared
#PBS -l cput=2800:00:00
#PBS -l walltime=5:00:00
#PBS -N WF_KL_ss_cos
#PBS -e /home/u17/jiachuanxu/output/
#PBS -o /home/u17/jiachuanxu/output/
#PBS -m bea
#PBS -M jiachuanxu@email.arizona.edu

### cput 8400:00:00 std
### walltime 15:00:00 std
cd $PBS_O_WORKDIR

module load python/2
module load mpich/ge/gcc/64/3.2.1
module load openmpi
VIRTUAL_ENV="/home/u17/jiachuanxu/python2_virtualenv"
export VIRTUAL_ENV

_OLD_VIRTUAL_PATH="$PATH"
PATH="$VIRTUAL_ENV/bin:$PATH"
export PATH
### run your executable program with begin and end date and time output
export MPI_DSM_DISTRIBUTE
echo $PBS_JOBNAME
date
#/usr/bin/time mpiexec -n 560 python runWFIRST_shear_shear_sys_opti.py
/usr/bin/time mpiexec -n 560 python runWFIRST_shear_shear_sys_KL.py
date
#echo "Your job $PBS_JOBID $PBS_JOBNAME is finished!" | mail -s "Your job $PBS_JOBID $PBS_JOBNAME is finished!" jiachuanxu@email.arizona.edu
