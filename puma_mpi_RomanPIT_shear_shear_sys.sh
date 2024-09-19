#!/bin/bash -l
#SBATCH --job-name=PIT
#SBATCH --output=/xdisk/timeifler/jiachuanxu/job_logs/RomanPIT-%A_%a.out

### 1. puma
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --partition=high_priority
#SBATCH --qos=user_qos_timeifler
#SBATCH --account=timeifler
#SBATCH --time=96:00:00

### 2. ocelote
###SBATCH --nodes=10
###SBATCH --ntasks-per-node=28
###SBATCH --cpus-per-task=1
###SBATCH --partition=standard
###SBATCH --qos=qual_qos_timeifler
###SBATCH --account=timeifler
###SBATCH --time=9:00:00

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jiachuanxu@arizona.edu

echo "SLURM_ARRAY_TASK_ID = ${SLURM_ARRAY_TASK_ID}"
echo "SLURM_NTASKS = ${SLURM_NTASKS}"

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
conda activate python2_ext
which mpiexec
which python
MPIEXEC="/opt/ohpc/pub/mpi/mpich-gnu8-ohpc/3.3.1/bin/mpiexec"
#export MPI_DSM_DISTRIBUTE
date

# hit = 0
# for i in $(seq 0 4)
# do
# 	for j in $(seq 0 3)
# 	do
# 		for k in $(seq 0 1)
# 		do
# 			if [ ${hit} -eq ${SLURM_ARRAY_TASK_ID} ]; then
# 				mpiexec -n ${MPI_NPROCESS} python runRomanPIT_shear_shear.py ${i} ${j} ${k}
# 			hit = $(( ${hit} + 1 ))
# 		done
# 	done
# done

hit=0
for i in $(seq 0 4)
do
	for k in $(seq 0 1)
	do
		#hit = $(( ${hit} + 1 ))
		((hit++))
		if [ "$hit" -eq "$SLURM_ARRAY_TASK_ID" ]; then
			${MPIEXEC} -n ${SLURM_NTASKS} python runRomanPIT_shear_shear.py ${i} 2 ${k}
		fi
	done
done

date
