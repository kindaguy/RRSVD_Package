#/bin/bash

#Change the walltime according to a rough estimate of the running time of the program.

#Change walltime if needed
#PBS -l walltime=0:10:00
#PBS -l select=1:ncpus=1:ngpus=1:mem=10GB

#PBS -o log_program_run.out

#Change this according to what you want to do:
# debug: < 30 minutes
# parallel < 6 hours
# longpar < 24 hours

#PBS -q debug

#Send a message to the user when the job  (a) is aborted (b) begins (e) ends
#PBS -m abe

#Replace this address with yours.

#PBS -M tamascelli@di.unimi.it

#PBS -A LI03p_SURGREEN_0

module load gnu 
module load mkl
module load intel
module load cuda

echo PBS: working directory is $PBS_O_WORKDIR

cd ${HOME}/RRSVD_CUDA_Release_20150317/Examples

make exec

cd ..

