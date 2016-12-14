#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m ae
#PBS -N test
#PBS -M lauren.taylor2@my.jcu.edu.au
#PBS -l walltime=1000:00:00
#PBS -l nodes=1:ppn=8
#PBS -l pmem=3gb

ncpu=`wc -l $PBS_NODEFILE`
echo "------------------------------------------------------"
echo " This job is allocated "$ncpu" CPU cores on "
cat $PBS_NODEFILE | uniq
echo "------------------------------------------------------"
echo "PBS: Submitted to $PBS_QUEUE@$PBS_O_HOST"
echo "PBS: Working directory is $PBS_O_WORKDIR"
echo "PBS: Job identifier is $PBS_JOBID"
echo "PBS: Job name is $PBS_JOBNAME"
echo "------------------------------------------------------"
