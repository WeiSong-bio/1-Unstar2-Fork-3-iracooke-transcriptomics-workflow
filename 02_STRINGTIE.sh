#!/bin/bash
#PBS -c s
#PBS -j oe
#Send a mail alert at abortion, beginning and end of execution
#PBS -m abe
#Set the name of the job
#PBS -N stringtie
#Send any mail to this address
#PBS -M lauren.taylor2@my.jcu.edu.au
#Allocate required amount of wall time
#PBS -l walltime=1000:00:00
#Set the number of nodes and processors
#PBS -l nodes=1:ppn=4
#Allocate required amount of memory
#PBS -l pmem=4gb

f=DUMMY

cd /shares/32/jc320251/untouched_data

/sw/stringtie/1.3.1c/stringtie $f -G ../Homo_sapiens.GRCh38.87.gtf -o ${f%.bam}stringtie.gtf -p 4 -v -l ${f%_R1.bam}
