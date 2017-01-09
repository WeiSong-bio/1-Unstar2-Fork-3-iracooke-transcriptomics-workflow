#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m abe
#PBS -N stringtie
#PBS -M username@my.jcu.edu.au
#PBS -l walltime=1000:00:00
#PBS -l nodes=1:ppn=4
#PBS -l pmem=4gb

f=DUMMY

cd /shares/32/jc320251/untouched_data

/sw/stringtie/1.3.1c/stringtie -e -B -p 4 -G stringtie_merged.gtf -o ballgown/${f%_R1.bam}/${f%_R1.bam}.gtf $f
