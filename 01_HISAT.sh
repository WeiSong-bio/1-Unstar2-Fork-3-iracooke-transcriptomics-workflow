#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m abe
#PBS -N hisat2
#PBS -M lauren.taylor2@my.jcu.edu.au
#PBS -l walltime=1000:00:00
#PBS -l nodes=1:ppn=4
#PBS -l pmem=4gb


f=DUMMY
cd /shares/32/jc320251/untouched_data
/sw/hisat2/2.0.5/hisat2 -p 4 -x ../grch38/genome -U $f -S ${f%.fastq.gz}.sam;
/sw/samtools/1.3/AMD/bin/samtools sort -@ 4 -o ${f%.fastq.gz}.bam ${f%.fastq.gz}.sam;
