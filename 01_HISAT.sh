#!/bin/bash
#PBS -c s
#PBS -j oe
#PBS -m abe
#PBS -N hisat2
#PBS -M lauren.taylor2@my.jcu.edu.au
#PBS -l walltime=1000:00:00
#PBS -l nodes=1:ppn=24
#PBS -l pmem=12gb

module load hisat2
module load samtools
for f in *.fastq.gz; do
  hisat2 -p 24 -x ../grch38/genome -1 $f -S ${f%.fastq.gz}.sam;
  samtools sort -@ 24 -o ${f%.fastq.gz}.bam ${f%.fastq.gz}.sam;
done
