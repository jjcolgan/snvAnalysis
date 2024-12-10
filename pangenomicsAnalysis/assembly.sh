#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=32G
#SBATCH --partition=blekhman
#SBATCH --output=assembly.out
#SBATCH --err=assembly.err
mkdir 02_ASSEMBLY
cat samples.txt | parallel -j 4 "megahit -1 01_QC/{}_R1_dehosted.fastq.gz \
    -2 01_QC/{}_R2_dehosted.fastq.gz \
    -t 12 \
    --min-contig-len 750 \
    -o 02_ASSEMBLY/{} > 02_ASSEMBLY/{}.out 2> 02_ASSEMBLY/{}.err"