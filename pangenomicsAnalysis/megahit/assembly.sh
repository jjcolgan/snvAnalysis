#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=32G
#SBATCH --partition=blekhman
#SBATCH --output=assembly.out
#SBATCH --err=assembly.err
mkdir 02_ASSEMBLY
cat samples.txt | parallel -j 4 "megahit -1 01_QC/{}_dehosted_R1.fastq.gz \
    -2 01_QC/{}_dehosted_R2.fastq.gz \
    -t 12 \
    --min-contig-len 750 \
    -o 02_ASSEMBLY/{} > 02_ASSEMBLY/{}.out 2> 02_ASSEMBLY/{}.err"