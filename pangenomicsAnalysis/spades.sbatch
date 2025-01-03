#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=40G
#SBATCH --partition=caslake
#SBATCH --output=spades/assembly.out
#SBATCH --err=spades/assembly.err

#used because genomes recovered with megahit assembly were low quality
mkdir spades/02_ASSEMBLY
for sample in `cat samples.txt`
do 
	mkdir spades/02_ASSEMBLY/$sample
	spades.py --meta -o spades/02_ASSEMBLY/$sample \
	-m 9 -t 12 \
	-1 01_QC/$sample_dehosted_R1.fastq.gz \
	-2 01_QC/$sample_dehosted_R2.fastq.gz > spades/02_ASSEMBLY/$sample/spades.out 2> spades/02_ASSEMBLY/$sample/spades.err
done
