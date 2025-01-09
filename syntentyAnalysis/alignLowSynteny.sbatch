#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --job-name=alignLowSynt
#SBATCH --err=alignLowSynteny.err
#SBATCH --output=alignLowSynteny.out
mkdir 04_ALIGNED
for sample in `cat samples.txt`
do
    bowtie2 -x 03_INDEX/lowSynteny \
    -1 /project/blekhman/jjcolgan/snv_data/cohort1RerurnCohort2Pangenomics/01_QC/${sample}_dehosted_R1.fastq.gz \
    -2 /project/blekhman/jjcolgan/snv_data/cohort1RerurnCohort2Pangenomics/01_QC/${sample}_dehosted_R2.fastq.gz \
    -p 4 -S 04_ALIGNED/${sample}.sam

    samtools view -bS 04_ALIGNED/${sample}.sam > 04_ALIGNED/${sample}.bam
    rm 04_ALIGNED/${sample}.sam

    samtools view -b -F 4 -q 1 04_ALIGNED/${sample}.bam > 04_ALIGNED/${sample}Filtered.bam
    rm 04_ALIGNED/${sample}.bam

    samtools sort 04_ALIGNED/${sample}Filtered.bam -o 04_ALIGNED/${sample}Sorted.bam
    rm 04_ALIGNED/${sample}Filtered.bam

    samtools index -@ 4 04_ALIGNED/${sample}Sorted.bam

done