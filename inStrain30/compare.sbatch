#!/bin/bash
#SBATCH --job-name=compare
#SBATCH --output=compare.out
#SBATCH --error=compare.err
#SBATCH --nodes=1 --account=pi-blekhman --partition=blekhman
#SBATCH --ntasks=6
#SBATCH --mem=60000

inStrain compare -i 04_instrain_minReads3_minFPoint005/20240801_Changlab_SSG_2678CoMGS_S2 \
04_instrain_minReads3_minFPoint005/20240801_Chang_lab_SSG_2678_DJ_MGS_S23 \
04_instrain_minReads3_minFPoint005/20240801_Changlab_SSG_2852CoMGS_S5 \
04_instrain_minReads3_minFPoint005/20240801_Chang_lab_SSG_2852_DJ_MGS_S25 \
04_instrain_minReads3_minFPoint005/20240801_Chang_lab_SSG_2868_Co_MGS_S26 \
04_instrain_minReads3_minFPoint005/20240801_Changlab_SSG_COMGSNT_4_S1 \
04_instrain_minReads3_minFPoint005/20240801_Chang_lab_SSG_CO_MGsS_8455_S24 \
04_instrain_minReads3_minFPoint005/20240801_Chang_lab_SSG_DJ_MGS_8455_S22 \
04_instrain_minReads3_minFPoint005/20240801_Chang_lab_SSG_DJ_MGS_NT_4_S21 \
04_instrain_minReads3_minFPoint005/20240801_Changlab_SSG_NT_3CoMGS_S3 \
04_instrain_minReads3_minFPoint005/20240801_Changlab_SSG_NT_3DJMGS_S4 \
04_instrain_minReads3_minFPoint005/240523-John-2439-co \
04_instrain_minReads3_minFPoint005/240523-John-2439-dj \
04_instrain_minReads3_minFPoint005/240523-John-2444-co \
04_instrain_minReads3_minFPoint005/240523-John-2444-dj \
04_instrain_minReads3_minFPoint005/240523-John-2445-co \
04_instrain_minReads3_minFPoint005/240523-John-2445-dj \
04_instrain_minReads3_minFPoint005/240523-John-2446-co \
04_instrain_minReads3_minFPoint005/240523-John-2446-dj \
04_instrain_minReads3_minFPoint005/240523-John-2447-dj \
04_instrain_minReads3_minFPoint005/240523-John-2448-co \
04_instrain_minReads3_minFPoint005/240523-John-2448-dj \
04_instrain_minReads3_minFPoint005/240523-John-2476-co \
04_instrain_minReads3_minFPoint005/240523-John-2476-dj \
-f .0005 \
--bams 01_QC/20240801_Chang_lab_SSG_2678_DJ_MGS_S23Aligned.sorted.bam \
01_QC/20240801_Changlab_SSG_2852CoMGS_S5Aligned.sorted.bam \
01_QC/20240801_Chang_lab_SSG_2852_DJ_MGS_S25Aligned.sorted.bam \
01_QC/20240801_Chang_lab_SSG_2868_Co_MGS_S26Aligned.sorted.bam \
01_QC/20240801_Changlab_SSG_COMGSNT_4_S1Aligned.sorted.bam \
01_QC/20240801_Chang_lab_SSG_CO_MGsS_8455_S24Aligned.sorted.bam \
01_QC/20240801_Chang_lab_SSG_DJ_MGS_8455_S22Aligned.sorted.bam \
01_QC/20240801_Chang_lab_SSG_DJ_MGS_NT_4_S21Aligned.sorted.bam \
01_QC/20240801_Changlab_SSG_NT_3CoMGS_S3Aligned.sorted.bam \
01_QC/20240801_Changlab_SSG_NT_3DJMGS_S4Aligned.sorted.bam \
01_QC/240523-John-2439-coAligned.sorted.bam \
01_QC/240523-John-2439-djAligned.sorted.bam \
01_QC/240523-John-2444-coAligned.sorted.bam \
01_QC/240523-John-2444-djAligned.sorted.bam \
01_QC/240523-John-2445-coAligned.sorted.bam \
01_QC/240523-John-2445-djAligned.sorted.bam \
01_QC/240523-John-2446-coAligned.sorted.bam \
01_QC/240523-John-2446-djAligned.sorted.bam \
01_QC/240523-John-2447-djAligned.sorted.bam \
01_QC/240523-John-2448-coAligned.sorted.bam \
01_QC/240523-John-2448-djAligned.sorted.bam \
01_QC/240523-John-2476-coAligned.sorted.bam \
01_QC/240523-John-2476-djAligned.sorted.bam \
-p 8 
-o 05_compare_minReads3_minFPoint005
