#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=4
#SBATCH --mem=12G
#SBATCH --partition=caslake
#SBATCH --output=exportKegg.out
#SBATCH --err=exportKegg.err

for sample in `cat samples.txt`
do
    anvi-export-functions --annotation-sources KOfam \
    -o 05_CONTIGS_DB/{$sample}/{$sample}Kofams.tsv \
    -c 05_CONTIGS_DB/{$sample/contigs.db
done