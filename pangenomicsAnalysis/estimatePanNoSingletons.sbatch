#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=32G
#SBATCH --partition=caslake
#SBATCH --output=estimatePan.out
#SBATCH --err=estimatePan.err

anvi-pan-genome -g genomeStorage-GENOMES.db \
                --project-name siVsLiPangenomicsNoSingletons \
                --num-threads 12 \
                --mcl-inflation 10 \
                -T 12 \
                --min-occurrence 2
