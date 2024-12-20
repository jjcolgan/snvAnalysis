#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=32G
#SBATCH --partition=caslake
#SBATCH --output=metabolism.out
#SBATCH --err=metabolism.err

anvi-estimate-metabolism -M metagenomes.txt \
--kegg-data-dir /project/blekhman/shared/anvioDbs/kegg/ \
--matrix-format

