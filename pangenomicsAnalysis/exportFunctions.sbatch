#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=32G
#SBATCH --partition=caslake
#SBATCH --output=estimatePan.out
#SBATCH --err=estimatePan.err

anvi-pan-genome -g genomeStorage-GENOMES.db \
                --project-name siVsLiPangenomics \
                --num-threads 12