#!/bin/bash
#SBATCH --account=pi-blekhman
#SBATCH --ntasks=12
#SBATCH --mem=32G
#SBATCH --partition=caslake
#SBATCH --output=synTracker.out
#SBATCH --err=synTracker.err

python syntracker.py -target 12_SYNTRACKER/assemblies/ \
  -ref 12_SYNTRACKER/reference \
  -out 12_SYNTRACKER/output -cores 12