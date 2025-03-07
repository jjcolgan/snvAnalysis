'''This is a wrapper script that submits instrain jobs to slurm and optimizes the parameters
to be identical to those used in https://www.nature.com/articles/nature11711, where snvs are called
if there are at least 4 reads supporting a call and their frequency is at least 1%. '''
import os
def calculateCandF(sampleCoverage):
    if sampleCoverage >= 400:
        c = 400
        f = 0.01
    else:
        c = str(sampleCoverage)
        f = str(4 / sampleCoverage)

    return c, f
def writeScript(c, f, sample):
    with open(sample+".sbatch", "w") as file:
        file.write("#!/bin/bash\n")
        file.write("#SBATCH --job-name="+sample+"\n")
        file.write("#SBATCH --output="+sample+".out\n")
        file.write("#SBATCH --error="+sample+".err\n")
        file.write("#SBATCH --nodes=1\n")
        file.write("#SBATCH --ntasks=8\n")
        file.write("#SBATCH --mem=30000")
        file.write("\n")
        file.write('inStrain profile 01_QC/'+sample+'Aligned.bam 02_CONTIGS/ref.fasta -g 02_CONTIGS/genes.fna -p 8 -f'+ f+' -c '+ c+'-o 04_instrain_schlossnig/'+sample+'/ > 04_instrain_schlossnig/'+sample+'.out 2> 04_instrain_schlossnig/'+sample+'.err')
        file.close()
        #os.system("sbatch "+sample+".sbatch")

inputFile = open('samples.txt', 'r')
samples = inputFile.read()
samples = samples.split('\n')
samples = samples[:-1]
inputFile.close()

inputFile = open('coverage.txt', 'r')
coverages = inputFile.read()
coverages = coverages.split('\n')
coverages = coverages[:-1]
inputFile.close()

for i in range(0, len(samples)):
    sampleCoverage = coverages[i]
    sample = samples[i]
    c, f = calculateCandF(sampleCoverage)
    writeScript(c, f, sample)
