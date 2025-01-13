input = open('samples.txt','r')
samples = input.read()
samples = samples.split('\n')
samples = samples[:-1]
rule all:
    input:
        expand('04_instrain/{sample}/done', sample = samples)
rule trim:
    resources:
        cpus_per_task = 1,
        mem_mb = 5000,
        tasks = 4,
        time = '2h',
        nodes = 1,
        account = 'pi-blekhman',
        partition = 'blekhman'
    threads :4
    input:
        R1 = '/project/blekhman/jjcolgan/snv_data/cohort1RerurnCohort2Pangenomics/01_QC/{sample}dehosted_R1.fastq.gz',
        R2= '/project/blekhman/jjcolgan/snv_data/cohort1RerurnCohort2Pangenomics/01_QC/{sample}dehosted_R2.fastq.gz'
    conda:'biobakery3'
    output:
        R1 = temp('01_QC/{sample}Dehosted_R1.fastq.gz'),
        R2 = temp('01_QC/{sample}Dehosted_R1.fastq.gz'),
        R1Orphaned = temp('01_QC/{sample}_dehosted_R1.fastq.gz'),
        R2Orphaned = temp('01_QC/{sample}_dehosted_R2.fastq.gz')
    log:
        out = '01_Q30_QC/{sample}Trim.out',
        err = '01_Q30_QC/{sample}Trim.err'
    shell:
        '''
        java -jar /project/blekhman/jjcolgan/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads {threads} -phred33 \
	    {input.R1} {input.R2} \
	    {output.R1} {output.12Orphaned} \
	    {output.R2} {output.R2Orphaned} \
	    ILLUMINACLIP:/project/blekhman/jjcolgan/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa/:2:30:10:8:true \
	    SLIDINGWINDOW:4:30 MINLEN:40 > {log.out} 2> {log.err}
        '''
rule alignRef:
    resources:
        cpus_per_task=1,
        mem_mb=5000,
        tasks=4,
        time='2h',
        nodes=1,
        account='pi-blekhman',
        partition='blekhman'
    threads:4
    conda:'biobakery3'
    input:
        R1 = '01_QC/{sample}Dehosted_R1.fastq.gz',
        R2 = '01_QC/{sample}Dehosted_R1.fastq.gz'
    output:
        sam = temp('01_QC/{sample}Aligned.sam')
    log:
        err = '01_QC/{sample}Aligned.err',
        out= '01_QC/{sample}Aligned.out'
    shell:
        '''
        bowtie2 -p {threads} \
        --no-unal \
        -x 03_INDEX/ref \
        -1 {input.R1} \
        -2 {input.R2} --very-sensitive \
        -S {output.sam} \
        2> {log.err} \
        > {log.out}
        '''
rule convertBam:
    resources:
        cpus_per_task=1,
        mem_mb=2000,
        tasks=1,
        time='2h',
        nodes=1,
        account='pi-blekhman',
        partition='blekhman'
    input:
        sam = '01_QC/{sample}Aligned.sam'
    output:
        bam = temp('01_QC/{sample}Aligned.bam')

    shell:
        '''
        samtools view -bS {input.sam} > {output.bam}
        '''
rule profile:
    resources:
        cpus_per_task=1,
        mem_mb=30000,
        tasks=8,
        time='2h',
        nodes=1,
        account='pi-blekhman',
        partition='blekhman'
    conda: 'instrain'
    input:
       bam= '01_QC/{sample}Aligned.bam'
    output:
        done ='04_instrain/{sample}/done'
    params:
        dir = '04_instrain/{sample}/'
    threads:8
    log:
        output = '04_instrain/{sample}.out',
        err = '04_instrain/{sample}.err'
    shell:
        '''
	    inStrain profile {input.bam} \
	    02_CONTIGS/reference.fasta \
	    -g 02_CONTIGS/genes.fna \
	    -p {threads} \
	    -f .05 \
	    -c 10 \
	    -o {params.dir} > {log.out} 2> {log.err}
	    touch {output.done}
        '''