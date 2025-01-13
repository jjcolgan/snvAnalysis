from Bio import SeqIO

sequences = open('syntentyAnalysis/sharedLowSyntenyRegions.tsv','r')
sequences = sequences.readlines()
fasta_file = "pangenomicsAnalysis/reference.fa"
fasta_out = 'syntentyAnalysis/lowSynTenyRegions.fasta'
sequences=sequences[1:]

for sequence in sequences:
    sequence = sequence.strip()
    sequence = sequence.split('\t')
    start = int(sequence[0])-2000
    stop = int(sequence[1])+2000

    with open(fasta_file, "r") as file:
        with open(fasta_out, "a") as output:  # Create or overwrite the output file
            for record in SeqIO.parse(file, "fasta"):
                # Extract the sequence from x to y
                new_id = record.id+'_'+str(start)+':'+str(stop)
                extracted_sequence = record.seq[start-1:stop]

                # Write the extracted sequence to the output file
                output.write(f">{new_id}_extracted\n{extracted_sequence}\n")