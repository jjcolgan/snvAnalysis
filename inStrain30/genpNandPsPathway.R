library(seqinr)      # For sequence operations (e.g., read.fasta, translate)
library(tidyverse)   # For data manipulation (e.g., read_tsv, filter)

# Function to calculate synonymous (sSites) and non-synonymous (nSites) sites for a given codon.
calculateCodonSandNSites = function(codon) {
  nucleotides = c('a', 't', 'g', 'c')  # List of possible nucleotides
  # Split the codon string into individual nucleotides
  codon <- strsplit(codon, split = "")[[1]]
  tempCodon = codon                 # Create a copy for mutation testing
  translatedCodon = translate(codon)  # Translate the original codon to its amino acid

  'this can remain the same for the entire gene'
  # NOTE: The above line is a string literal, not a comment. To comment it out, use #.

  nSites = 0  # Initialize counter for non-synonymous sites
  sSites = 0  # Initialize counter for synonymous sites

  # Loop over each position in the codon
  for (i in 1:length(codon)) {
    # Try substituting each nucleotide at position i
    for (n in nucleotides) {
      if (codon[i] != n) {           # Only change if the nucleotide is different
        tempCodon[i] = n             # Substitute nucleotide at position i
        mutated = translate(tempCodon)  # Translate the mutated codon
        tempCodon = codon            # Reset tempCodon back to the original codon

        # Compare the mutated amino acid to the original amino acid
        if (mutated != translatedCodon) {
          nSites = nSites + 1        # Count as non-synonymous change
        } else {
          sSites = sSites + 1        # Count as synonymous change
        }
      }
    }
  }
  # Return a vector: first element is synonymous sites, second is non-synonymous sites.
  return(c(sSites, nSites))
}

# Function to calculate the total number of non-synonymous and synonymous sites for a gene sequence.
calcGeneNandSSites = function(sequence) {
  start = 1
  stop = 3
  nSitesGene = 0  # Total non-synonymous sites for the gene
  sSitesGene = 0  # Total synonymous sites for the gene

  # Loop over the sequence codon-by-codon.
  # NOTE: This loop does not check if the sequence length is a multiple of 3.
  while (start < nchar(sequence)) {
    codonIn = substr(sequence, start, stop)  # Extract a codon from the sequence
    start = start + 3
    stop = stop + 3
    pAndNSites = calculateCodonSandNSites(codon = codonIn)
    # pAndNSites[1] is synonymous sites; pAndNSites[2] is non-synonymous sites.
    nSitesGene = pAndNSites[2] + nSitesGene
    sSitesGene = pAndNSites[1] + sSitesGene
  }
  # Return a vector: first element is total non-synonymous sites, second is total synonymous sites.
  return(c(nSitesGene, sSitesGene))
}

# Read pathways data from a TSV file.
# NOTE: The first read_tsv call is redundant since it is overwritten by the second call.
pathways = read_tsv('inStrain30/metabolism_modules.tsv')
pathways <- read_tsv('inStrain30/metabolism_modules.tsv',
                     col_types = cols(
                       .default = col_guess(),                    # Automatically guess types for other columns
                       gene_caller_ids_in_module = col_character()  # Explicitly set gene_caller_ids_in_module as character
                     ))

# Ensure the gene_caller_ids_in_module column is of type character.
pathways$gene_caller_ids_in_module = as.character(pathways$gene_caller_ids_in_module)

# Filter for pathways where either pathwise_module_is_complete or stepwise_module_is_complete is TRUE.
pathwaysComplete = pathways %>%
  filter(pathwise_module_is_complete == T | stepwise_module_is_complete == T)

# Read gene sequences from a FASTA file.
genes = read.fasta('genes.fna')

# Create an empty data frame to store pathway names and their summed synonymous and non-synonymous sites.
pathwaySandNsites = data.frame('pathway' = character(),
                               'pathwaySsites' = numeric(),
                               'pathwayNsites' = numeric())

# Loop over each complete pathway.
for (p in 1:nrow(pathwaysComplete)) {
  temp = pathwaysComplete[p, ]
  enzymes = temp$gene_caller_ids_in_module
  # Split the enzyme IDs by comma into a vector.
  enzymes = str_split_1(enzymes, pattern = ',')

  pathwayNSites = 0  # Initialize pathway non-synonymous site counter
  pathwaySSites = 0  # Initialize pathway synonymous site counter

  # Loop over each enzyme in the pathway.
  for (i in 1:length(enzymes)) {
    # Construct the full gene ID by appending a prefix.
    enzymes[i] = paste0('c_000000000001_', enzymes[i])

    # Retrieve the gene sequence vector from the genes list.
    g = genes[[enzymes[i]]]
    # Collapse the gene sequence vector into a single string.
    gene_seq_string <- paste(g, collapse = "")

    # Calculate the gene's non-synonymous and synonymous sites.
    geneNandSSites = calcGeneNandSSites(gene_seq_string)

    # geneNandSSites[1] is non-synonymous; geneNandSSites[2] is synonymous.
    pathwayNSites = pathwayNSites + geneNandSSites[1]
    pathwaySSites = pathwaySSites + geneNandSSites[2]
  }

  # Append the pathway's data to the results data frame.
  pathwaySandNsites = rbind(pathwaySandNsites, data.frame('pathway' = temp$module_name,
                                                          'pathwaySsites' = pathwaySSites,
                                                          'pathwayNsites' = pathwayNSites))
}



