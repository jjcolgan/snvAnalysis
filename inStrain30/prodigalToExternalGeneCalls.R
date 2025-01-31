library(seqinr)
library(tidyverse)
fasta = read.fasta('inStrain30/genes.faa',as.string = T, seqtype = 'AA')
output = data.frame('gene_callers_id' = character(),
                 'contig' = character(),
                 'start' = numeric(),
                 'stop' = numeric(),
                 'direction' = character(),
                 'partial' = character(),
                 'call_type' = numeric(),
                 'source' = character(),
                 'version' = character())
for (i in 1:length(fasta)){
  temp = fasta[i]
  #list containing elements needed to create the the columns gene caller id, contig, start, stop, and strand
  nameStartStopStrand=getAnnot(temp)%>%
    str_split(pattern = '#')%>%
    .[1]
  name = nameStartStopStrand[[1]][1]
  name = str_split(name, '_')
  geneCaller = name[[1]][3]
  #this can be hard coded to just be 'c_000000000001
  contig = 'c_000000000001'
  start = nameStartStopStrand[[1]][2]
  #anvio is python and base 0 cause biology is dumb ig
  start = as.numeric(start) -1
  stop = nameStartStopStrand[[1]][3]
  stop = as.numeric(stop)-1
  'Put strand in anvio format'
  strand = nameStartStopStrand[[1]][4]
  if (strand == ' 1 ') {
    strand <- 'f'
  } else {
    strand <- 'r'
  }
  #gucci, now just need the partial / non partial, extract the ; delimited part
  genBankFormatInfo = nameStartStopStrand[[1]][5]
  genBankFormatInfo = genBankFormatInfo%>%
    str_split(pattern = ';')
  partial = genBankFormatInfo[[1]][2]
  if (grepl(pattern = '00', partial) == T){
    partial = '0'
  }else{
    partial = 1
  }

  #ok put the whole thing together
  row = data.frame('gene_callers_id' = geneCaller,
                   'contig' = contig,
                   'start' = start,
                   'stop' = stop,
                   'direction' = strand,
                   'partial' = partial,
                   'call_type' = 1,
                   'source' = 'prodigal',
                   'version' = 'V2.6.3')
  output = rbind(output, row)
}
write_tsv(output, 'externalGeneCalls.tsv')



