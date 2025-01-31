library(tidyverse)
samples <- read_tsv('inStrain30/04_instrain/samples.txt', col_names = FALSE)
genomeInfo <- read_tsv(
paste0('inStrain30/04_instrain/04_instrain/',samples$X1[1],'/output/',samples$X1[1],'_genome_info.tsv'))
genomeInfo$sample <- samples$X1[1]
genesInfo<- read_tsv(
paste0('inStrain30/04_instrain/04_instrain/',samples$X1[1],'/output/',samples$X1[1],'_gene_info.tsv'))
genesInfo$sample <- samples$X1[1]


for (i in 2:length(samples$X1)){
temp <- read_tsv(
paste0('inStrain30/04_instrain/04_instrain/',samples$X1[i],'/output/',samples$X1[i],'_genome_info.tsv'))
temp$sample <- samples$X1[i]
genomeInfo<- rbind(genomeInfo, temp)
temp<- read_tsv(
paste0('inStrain30/04_instrain/04_instrain/',samples$X1[i],'/output/',samples$X1[i],'_gene_info.tsv'))
temp$sample <- samples$X1[i]
genesInfo <- rbind(genesInfo, temp)

}

genesInfo$tissue = 'co'
genesInfo[grepl(genesInfo$sample, pattern = 'dj', ignore.case = T),]$tissue = 'dj'
genesInfo$cohort = '2'
genesInfo[grepl(genesInfo$sample, pattern ='john', ignore.case = T),]$cohort = '1'

genomeInfo$tissue = 'co'
genomeInfo[grepl(genomeInfo$sample, pattern = 'dj', ignore.case = T),]$tissue = 'dj'
genomeInfo$cohort = '2'
genomeInfo[grepl(genomeInfo$sample, pattern ='john', ignore.case = T),]$cohort = '1'

genesInfo%>%
  filter(is.na(pNpS_variants)==F)%>%
  select(sample)%>%
  dim()

library(readr)

pathways <- read_tsv('inStrain30/metabolism_modules.tsv',
                     col_types = cols(
                       .default = col_guess(),   # Automatically guess other columns
                       gene_caller_ids_in_module = col_character()    # Set the 16th column as character
                     ))

pathways$gene_caller_ids_in_module = as.character(pathways$gene_caller_ids_in_module)
pathwaysComplete=pathways%>%
  filter(pathwise_module_is_complete == T | stepwise_module_is_complete == T)

temp=pathwaysComplete[1,]