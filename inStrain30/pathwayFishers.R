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
results = data.frame('pathway' = character(),
                     'pvalue' = double())
for (p in 1:nrow(pathwaysComplete)){
  temp=pathwaysComplete[p,]
  enzymes = temp$gene_caller_ids_in_module
  enzymes = str_split_1(enzymes, pattern = ',')
  for (i in 1:length(enzymes)){
    enzymes[i] = paste0('c_000000000001_', enzymes[i])
  }
#store the pnps scores for the gene in the pathway
  contigencyTable = data.frame('tissue' = character(),
                             'above1'= integer(),
                             'below1' = integer())
  tissues = unique(genesInfo$tissue)
  for (t in 1:length(tissues)){
    pathwayPnPs=genesInfo%>%
    filter(tissue == tissues[t],
         gene %in% enzymes)
    nGenesAbove1 = pathwayPnPs %>%
    filter(is.na(pNpS_variants)==F | SNV_N_count >0)%>%
    filter(pNpS_variants > 1)%>%
    nrow()

    nGenesBelow1 = pathwayPnPs %>%
    filter(is.na(pNpS_variants)==F)%>%
    filter(pNpS_variants < 1)%>%
    nrow()
    row = data.frame('tissue' = tissues[t],
                    'above1'= nGenesAbove1,
                   'below1' = nGenesBelow1)
    contigencyTable = rbind(contigencyTable,
                          row)


}
  print(contigencyTable)
  fishersOut= contigencyTable %>%
  column_to_rownames('tissue')%>%
  fisher.test()
  resAdd = data.frame('pathway' = pathwaysComplete[p,]$module,
                    'pvalue' = fishersOut$p.value)
  results = rbind(results,resAdd)
}
