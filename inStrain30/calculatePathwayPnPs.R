library(tidyverse)
library(ggrepel)
pathwayPnPs = read_tsv('pathwaySandNCounts.tsv')
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
samples <- read_tsv('inStrain30/samples.txt', col_names = FALSE)
snvs<- read_tsv(
  paste0('inStrain30/',samples$X1[1],'_SNVs.tsv'))
snvs$sample <- samples$X1[1]

#metadata = read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData/cohort1and2metadata.txt')


for (i in 2:length(samples$X1)){
  print(samples$X1[i])
  temp<- read_tsv(
    paste0('inStrain30/',samples$X1[i],'_SNVs.tsv'))
  temp$sample <- samples$X1[i]
  snvs <- rbind(snvs, temp)

}

snvs$tissue = 'co'
snvs[grepl(snvs$sample, pattern = 'dj', ignore.case = T),]$tissue = 'dj'
snvs$cohort = '2'
snvs[grepl(snvs$sample, pattern ='john', ignore.case = T),]$cohort = '1'

# Create an empty data frame to store output information.
output = data.frame('sample' = character(),
                    'pathway' = character(),
                    'nSites' = numeric(),
                    'sSites' = numeric())

# Loop over each sample (make sure the variable 'samples' is defined beforehand).
output = data.frame('sample' = character(),
                    'pathway' = character(),
                    'nSites' = numeric(),
                    'sSites' = numeric())

for (s in samples) {

  for (i in 1:nrow(pathwaysComplete)) {
    temp = pathwaysComplete[i, ]
    enzymes = temp$gene_caller_ids_in_module
    # Split the enzyme IDs by comma into a vector.
    enzymes = str_split_1(enzymes, pattern = ',')

    # Use a different loop variable here.
    for (j in 1:length(enzymes)) {
      enzymes[j] = paste0('c_000000000001_', enzymes[j])
    }

    sampleSNVs = snvs %>%
      filter(sample == s) %>%
      filter(gene %in% enzymes)

    nSample = sampleSNVs %>%
      filter(mutation_type == 'N') %>%
      nrow()

    sSample = sampleSNVs %>%
      filter(mutation_type == 'S') %>%
      nrow()

    output = rbind(output,
                   data.frame('sample' = s,
                              'pathway' = temp$module_name[1],
                              'nSites' = nSample,
                              'sSites' = sSample))
  }
}



output = data.frame('sample' = character(),
                    'pathway' = character(),
                    'nSites' = numeric(),
                    'sSites' = numeric(),
                    'pathwaySnvs' = numeric(),
                    'samplePn' = numeric(),
                    'samplePs' = numeric(),
                    'pathwayPnPs' = numeric())
for (s in 1:nrow(samples)){
  for(i in 1:nrow(pathwaysComplete)){
    temp = pathwaysComplete[i, ]
    enzymes = temp$gene_caller_ids_in_module
    # Split the enzyme IDs by comma into a vector.
    enzymes = str_split_1(enzymes, pattern = ',')
    for( i in 1:length(enzymes)){
      enzymes[i] = paste0('c_000000000001_', enzymes[i])
    }
    sampleSNVs = snvs %>%
      filter(sample == samples$X1[s])%>%
      filter(gene %in% enzymes)
    nSample = sampleSNVs %>%
      filter(mutation_type == 'N')%>%
      nrow()
    sSample = sampleSNVs %>%
      filter(mutation_type == 'S')%>%
      nrow()
    if (nSample > 0 & sSample < 1){
      sSample = 1
    }
    pathwayPnPsTemp = pathwayPnPs%>%
      filter(pathway == temp$module_name[1])
    pNSample = nSample / pathwayPnPsTemp$pathwayNsites
    pSSample = sSample / pathwayPnPsTemp$pathwaySsites
    output= rbind(output,
                  data.frame('sample' = samples$X1[s],
                             'pathway' = temp$module_name[1],
                             'nSites' = nSample,
                             'sSites' = sSample,
                             'pathwaySnvs' = nSample+sSample,
                             'samplePn' = pNSample,
                             'samplePs' = pSSample,
                             'pathwayPnPs' = pNSample/pSSample))
  }
}

output %>%
  pivot_wider(names_from = 'pathway',
              id_cols = 'sample',
              values_from = 'pathwayPnPs')%>%
  view()
output$tissue = 'co'
output[grepl(output$sample, pattern = 'dj', ignore.case = T),]$tissue = 'dj'
meanPathwaypNpS=output %>%
  group_by(tissue, pathway)%>%
  summarise('meanTissuePnPS' = mean(pathwayPnPs, na.rm = T))


meanPathwaypNpSWide=meanPathwaypNpS%>%
  pivot_wider(names_from = tissue, id_cols = pathway, values_from = meanTissuePnPS)

meanPathwaypNpSWide$selection = 'negative'
meanPathwaypNpSWide$selection[meanPathwaypNpSWide$dj > 1 & meanPathwaypNpSWide$co < 1] = 'adaptive dj'
meanPathwaypNpSWide$selection[meanPathwaypNpSWide$dj < 1 & meanPathwaypNpSWide$co > 1] = 'adaptive co'
meanPathwaypNpSWide$selection[meanPathwaypNpSWide$dj > 1 & meanPathwaypNpSWide$co > 1] = 'adaptive both'
meanPathwaypNpSWide$label = NA
meanPathwaypNpSWide$delta = abs(meanPathwaypNpSWide$co-meanPathwaypNpSWide$dj)
meanPathwaypNpSWide$label[meanPathwaypNpSWide$delta > 0.3 ] <- meanPathwaypNpSWide$pathway[meanPathwaypNpSWide$delta > 0.3]
meanPathwaypNpSWide%>%
  ggplot(aes(x = dj,
             y = co,
             label = label,
             col = selection))+
  geom_point(alpha = .5)+
  geom_text_repel(size = 2)

meanPathwaypNpSWide %>%
  filter(selection == 'adaptive dj')%>%
  arrange(desc(delta))%>%
  head(5)

meanPathwaypNpSWide %>%
  filter(selection == 'adaptive co')%>%
  arrange(desc(delta))

meanPathwaypNpSWide %>%
  filter(selection == 'negative')

meanPathwaypNpSWide %>%
  filter(selection == 'adaptive both')


outputWide=output %>%
  pivot_wider(names_from = 'pathway',
              id_cols = 'sample',
              values_from = 'pathwayPnPs')
outputWide[is.na(outputWide)]=-1

pcaOut=outputWide %>%
  column_to_rownames('sample')%>%
  prcomp(center = T,scale = T)

metadata=output %>%
  select(sample, tissue)%>%
  distinct()

summary(pcaOut)
pcaOut$x %>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(metadata, by = 'sample')%>%
  ggplot(aes(x = PC2,
             col = tissue,
             y =PC3))+
  geom_point(alpha = .25)