library(tidyverse)
samples <- read_tsv('inStrain30/samples.txt', col_names = FALSE)
genomeInfo <- read_tsv(
  paste0('inStrain30/',samples$X1[1],'_genome_info.tsv'))
genomeInfo$sample <- samples$X1[1]
genesInfo<- read_tsv(
  paste0('inStrain30/',samples$X1[1],'_gene_info.tsv'))
genesInfo$sample <- samples$X1[1]


for (i in 2:length(samples$X1)){
  print(samples$X1[i])
  temp <- read_tsv(
    paste0('inStrain30/',samples$X1[i],'_genome_info.tsv'))
  temp$sample <- samples$X1[i]
  genomeInfo<- rbind(genomeInfo, temp)
  temp<- read_tsv(
    paste0('inStrain30/',samples$X1[i],'_gene_info.tsv'))
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


meanTissuePnPS=genesInfo %>%
  filter(SNV_count >= 3 )%>%
  group_by(tissue, gene)%>%
  summarise('meanTissuePnPS' = mean(pNpS_variants, na.rm = T))


meanTissuePnPsWide=meanTissuePnPS%>%
  pivot_wider(names_from = tissue, id_cols = gene, values_from = meanTissuePnPS)

meanTissuePnPsWide$selection = 'negative'
meanTissuePnPsWide$selection[meanTissuePnPsWide$dj > 1 & meanTissuePnPsWide$co < 1] = 'adaptive dj'
meanTissuePnPsWide$selection[meanTissuePnPsWide$dj < 1 & meanTissuePnPsWide$co > 1] = 'adaptive co'
meanTissuePnPsWide$selection[meanTissuePnPsWide$dj > 1 & meanTissuePnPsWide$co > 1] = 'adaptive both'

meanTissuePnPsWide %>%
  na.omit()%>%
  ggplot(aes(x = selection))+
  geom_bar()

meanTissuePnPsWide%>%
  ggplot(aes(x = dj,
             y = co,
             col = selection))+
  geom_point(alpha = .5)+
  labs(title = 'mean pNpS no SNV count > 30')

meanTissuePnPsPrevelanceFiltered=genesInfo %>%
  filter(SNV_count >= 3) %>%
  group_by(sample, SNV_N_count, tissue) %>%  # Group by tissue
  mutate(pNpS_variants = ifelse(is.na(pNpS_variants), mean(pNpS_variants, na.rm = TRUE), pNpS_variants)) %>%
  pivot_wider(names_from = gene, id_cols = c(sample, tissue), values_from = pNpS_variants) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  # Group by tissue before filtering
  select(sample, tissue, where(~ mean(. == 0, na.rm = TRUE) <= 0.75)) %>%
  pivot_longer(c(3:ncol(.)), values_to = 'pNpS_variants', names_to = 'gene')%>%
  group_by(tissue, gene) %>%
  summarise('meanTissuePnPS' = mean(pNpS_variants, na.rm = T))

meanTissuePnPsPrevelanceFilteredWide=meanTissuePnPsPrevelanceFiltered%>%
  pivot_wider(names_from = tissue, id_cols = gene, values_from = meanTissuePnPS)

meanTissuePnPsPrevelanceFilteredWide$selection = 'negative'
meanTissuePnPsPrevelanceFilteredWide$selection[meanTissuePnPsPrevelanceFilteredWide$dj > 1 & meanTissuePnPsPrevelanceFilteredWide$co < 1] = 'adaptive dj'
meanTissuePnPsPrevelanceFilteredWide$selection[meanTissuePnPsPrevelanceFilteredWide$dj < 1 & meanTissuePnPsPrevelanceFilteredWide$co > 1] = 'adaptive co'
meanTissuePnPsPrevelanceFilteredWide$selection[meanTissuePnPsPrevelanceFilteredWide$dj > 1 & meanTissuePnPsPrevelanceFilteredWide$co > 1] = 'adaptive both'

meanTissuePnPsPrevelanceFilteredWide%>%
  ggplot(aes(x = dj,
             y = co,
             col = selection))+
  geom_point(alpha = .5)

intersect(meanTissuePnPsPrevelanceFilteredWide$gene, meanTissuePnPsWide$gene)

genesInfo %>%
  filter(SNV_count >= 3) %>%
  group_by(sample, SNV_N_count, tissue) %>%  # Group by tissue
  mutate(pNpS_variants = ifelse(is.na(pNpS_variants), mean(pNpS_variants, na.rm = TRUE), pNpS_variants)) %>%
  pivot_wider(names_from = gene, id_cols = c(sample, tissue), values_from = pNpS_variants) %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  # Group by tissue before filtering
  select(sample, tissue, where(~ mean(. == 0, na.rm = TRUE) <= 0.75)) %>%
  pivot_longer(c(3:ncol(.)), values_to = 'pNpS_variants', names_to = 'gene')%>%
  view()


