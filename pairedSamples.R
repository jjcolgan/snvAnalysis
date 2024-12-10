setwd('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData')
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(FactoMineR)
library(ggExtra)

'Read in data and filter, remove samples with a min coverage less than 45, this removes the problematic samples, identified
 via ord plots'
genesInfo<-read_csv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData/bothCohortsGenesInfo.csv')
genomeInfo<- read_csv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData/bothCohortsGenomeInfo.csv')

'Going to try with a less strict coverage regarding coverage'
passedFiltered<-genomeInfo %>%
  filter(coverage_median > 45) %>%
  .$sample

metadata <- read.csv('snvMeta.csv')
metadata<-metadata %>%
  select(-1)

pairedSamples<-metadata %>%
  group_by(Mouse)%>%
  filter(sample %in% passedFiltered) %>%
  summarise('nPerMouse' = n())%>%
  filter(nPerMouse > 1)%>%
  .$Mouse
pairedSamples

pairedGenes <- genesInfo %>%
  merge(metadata, by = 'sample')%>%
  filter(Mouse %in% pairedSamples)

'Same filtering scheme from previous analysis'
nDj<-pairedGenes %>%
  filter(tissue == 'dj')%>%
  select(sample)%>%
  distinct()%>%
  nrow()
cutoffDj<-nDj*.75
passedDj<-pairedGenes %>%
  filter(SNV_count > 2,
         tissue == 'dj') %>%
  group_by(gene) %>%
  summarise(n()) %>%
  filter(`n()`>= cutoffDj) %>%
  .$gene

nColon<-pairedGenes %>%
  filter(tissue == 'colon')%>%
  select(sample)%>%
  distinct()%>%
  nrow()
cutoffColon<-nColon*.75

passedColon<-pairedGenes %>%
  filter(SNV_count > 2,
         tissue == 'colon') %>%
  group_by(gene) %>%
  summarise(n()) %>%
  filter(`n()`>= cutoffColon) %>%
  .$gene
passedBoth<-intersect(passedColon,
                      passedDj)
'This yields an ass load more genes'
genesPassedPrevelanceFilter<-pairedGenes%>%
  filter(gene %in% passedBoth)

'Quick PCA / MCA'
kofams<-read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData/kofams.tsv')
functions<-kofams %>%
  group_by(gene) %>%
  filter(e_value == min(e_value))

genesmeanNA<-pairedGenes %>%
  filter(SNV_count > 2)

genesmeanNA<- genesmeanNA %>%
  group_by(sample) %>%
  mutate(pNpS_variants = ifelse(is.na(pNpS_variants) & SNV_count > 2, mean(pNpS_variants, na.rm = TRUE), pNpS_variants)) %>%
  ungroup()


genesmeanNA$selection <- 'none'
genesmeanNA$selection[genesmeanNA$pNpS_variants > 1 ] <- 'weak adaptive'
genesmeanNA$selection[genesmeanNA$pNpS_variants < 1 ] <- 'weak purifying'
#genesmeanNA$selection[genesmeanNA$pNpS_variants > 1.25 ] <- 'adaptive'
#genesmeanNA$selection[genesmeanNA$pNpS_variants < .8 ] <- 'purifying'
#genesmeanNA$selection[genesmeanNA$pNpS_variants > 2.5] <- 'Strong adaptive'
#genesmeanNA$selection[genesmeanNA$pNpS_variants < .364] <- 'Strong purifying'


genesmeanNA %>%
  filter(selection == 'none')
genesmeanNAMca<-genesmeanNA %>%
  filter()%>%
  pivot_wider(id_cols = sample,
              names_from = gene,
              values_from = selection)

genesmeanNAMca[is.na(genesmeanNAMca)]<- 'none'


mcaOut<-genesmeanNAMca %>%
  column_to_rownames('sample') %>%
  mutate_if(is.character, as.factor) %>%
  MCA(graph = FALSE)
mcaOut$eig

plotData<- mcaOut$ind[1] %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  merge(metadata, by = 'sample')

'This did not work, lots of variability, but I think I can actually resolve differences pretty well'
ggplot(data = plotData,
          aes(x = coord.Dim.1,
              y = coord.Dim.2,
              col = tissue,
              label = sample))+
  geom_point()+
  labs(title = 'MCA mean paired')

genesMeanPCa<-genesmeanNA %>%
  filter()%>%
  pivot_wider(id_cols = sample,
              names_from = gene,
              values_from = pNpS_variants)

genesMeanPCa[is.na(genesMeanPCa)]<- -1

genesMeanPCaScores<-genesMeanPCa %>%
  column_to_rownames('sample')%>%
  prcomp(center=TRUE,
         scale =TRUE)

genesMeanPCaScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  merge(metadata,
        by = 'sample') %>%
  merge(genomeInfo, by = 'sample') %>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = tissue,
             label = sample))+
  geom_point()+labs('PCA meanNA paired')

genesMeanPCaScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  merge(metadata,
        by = 'sample') %>%
  merge(genomeInfo, by = 'sample') %>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = log10(coverage),
             label = sample))+
  geom_text_repel()+
  geom_point()+labs('PCA meanNA paired')

'This is not working, need to figure it out'
output <- data.frame('gene' = character(),
                     'pvalue' = numeric())
for (g in 1:length(passedBoth)){
  geneBeingTested <- passedBoth[g]
  temp <- genesmeanNA %>%
    filter(gene==geneBeingTested) %>%
    filter(is.na(pNpS_variants)== FALSE)
  pairedTemp<-temp %>%
    group_by(Mouse)%>%
    select(Mouse)%>%
    summarise('nPerMouse' = n())%>%
    filter(nPerMouse > 1)%>%
    .$Mouse
  print(print(length(pairedTemp)))
  temp<-temp %>%
    filter(Mouse %in% pairedTemp)
  dj <- temp %>%
    filter(tissue == 'dj')%>%
    .$pNpS_variants
  colon <- temp %>%
    filter(tissue == 'colon')%>%
    .$pNpS_variants
  wilcoxOut<- wilcox.test(dj, colon, paired = TRUE, exact = FALSE)
  tempOut <- data.frame('gene' = geneBeingTested,
                        'pvalue' = wilcoxOut$p.value)
  output <- rbind(output, tempOut)
}


output$padj <- p.adjust(output$pvalue, method = 'BH')

output %>%
  filter(pvalue < .05)
