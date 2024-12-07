setwd('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData')
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(FactoMineR)
library(ggExtra)
library(vegan)

'Read in data and filter, remove samples with a min coverage less than 45, this removes the problematic samples and
results in serpation via PCA. '
genesInfo<-read_csv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData/bothCohortsGenesInfo.csv')
genomeInfo<- read_csv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData/bothCohortsGenomeInfo.csv')

passedFiltered<-genomeInfo %>%
  filter(coverage_median > 45) %>%
  .$sample
genesFiltered<- filteredGenesInfo<-genesInfo %>%
  filter(sample %in% passedFiltered) %>%
  filter(breadth_minCov >.9)

metadata<-filteredGenesInfo %>%
  dplyr::select(sample, cohort) %>%
  distinct()

metadata$tissue <- 'colon'
metadata$tissue[grepl(metadata$sample, pattern = 'dj', ignore.case = TRUE)]<- 'dj'

genesFiltered$tissue <- 'colon'
genesFiltered$tissue[grepl(genesFiltered$sample, pattern = 'dj', ignore.case = TRUE)]<- 'dj'

'Only consider genes which are present in at least 50% of samples with 3 or more snvs'
nDj<-genesFiltered %>%
  filter(tissue == 'dj')%>%
  select(sample)%>%
  distinct()%>%
  nrow()
cutoffDj<-nDj*.7
passedDj<-genesFiltered %>%
  filter(SNV_count > 2,
         tissue == 'dj') %>%
  group_by(gene) %>%
  summarise(n()) %>%
  filter(`n()`>= cutoffDj) %>%
  .$gene

view(genesFiltered)

nColon<-genesFiltered %>%
  filter(tissue == 'colon')%>%
  select(sample)%>%
  distinct()%>%
  nrow()
cutoffColon<-nColon*.7

passedColon<-genesFiltered %>%
  filter(SNV_count > 2,
         tissue == 'colon') %>%
  group_by(gene) %>%
  summarise(n()) %>%
  filter(`n()`>= cutoffColon) %>%
  .$gene
passedBoth<-intersect(passedColon,
                      passedDj)

genesPassedPrevelanceFilter<-genesFiltered%>%
  filter(gene %in% passedBoth)

kofams<-read_tsv('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/b. breve/b. breve wt transplants/firstCohortRerunAndSecondCohortSnvData/kofams.tsv')
functions<-kofams %>%
  group_by(gene) %>%
  filter(e_value == min(e_value))

genesmeanNA<-genesPassedPrevelanceFilter %>%
  filter(SNV_count > 2) %>%
  group_by(sample) %>%
  mutate(pNpS_variants = ifelse(is.na(pNpS_variants) & SNV_count > 2, mean(pNpS_variants, na.rm = TRUE), pNpS_variants)) %>%
  ungroup()


genesmeanNA$selection <- 'none'
genesmeanNA$selection[genesmeanNA$pNpS_variants > 1 ] <- 'weak adaptive'
genesmeanNA$selection[genesmeanNA$pNpS_variants < 1 ] <- 'weak purifying'
genesmeanNA$selection[genesmeanNA$pNpS_variants > 1.25 ] <- 'adaptive'
genesmeanNA$selection[genesmeanNA$pNpS_variants < .8 ] <- 'purifying'
genesmeanNA$selection[genesmeanNA$pNpS_variants > 2.5] <- 'Strong adaptive'
genesmeanNA$selection[genesmeanNA$pNpS_variants < .364] <- 'Strong purifying'


genesmeanNA %>%
  filter(selection == 'none')
genesmeanNAMca<-genesmeanNA %>%
  filter()%>%
  pivot_wider(id_cols = sample,
              names_from = gene,
              values_from = selection)

genesmeanNAMca[is.na(genesmeanNAMca)]<- 'no evidence'


mcaOut<-genesmeanNAMca %>%
  column_to_rownames('sample') %>%
  mutate_if(is.character, as.factor) %>%
  MCA(graph = FALSE)
mcaOut$eig

plotData<- mcaOut$ind[1] %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  merge(metadata, by = 'sample')


p<-ggplot(data = plotData,
          aes(x = coord.Dim.1,
              y = coord.Dim.2,
              col = tissue,
              label = sample))+
  geom_point()
ggMarginal(p,groupFill = TRUE,type = 'boxplot')



ggplot(data = plotData,
       aes(x = coord.Dim.3,
           y = coord.Dim.4,
           col = as.factor(tissue),
           label = sample))+
  geom_point()
'Probably the best ordination with the mean imputed data
'
ggplot(data = plotData,
       aes(x = coord.Dim.2,
           y = coord.Dim.3,
           col = as.factor(tissue),
           label = sample))+
  geom_point()+stat_ellipse(level = .9)+stat_stars()
ggplot(data = plotData,
       aes(x = coord.Dim.2,
           y = coord.Dim.3,
           col = as.factor(tissue),
           label = sample))+
  geom_point()+stat_ellipse(level = .95)


genesMeanPCa<-genesmeanNA %>%
  filter()%>%
  pivot_wider(id_cols = sample,
              names_from = gene,
              values_from = pNpS_variants)

genesMeanPCa[is.na(genesMeanPCa)]<- -.0000001

genesMeanPCaScores<-genesMeanPCa %>%
  column_to_rownames('sample')%>%
  prcomp(center=TRUE,
         scale =TRUE)

summary(genesMeanPCaScores)

genesMeanPCaScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  merge(metadata,
        by = 'sample') %>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = tissue))+
  geom_point()
"Best, pca I think"
genesMeanPCaScores$x %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  merge(metadata,
        by = 'sample') %>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = tissue))+
  geom_point()


metafull <- read.csv('snvMeta.csv')
metafull<-metafull %>%
  select(-1)%>%
  select(sample,cage,Mouse,Cohort)
metadata<-metadata %>%
  rename('run' = cohort)
adonisData <- genesMeanPCa %>%
  column_to_rownames('sample')
adonis2(adonisData~tissue+cage+Mouse+Cohort+run, data = merge(metadata,
                                             metafull, by = 'sample'), method = 'bray')
adonis2(adonisData~tissue+cage+Mouse, data = merge(metadata,
                                        metafull, by = 'sample'), method = 'man')
adonis2(adonisData~tissue+cage+Mouse+run+Cohort, data = merge(metadata,
                                             metafull, by = 'sample'), method = 'jac')
'so far the best, not significant. Best with all terms. '
adonis2(adonisData~tissue+cage+Mouse+run+Cohort, data = merge(metadata,
                                             metafull, by = 'sample'), method = 'canberra')
adonis2(adonisData~tissue+cage, data = merge(metadata,
                                                              metafull, by = 'sample'), method = 'gower')
adonis2(adonisData~tissue+cage, data = merge(metadata,
                                                              metafull, by = 'sample'), method = 'altGower')
