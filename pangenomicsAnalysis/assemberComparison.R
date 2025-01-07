'Summary:
Looks like from a gene content POV, megahit did a better job at assembling.
The amount of genes is still in excess of the reference but there are less genes
overall relative to the megahit assembly (Fig 1), less partial genes (Fig 2)
and less complete genes (Fig 3)'



library(tidyverse)
library(ggplot2)
library(dtplyr)
library(dplyr)
library(readr)
library(ggpubr)
assemblyStats<-read_tsv(file = 'pangenomicsAnalysis/assemblerStats.tsv')
ggplot(data = assemblyStats,
       aes(x = assembly,
           y = totalGenes))+
  geom_boxplot()

ggplot(data = assemblyStats,
       aes(x = assembly,
           y = completeGenes))+
  geom_boxplot()
ggplot(data = assemblyStats,
       aes(y = partialGenes,
           x = assembly))+
  geom_boxplot()
'In general it seems like the megahit assembly is closer to the gene content of the reference'

totalGenesReference = 2051
partialGenesReference = 1
completeGenesReference = 2050
assemblyStats$totalSurplusGenes = assemblyStats$totalGenes - totalGenesReference
assemblyStats$surplusPartialGenes = assemblyStats$partialGenes - partialGenesReference
assemblyStats$surplusCompleteGenes = assemblyStats$completeGenes - completeGenesReference
'Statistically more excess genes in the spades assembly'
assemblyStats %>%
  filter(assembly != 'Reference') %>%
  ggplot(aes(x = assembly,
             y = abs(totalSurplusGenes)))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = 'fig 1')
'More partial gene callers, significant but the difference is not as bad.'
assemblyStats %>%
  filter(assembly != 'Reference') %>%
  ggplot(aes(x = assembly,
             y = abs(surplusPartialGenes)))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = 'fig 2')
assemblyStats %>%
  filter(assembly != 'Reference') %>%
  ggplot(aes(x = assembly,
             y = abs(surplusCompleteGenes)))+
  geom_boxplot()+
  stat_compare_means()+
  labs(title = 'Fig 3')
'Way more complete genes as well'
