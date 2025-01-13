library(tidyverse)
library(readr)
library(dplyr)
spades<-read.csv('pangenomicsAnalysis/spades/12_SYNTRACKER/output/summary_output/synteny_scores_per_region.csv')
megahit <-read.csv('pangenomicsAnalysis/megahit/12_SYNTRACKER/output/summary_output/synteny_scores_per_region.csv')

lowSyntenyRegionsSpadesAssembly<-spades %>%
  as.data.frame()%>%
  filter(Synteny_score <.955) %>%
  group_by(Region) %>%
  summarise(n()) %>%
  dplyr::select(Region)

lowSyntenyRegionsMegahitAssembly<-megahit %>%
  as.data.frame()%>%
  filter(Synteny_score <.955) %>%
  group_by(Region) %>%
  summarise(n()) %>%
  dplyr::select(Region)
sharedLowSyntenyRegions<-intersect(lowSyntenyRegionsSpadesAssembly$Region,
      lowSyntenyRegionsMegahitAssembly$Region)
'6 regions shared between assemnly methods with synteny scores below the "same strain" threshold'
startStop <- do.call(rbind, strsplit(sharedLowSyntenyRegions, "_"))
startStop<-startStop[,3:4]
colnames(startStop) <-c('positonStart', 'positionStop')
write_tsv(as.data.frame(startStop), file = 'syntentyAnalysis/sharedLowSyntenyRegions.tsv')