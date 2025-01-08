library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
apssStrainCutoff = .955

data = read_csv('pangenomicsAnalysis/megahit/12_SYNTRACKER/output/summary_output/avg_synteny_scores_40_regions.csv')
'No different strains at 40 regions'
data %>%
  filter(APSS < apssStrainCutoff)

data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .01)

data$sample1Tissue <- 'co'
data$sample1Tissue[grepl(data$Sample1, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$sample2Tissue <- 'co'
data$sample2Tissue[grepl(data$Sample2, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$tissueComparison <- 'differentTissue'
data$tissueComparison[data$sample1Tissue =='dj' & data$sample2Tissue == 'dj'] <- 'dj'
data$tissueComparison[data$sample1Tissue =='co' & data$sample2Tissue == 'co'] <- 'co'

data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))+
  labs(title = 'APSS 40 regions')

data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .005)+
  facet_wrap(~tissueComparison,
             ncol =1)

data %>%
  group_by(tissueComparison) %>%
  summarise(n())

data = read_csv('pangenomicsAnalysis/megahit/12_SYNTRACKER/output/summary_output/avg_synteny_scores_60_regions.csv')
'No different strains at 60 regions'
data %>%
  filter(APSS < apssStrainCutoff)

data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .01)

data$sample1Tissue <- 'co'
data$sample1Tissue[grepl(data$Sample1, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$sample2Tissue <- 'co'
data$sample2Tissue[grepl(data$Sample2, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$tissueComparison <- 'differentTissue'
data$tissueComparison[data$sample1Tissue =='dj' & data$sample2Tissue == 'dj'] <- 'dj'
data$tissueComparison[data$sample1Tissue =='co' & data$sample2Tissue == 'co'] <- 'co'

data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))


data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))+
  labs(title = 'APSS 60 regions')
data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .005)+
  facet_wrap(~tissueComparison,
             ncol =1)

data = read_csv('pangenomicsAnalysis/megahit/12_SYNTRACKER/output/summary_output/avg_synteny_scores_80_regions.csv')
'No different strains at 60 regions'
data %>%
  filter(APSS < apssStrainCutoff)

data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .01)

data$sample1Tissue <- 'co'
data$sample1Tissue[grepl(data$Sample1, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$sample2Tissue <- 'co'
data$sample2Tissue[grepl(data$Sample2, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$tissueComparison <- 'differentTissue'
data$tissueComparison[data$sample1Tissue =='dj' & data$sample2Tissue == 'dj'] <- 'dj'
data$tissueComparison[data$sample1Tissue =='co' & data$sample2Tissue == 'co'] <- 'co'

data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  labs(title = 'APSS 80 regions')+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
'No longer significant differences'
data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .005)+
  facet_wrap(~tissueComparison,
             ncol =1)

data = read_csv('pangenomicsAnalysis/megahit/12_SYNTRACKER/output/summary_output/avg_synteny_scores_100_regions.csv')

data %>%
  filter(APSS < apssStrainCutoff)

data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .01)

data$sample1Tissue <- 'co'
data$sample1Tissue[grepl(data$Sample1, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$sample2Tissue <- 'co'
data$sample2Tissue[grepl(data$Sample2, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$tissueComparison <- 'differentTissue'
data$tissueComparison[data$sample1Tissue =='dj' & data$sample2Tissue == 'dj'] <- 'dj'
data$tissueComparison[data$sample1Tissue =='co' & data$sample2Tissue == 'co'] <- 'co'

data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  labs(title = 'APSS 100 Regions')+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))

data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .005)+
  facet_wrap(~tissueComparison,
             ncol =1)

data = read_csv('pangenomicsAnalysis/megahit/12_SYNTRACKER/output/summary_output/avg_synteny_scores_200_regions.csv')
'No different strains at 60 regions'
data %>%
  filter(APSS < apssStrainCutoff)

data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .01)

data$sample1Tissue <- 'co'
data$sample1Tissue[grepl(data$Sample1, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$sample2Tissue <- 'co'
data$sample2Tissue[grepl(data$Sample2, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$tissueComparison <- 'differentTissue'
data$tissueComparison[data$sample1Tissue =='dj' & data$sample2Tissue == 'dj'] <- 'dj'
data$tissueComparison[data$sample1Tissue =='co' & data$sample2Tissue == 'co'] <- 'co'

data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))+
  labs(title = 'APSS 200 regions')+
  geom_jitter()

data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .005)+
  facet_wrap(~tissueComparison,
             ncol =1)
data = read_csv('pangenomicsAnalysis/megahit/12_SYNTRACKER/output/summary_output/avg_synteny_scores_all_regions.csv')

data %>%
  filter(APSS < apssStrainCutoff)

data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .01)

data$sample1Tissue <- 'co'
data$sample1Tissue[grepl(data$Sample1, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$sample2Tissue <- 'co'
data$sample2Tissue[grepl(data$Sample2, pattern = 'dj', ignore.case = TRUE)]<-'dj'
data$tissueComparison <- 'differentTissue'
data$tissueComparison[data$sample1Tissue =='dj' & data$sample2Tissue == 'dj'] <- 'dj'
data$tissueComparison[data$sample1Tissue =='co' & data$sample2Tissue == 'co'] <- 'co'

data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))+
  labs(title = 'APSS all regions')
'DJ has statistically significantly lower APSS than different tissues comparisons
but still not different strains at 60 regions'
data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_boxplot()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = tissueComparison,
             y = APSS))+
  geom_violin()+
  geom_jitter()+
  stat_compare_means(comparisons = list(c('co','differentTissue'),
                                        c('co', 'dj'),
                                        c('differentTissue', 'dj')))
data %>%
  ggplot(aes(x = APSS))+
  geom_histogram(binwidth = .005)+
  facet_wrap(~tissueComparison,
             ncol =1)
data = read_csv('pangenomicsAnalysis/megahit/12_SYNTRACKER/output/summary_output/synteny_scores_per_region.csv')
data %>%
  filter(Synteny_score < apssStrainCutoff)%>%
  group_by(Region) %>%
  summarise(n())
'This just looks like a region that was hard to assemble. All samples are very different from each other. In general it looks like this is just do to some weirdness in the samples
'
'This regions a lot of rearrangements in both assemblies'
data%>%
  filter(Region == 'NC_020517.1_1015000_1016000',
         Synteny_score < apssStrainCutoff)%>%
  view()

data%>%
  filter(Region == 'NC_020517.1_1180000_1181000',
         Synteny_score < apssStrainCutoff)%>%
  view()

data%>%
  filter(Region == 'NC_020517.1_1850000_1851000',
         Synteny_score < apssStrainCutoff)%>%
  view()