

# Supp Figure 1 - BOXPLOTS PD #


library("conflicted")
suppressPackageStartupMessages(library("tidyverse"))
conflict_prefer("filter", "dplyr")
#> [conflicted] Will prefer dplyr::filter over any other package
library("broom")
library("kableExtra")
library(plyr)
library(dplyr)
library(stringr)
library(useful)
library("data.table")
conflict_prefer("mutate", "dplyr")
conflict_prefer("rename", "dplyr")
library(epiDisplay)
library(ggplot2)
library(janitor)

source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level Arabidopsis.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Function frequency of site.R')


Ara_PXD008355 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD008355.csv')
Hum_PXD000612 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD000612.csv')

Ara_PXD008355$dataset <- "Ara_PXD008355"
Hum_PXD000612$dataset <- "Hum_PXD000612"

L_Ara_PXD008355 <- binAdjustAra(Ara_PXD008355)
L_Hum_PXD000612 <- binAdjustPSM(Hum_PXD000612)

All_data_PD <- dplyr::bind_rows(L_Ara_PXD008355, L_Hum_PXD000612)

GraphData_PD <- subset(All_data_PD,select=c(dataset, Amino, PROTEIN_LOC, PTM_final_prob))

GraphData_PD$Match <- ifelse(GraphData_PD$Amino=="A","pAla","pSTY")

GraphData_count_PD <- GraphData_PD %>% dplyr::count(dataset, PROTEIN_LOC, Match)

ggplot(GraphData_count_PD, aes(Match, n)) +
  geom_boxplot() + ggtitle("All Phospho-sites") +
  xlab("Type of Match") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)

ggplot(GraphData_count_PD, aes(Match, n)) +
  geom_boxplot() + scale_y_continuous(limits=c(0,50)) + ggtitle("Phospho-sites observed up to 50 times") +
  xlab("Type of Match") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)


GraphData_95_PD <- GraphData_PD[GraphData_PD$PTM_final_prob>0.95,]

GraphData_95count_PD <- GraphData_95_PD %>% dplyr::count(dataset, PROTEIN_LOC, Match)


ggplot(GraphData_95count_PD, aes(Match, n)) +
  geom_boxplot() + scale_y_continuous(limits=c(0,50)) + ggtitle("Phospho-sites observed up to 50 times and Probability>0.95") +
  xlab("Type of Match") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)


