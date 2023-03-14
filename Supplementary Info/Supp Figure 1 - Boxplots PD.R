

# Supp Figure 1 - BOXPLOTS PD #


library("conflicted")
suppressPackageStartupMessages(library("tidyverse"))
conflict_prefer("filter", "dplyr")
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

source('D:/Paper submission/R Functions/Binomial Scores Function PSM level.R')
source('D:/Paper submission/R Functions/Function frequency of site.R')


Ara_PXD008355 <- read.csv(file = 'D:/Paper submission/Data PD pipeline/PD_PXD008355.csv')
Hum_PXD000612 <- read.csv(file = 'D:/Paper submission/Data PD pipeline/PD_PXD000612.csv')

Ara_PXD008355$dataset <- "Ara_PXD008355"
Hum_PXD000612$dataset <- "Hum_PXD000612"

Ara_PXD008355 <- rename(Ara_PXD008355, PTM.positions = PTM_positions)
Ara_PXD008355 <- rename(Ara_PXD008355, PTM.Score = PTM_score)
Hum_PXD000612 <- rename(Hum_PXD000612, PTM.positions = PTM_positions)
Hum_PXD000612 <- rename(Hum_PXD000612, PTM.Score = PTM_score)

L_Ara_PXD008355 <- binAdjustPSM(Ara_PXD008355)
L_Hum_PXD000612 <- binAdjustPSM(Hum_PXD000612)

All_data_PD <- dplyr::bind_rows(L_Ara_PXD008355, L_Hum_PXD000612)

All_data_PD$Amino <- str_sub(substr(All_data_PD$Peptide,1,All_data_PD$PTM.positions),-1)

GraphData_PD <- subset(All_data_PD,select=c(dataset, Amino, PROTEIN_LOC, PTM_final_prob))

GraphData_count_PD <- GraphData_PD %>% dplyr::count(dataset, PROTEIN_LOC, Amino)

ggplot(GraphData_count_PD, aes(Amino, n)) +
  geom_boxplot() + ggtitle("All Phospho-sites") +
  xlab("Amino Acid") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)+ theme(axis.text = element_text(size = 12))

ggplot(GraphData_count_PD, aes(Amino, n)) +
  geom_boxplot() + scale_y_continuous(limits=c(0,50)) + ggtitle("Phospho-sites observed up to 50 times") +
  xlab("Amino Acid") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)+ theme(axis.text = element_text(size = 12))


GraphData_95_PD <- GraphData_PD[GraphData_PD$PTM_final_prob>0.95,]

GraphData_95count_PD <- GraphData_95_PD %>% dplyr::count(dataset, PROTEIN_LOC, Amino)


ggplot(GraphData_95count_PD, aes(Amino, n)) +
  geom_boxplot() + scale_y_continuous(limits=c(0,50)) + ggtitle("Phospho-sites observed up to 50 times and Probability>0.95") +
  xlab("Amino Acid") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)+ theme(axis.text = element_text(size = 12))


