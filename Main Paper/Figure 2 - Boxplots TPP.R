
# Figure 2 #
# BOXPLOTS #


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

source('D:/Paper submission/R Functions/Binomial Scores Function PSM level.R')
source('D:/Paper submission/R Functions/Function frequency of site.R')

Ara_PXD008355 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD008355.csv')
Hum_PXD000612 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD000612.csv')
Ric_PXD000923 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD000923.csv')
Ric_PXD002222 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD002222.csv')
Ric_PXD002756 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD002756.csv')
Ric_PXD004705 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD004705.csv')
Ric_PXD004939 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD004939.csv')
Ric_PXD005241 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD005241.csv')
Ric_PXD012764 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD012764.csv')
Ric_PXD019291 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD019291.csv')

Ara_PXD008355$dataset <- "Ara_PXD008355"
Hum_PXD000612$dataset <- "Hum_PXD000612"
Ric_PXD000923$dataset <- "Ric_PXD000923"
Ric_PXD002222$dataset <- "Ric_PXD002222"
Ric_PXD002756$dataset <- "Ric_PXD002756"
Ric_PXD004705$dataset <- "Ric_PXD004705"
Ric_PXD004939$dataset <- "Ric_PXD004939"
Ric_PXD005241$dataset <- "Ric_PXD005241"
Ric_PXD012764$dataset <- "Ric_PXD012764"
Ric_PXD019291$dataset <- "Ric_PXD019291"

Ara_PXD008355 <- rename(Ara_PXD008355, PTM.positions = PTM_positions)
Hum_PXD000612 <- rename(Hum_PXD000612, PTM.positions = PTM_positions)
Ara_PXD008355 <- rename(Ara_PXD008355, PTM.Score = PTM_score)
Hum_PXD000612 <- rename(Hum_PXD000612, PTM.Score = PTM_score)

L_Ara_PXD008355 <- binAdjustPSM(Ara_PXD008355)
L_Hum_PXD000612 <- binAdjustPSM(Hum_PXD000612)
L_Ric_PXD000923 <- binAdjustPSM(Ric_PXD000923)
L_Ric_PXD002222 <- binAdjustPSM(Ric_PXD002222)
L_Ric_PXD002756 <- binAdjustPSM(Ric_PXD002756)
L_Ric_PXD004705 <- binAdjustPSM(Ric_PXD004705)
L_Ric_PXD004939 <- binAdjustPSM(Ric_PXD004939)
L_Ric_PXD005241 <- binAdjustPSM(Ric_PXD005241)
L_Ric_PXD012764 <- binAdjustPSM(Ric_PXD012764)
L_Ric_PXD019291 <- binAdjustPSM(Ric_PXD019291)

All_data <- dplyr::bind_rows(L_Ara_PXD008355, L_Hum_PXD000612, L_Ric_PXD000923, L_Ric_PXD002222,
                             L_Ric_PXD002756, L_Ric_PXD004705, L_Ric_PXD004939, L_Ric_PXD005241,
                             L_Ric_PXD012764, L_Ric_PXD019291)

All_data$Amino <- str_sub(substr(All_data$Peptide,1,All_data$PTM.positions),-1)

GraphData <- subset(All_data,select=c(dataset, Amino, PROTEIN_LOC, PTM_final_prob))

GraphData$Match <- ifelse(GraphData$Amino=="A","pAla","pSTY")

GraphData_count <- GraphData %>% dplyr::count(dataset, PROTEIN_LOC, Match)

ggplot(GraphData_count, aes(Match, n)) +
  geom_boxplot() + ggtitle("All Phosphosites") +
  xlab("Type of Match") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)

ggplot(GraphData_count, aes(Match, n)) +
  geom_boxplot() + scale_y_continuous(limits=c(0,50)) + ggtitle("Phosphosites observed up to 50 times") +
  xlab("Type of Match") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)


GraphData_95 <- GraphData[GraphData$PTM_final_prob>0.95,]

GraphData_95count <- GraphData_95 %>% dplyr::count(dataset, PROTEIN_LOC, Match)


ggplot(GraphData_95count, aes(Match, n)) +
  geom_boxplot() + scale_y_continuous(limits=c(0,50)) + ggtitle("Phosphosites observed up to 50 times and Probability>0.95") +
  xlab("Type of Match") + ylab("Number of PSMs") + facet_wrap(~dataset, ncol = 5)


