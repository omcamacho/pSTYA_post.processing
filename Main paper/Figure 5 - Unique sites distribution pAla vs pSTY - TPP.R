
# Figure 7

source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function Peptidoform level.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Function frequency of site.R')


library("conflicted")
suppressPackageStartupMessages(library("tidyverse"))
conflict_prefer("filter", "dplyr")
#> [conflicted] Will prefer dplyr::filter over any other package
library("betareg")
library("broom")
library("kableExtra")
library(plyr)
library(dplyr)
library(stringr)
library(useful)
library("data.table")
conflict_prefer("mutate", "dplyr")
library(reshape2)
library(epiDisplay)


PXD000923 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD000923.csv')
PXD002222 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD002222.csv')
PXD002756 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD002756.csv')
PXD004705 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD004705.csv')
PXD004939 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD004939.csv')
PXD005241 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD005241.csv')
PXD012764 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD012764.csv')
PXD019291 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD019291.csv')

PXD000923$PTM_score <- PXD000923$PTM.Score
PXD002222$PTM_score <- PXD002222$PTM.Score
PXD002756$PTM_score <- PXD002756$PTM.Score
PXD004705$PTM_score <- PXD004705$PTM.Score
PXD004939$PTM_score <- PXD004939$PTM.Score
PXD005241$PTM_score <- PXD005241$PTM.Score
PXD012764$PTM_score <- PXD012764$PTM.Score
PXD019291$PTM_score <- PXD019291$PTM.Score

PXD000923$PTM_positions <- PXD000923$PTM.positions
PXD002222$PTM_positions <- PXD002222$PTM.positions
PXD002756$PTM_positions <- PXD002756$PTM.positions
PXD004705$PTM_positions <- PXD004705$PTM.positions
PXD004939$PTM_positions <- PXD004939$PTM.positions
PXD005241$PTM_positions <- PXD005241$PTM.positions
PXD012764$PTM_positions <- PXD012764$PTM.positions
PXD019291$PTM_positions <- PXD019291$PTM.positions

Bin_Adjusted_PXD000923 <- binAdjustPform(PXD000923)
Bin_Adjusted_PXD002222 <- binAdjustPform(PXD002222)
Bin_Adjusted_PXD002756 <- binAdjustPform(PXD002756)
Bin_Adjusted_PXD004705 <- binAdjustPform(PXD004705)
Bin_Adjusted_PXD004939 <- binAdjustPform(PXD004939)
Bin_Adjusted_PXD005241 <- binAdjustPform(PXD005241)
Bin_Adjusted_PXD012764 <- binAdjustPform(PXD012764)
Bin_Adjusted_PXD019291 <- binAdjustPform(PXD019291)


Bin_Adjusted_PXD000923$Bin_Adjusted_Score <- Bin_Adjusted_PXD000923$NewScore3
Bin_Adjusted_PXD002222$Bin_Adjusted_Score <- Bin_Adjusted_PXD002222$NewScore3
Bin_Adjusted_PXD002756$Bin_Adjusted_Score <- Bin_Adjusted_PXD002756$NewScore3
Bin_Adjusted_PXD004705$Bin_Adjusted_Score <- Bin_Adjusted_PXD004705$NewScore3
Bin_Adjusted_PXD004939$Bin_Adjusted_Score <- Bin_Adjusted_PXD004939$NewScore3
Bin_Adjusted_PXD005241$Bin_Adjusted_Score <- Bin_Adjusted_PXD005241$NewScore3
Bin_Adjusted_PXD012764$Bin_Adjusted_Score <- Bin_Adjusted_PXD012764$NewScore3
Bin_Adjusted_PXD019291$Bin_Adjusted_Score <- Bin_Adjusted_PXD019291$NewScore3


Bin_Adjusted_PXD000923 <- FLR_Adj(Bin_Adjusted_PXD000923)
Bin_Adjusted_PXD002222 <- FLR_Adj(Bin_Adjusted_PXD002222)
Bin_Adjusted_PXD002756 <- FLR_Adj(Bin_Adjusted_PXD002756)
Bin_Adjusted_PXD004705 <- FLR_Adj(Bin_Adjusted_PXD004705)
Bin_Adjusted_PXD004939 <- FLR_Adj(Bin_Adjusted_PXD004939)
Bin_Adjusted_PXD005241 <- FLR_Adj(Bin_Adjusted_PXD005241)
Bin_Adjusted_PXD012764 <- FLR_Adj(Bin_Adjusted_PXD012764)
Bin_Adjusted_PXD019291 <- FLR_Adj(Bin_Adjusted_PXD019291)

Bin_Adjusted_PXD000923$dataset <- "PXD000923"
Bin_Adjusted_PXD002222$dataset <- "PXD002222"
Bin_Adjusted_PXD002756$dataset <- "PXD002756"
Bin_Adjusted_PXD004705$dataset <- "PXD004705"
Bin_Adjusted_PXD004939$dataset <- "PXD004939"
Bin_Adjusted_PXD005241$dataset <- "PXD005241"
Bin_Adjusted_PXD012764$dataset <- "PXD012764"
Bin_Adjusted_PXD019291$dataset <- "PXD019291"


AllRice_Pform<-dplyr::bind_rows(Bin_Adjusted_PXD000923, Bin_Adjusted_PXD002222, Bin_Adjusted_PXD002756, Bin_Adjusted_PXD004705,
                                      Bin_Adjusted_PXD004939, Bin_Adjusted_PXD005241, Bin_Adjusted_PXD012764, Bin_Adjusted_PXD019291)


# GRAPH #

AllRice_Pform$PepwithA<- as.character(gregexpr(pattern ='A\\[Phospho',AllRice_Pform$Peptide_mod))


AllRice_Pform$Class <-ifelse(AllRice_Pform$PepwithA=='-1',"Target","Decoy")


ggplot(AllRice_Pform,aes(x=NewScore3,fill=Class))+
  geom_histogram(aes(y=..density..),
                 alpha=0.5,position='identity',binwidth=0.02)+
  facet_wrap(~dataset,nrow=2) + labs(x = "Scores", y = "Density")

