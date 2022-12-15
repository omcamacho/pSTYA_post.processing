
# Figure 5 #
############

source('D:/Paper submission/R Functions/Binomial Scores Function Peptidoform level.R')
source('D:/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('D:/Paper submission/R Functions/Function frequency of site.R')


library("conflicted")
suppressPackageStartupMessages(library("tidyverse"))
conflict_prefer("filter", "dplyr")

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


PXD000923 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD000923.csv')
PXD002222 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD002222.csv')
PXD002756 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD002756.csv')
PXD004705 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD004705.csv')
PXD004939 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD004939.csv')
PXD005241 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD005241.csv')
PXD012764 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD012764.csv')
PXD019291 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD019291.csv')


Bin_Adjusted_PXD000923 <- binAdjustPform(PXD000923)
Bin_Adjusted_PXD002222 <- binAdjustPform(PXD002222)
Bin_Adjusted_PXD002756 <- binAdjustPform(PXD002756)
Bin_Adjusted_PXD004705 <- binAdjustPform(PXD004705)
Bin_Adjusted_PXD004939 <- binAdjustPform(PXD004939)
Bin_Adjusted_PXD005241 <- binAdjustPform(PXD005241)
Bin_Adjusted_PXD012764 <- binAdjustPform(PXD012764)
Bin_Adjusted_PXD019291 <- binAdjustPform(PXD019291)


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


ggplot(AllRice_Pform,aes(x=Bin_Adjusted_Score,fill=Class))+
  geom_histogram(aes(y=..density..),
                 alpha=0.5,position='identity',binwidth=0.02)+
  facet_wrap(~dataset,nrow=2) + labs(x = "Scores", y = "Density")

