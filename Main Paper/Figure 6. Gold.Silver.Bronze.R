
# Figure 6#
###########


library(plyr)
library(dplyr)
library(stringr)
library(useful)
library("data.table")
library("conflicted")
library(reshape2)
conflict_prefer("mutate", "dplyr")
suppressPackageStartupMessages(library("tidyverse"))
conflict_prefer("filter", "dplyr")
library(ggplot2)
library(dplyr)
library(epiDisplay)
library(gmodels)
library(webr)
library(ggsunburst)

source('D:/Paper submission/R Functions/Binomial Scores Function Peptidoform level.R')
source('D:/Paper submission/R Functions/GSB_Function.R')
source('D:/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('D:/Paper submission/R Functions/Function frequency of site.R')


# Peptidoform level using maximum #

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

Adjusted_PXD000923 <- FLR_Adj(Bin_Adjusted_PXD000923)
Adjusted_PXD002222 <- FLR_Adj(Bin_Adjusted_PXD002222)
Adjusted_PXD002756 <- FLR_Adj(Bin_Adjusted_PXD002756)
Adjusted_PXD004705 <- FLR_Adj(Bin_Adjusted_PXD004705)
Adjusted_PXD004939 <- FLR_Adj(Bin_Adjusted_PXD004939)
Adjusted_PXD005241 <- FLR_Adj(Bin_Adjusted_PXD005241)
Adjusted_PXD012764 <- FLR_Adj(Bin_Adjusted_PXD012764)
Adjusted_PXD019291 <- FLR_Adj(Bin_Adjusted_PXD019291)

Adjusted_PXD000923$dataset <- "PXD000923"
Adjusted_PXD002222$dataset <- "PXD002222"
Adjusted_PXD002756$dataset <- "PXD002756"
Adjusted_PXD004705$dataset <- "PXD004705"
Adjusted_PXD004939$dataset <- "PXD004939"
Adjusted_PXD005241$dataset <- "PXD005241"
Adjusted_PXD012764$dataset <- "PXD012764"
Adjusted_PXD019291$dataset <- "PXD019291"


Adjusted_PXD000923_01 <- Adjusted_PXD000923[1:max(which(Adjusted_PXD000923$FLR_Adj_Score<=0.01)),]
Adjusted_PXD002222_01 <- Adjusted_PXD002222[1:max(which(Adjusted_PXD002222$FLR_Adj_Score<=0.01)),]
Adjusted_PXD002756_01 <- Adjusted_PXD002756[1:max(which(Adjusted_PXD002756$FLR_Adj_Score<=0.01)),]
Adjusted_PXD004705_01 <- Adjusted_PXD004705[1:max(which(Adjusted_PXD004705$FLR_Adj_Score<=0.01)),]
Adjusted_PXD004939_01 <- Adjusted_PXD004939[1:max(which(Adjusted_PXD004939$FLR_Adj_Score<=0.01)),]
Adjusted_PXD005241_01 <- Adjusted_PXD005241[1:max(which(Adjusted_PXD005241$FLR_Adj_Score<=0.01)),]
Adjusted_PXD012764_01 <- Adjusted_PXD012764[1:max(which(Adjusted_PXD012764$FLR_Adj_Score<=0.01)),]
Adjusted_PXD019291_01 <- Adjusted_PXD019291[1:max(which(Adjusted_PXD019291$FLR_Adj_Score<=0.01)),]

Adjusted_PXD000923_05 <- Adjusted_PXD000923[1:max(which(Adjusted_PXD000923$FLR_Adj_Score<=0.05)),]
Adjusted_PXD002222_05 <- Adjusted_PXD002222[1:max(which(Adjusted_PXD002222$FLR_Adj_Score<=0.05)),]
Adjusted_PXD002756_05 <- Adjusted_PXD002756[1:max(which(Adjusted_PXD002756$FLR_Adj_Score<=0.05)),]
Adjusted_PXD004705_05 <- Adjusted_PXD004705[1:max(which(Adjusted_PXD004705$FLR_Adj_Score<=0.05)),]
Adjusted_PXD004939_05 <- Adjusted_PXD004939[1:max(which(Adjusted_PXD004939$FLR_Adj_Score<=0.05)),]
Adjusted_PXD005241_05 <- Adjusted_PXD005241[1:max(which(Adjusted_PXD005241$FLR_Adj_Score<=0.05)),]
Adjusted_PXD012764_05 <- Adjusted_PXD012764[1:max(which(Adjusted_PXD012764$FLR_Adj_Score<=0.05)),]
Adjusted_PXD019291_05 <- Adjusted_PXD019291[1:max(which(Adjusted_PXD019291$FLR_Adj_Score<=0.05)),]


AllRice_pform_Max<-dplyr::bind_rows(Adjusted_PXD000923, Adjusted_PXD002222, Adjusted_PXD002756, Adjusted_PXD004705,
                                      Adjusted_PXD004939, Adjusted_PXD005241, Adjusted_PXD012764, Adjusted_PXD019291)

AllRice_pform_Max_01<-dplyr::bind_rows(Adjusted_PXD000923_01, Adjusted_PXD002222_01, Adjusted_PXD002756_01, Adjusted_PXD004705_01,
                                    Adjusted_PXD004939_01, Adjusted_PXD005241_01, Adjusted_PXD012764_01, Adjusted_PXD019291_01)

AllRice_pform_Max_05<-dplyr::bind_rows(Adjusted_PXD000923_05, Adjusted_PXD002222_05, Adjusted_PXD002756_05, Adjusted_PXD004705_05,
                                       Adjusted_PXD004939_05, Adjusted_PXD005241_05, Adjusted_PXD012764_05, Adjusted_PXD019291_05)


AllRice_pform_Max$New_FLR_PEP <-AllRice_pform_Max$FLR_Adj_Score
AllRice_pform_Max_01$New_FLR_PEP <-AllRice_pform_Max_01$FLR_Adj_Score
AllRice_pform_Max_05$New_FLR_PEP <-AllRice_pform_Max_05$FLR_Adj_Score

AllRice_pform_Max_Final <- GSB_Function(AllRice_pform_Max_01,AllRice_pform_Max_05)

str(AllRice_pform_Max_Final)

AllRice_pform_Max_Final$PROTEIN_POS<-NULL
AllRice_pform_Max_Final$PRO_pos_list<-NULL
AllRice_pform_Max_Final$PTM_length<-NULL
AllRice_pform_Max_Final$PTM_beg2<-NULL
AllRice_pform_Max_Final$PTM_end2<-NULL
AllRice_pform_Max_Final$PTM_End<-NULL
AllRice_pform_Max_Final$PTM_Beginning<-NULL

RiceY <- AllRice_pform_Max_Final[(AllRice_pform_Max_Final$Amino=="Y") & (AllRice_pform_Max_Final$cat != "Bronze"),]


CrossTable(AllRice_pform_Max_Final$cat, AllRice_pform_Max_Final$Amino)

df2 <- data.frame(Amino_Acid=rep(c("S", "T", "Y", "A"), each=3),
                  Level=rep(c("Gold:2114", "Silver:4254", "Bronze:5460"),4),
                  Unique_sites=c(1887, 3760, 4512, 213, 439, 603, 9, 22, 49, 5, 33, 296))
head(df2)

ggplot(data=df2, aes(x=Level, y=Unique_sites, fill=Amino_Acid)) +
  geom_bar(stat="identity")

p <- df2 %>%
  dplyr::arrange(Unique_sites) %>%
  mutate(Level = factor(Level, levels=c("Gold:2114", "Silver:4254", "Bronze:5460"))) %>%
  ggplot(aes(x=Level, y=Unique_sites, fill=Amino_Acid)) +
  geom_bar(stat="identity") + theme(axis.text = element_text(size = 12))  +
  xlab("") 

p + geom_text(aes(label = Unique_sites), position = position_stack(vjust = 0.9)) 

