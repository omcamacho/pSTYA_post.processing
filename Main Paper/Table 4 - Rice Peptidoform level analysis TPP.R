

# Table 4


source('D:/Paper submission/R Functions/Binomial Scores Function PSM level.R')
source('D:/Paper submission/R Functions/FLR function Unadjusted.R')
source('D:/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('D:/Paper submission/R Functions/Function frequency of site.R')


# We want to assess if there are any differences between the Max and MM collapsing methods using all rice data sets #

#> [conflicted] Will prefer dplyr::filter over any other package
suppressPackageStartupMessages(library("tidyverse"))

library(dplyr)
library(stringr)
library(useful)
library(MASS)
library(reshape2)
library(epiDisplay)

# First we calculate all binomial adjusted data #
#################################################

PXD000923 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD000923.csv')
PXD002222 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD002222.csv')
PXD002756 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD002756.csv')
PXD004705 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD004705.csv')
PXD004939 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD004939.csv')
PXD005241 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD005241.csv')
PXD012764 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD012764.csv')
PXD019291 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD019291.csv')

PXD000923$dataset <- "PXD000923"
PXD002222$dataset <- "PXD002222"
PXD002756$dataset <- "PXD002756"
PXD004705$dataset <- "PXD004705"
PXD004939$dataset <- "PXD004939"
PXD005241$dataset <- "PXD005241"
PXD012764$dataset <- "PXD012764"
PXD019291$dataset <- "PXD019291"


PXD000923$Peptidoform <- paste0(PXD000923$Peptide_mod,"_",PXD000923$PTM.positions)
PXD002222$Peptidoform <- paste0(PXD002222$Peptide_mod,"_",PXD002222$PTM.positions)
PXD002756$Peptidoform <- paste0(PXD002756$Peptide_mod,"_",PXD002756$PTM.positions)
PXD004705$Peptidoform <- paste0(PXD004705$Peptide_mod,"_",PXD004705$PTM.positions)
PXD004939$Peptidoform <- paste0(PXD004939$Peptide_mod,"_",PXD004939$PTM.positions)
PXD005241$Peptidoform <- paste0(PXD005241$Peptide_mod,"_",PXD005241$PTM.positions)
PXD012764$Peptidoform <- paste0(PXD012764$Peptide_mod,"_",PXD012764$PTM.positions)
PXD019291$Peptidoform <- paste0(PXD019291$Peptide_mod,"_",PXD019291$PTM.positions)


PXD000923$Unadjusted_Score <- PXD000923$PTM_final_prob
PXD002222$Unadjusted_Score <- PXD002222$PTM_final_prob
PXD002756$Unadjusted_Score <- PXD002756$PTM_final_prob
PXD004705$Unadjusted_Score <- PXD004705$PTM_final_prob
PXD004939$Unadjusted_Score <- PXD004939$PTM_final_prob
PXD005241$Unadjusted_Score <- PXD005241$PTM_final_prob
PXD012764$Unadjusted_Score <- PXD012764$PTM_final_prob
PXD019291$Unadjusted_Score <- PXD019291$PTM_final_prob




# Unadjusted data collapsed using Max method #

PXD000923_Max <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD000923, max)
PXD002222_Max <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD002222, max)
PXD002756_Max <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD002756, max)
PXD004705_Max <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD004705, max)
PXD004939_Max <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD004939, max)
PXD005241_Max <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD005241, max)
PXD012764_Max <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD012764, max)
PXD019291_Max <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD019291, max)

FLR_PXD000923_Max <- FLR_NotAdj(PXD000923_Max)
FLR_PXD002222_Max <- FLR_NotAdj(PXD002222_Max)
FLR_PXD002756_Max <- FLR_NotAdj(PXD002756_Max)
FLR_PXD004705_Max <- FLR_NotAdj(PXD004705_Max)
FLR_PXD004939_Max <- FLR_NotAdj(PXD004939_Max)
FLR_PXD005241_Max <- FLR_NotAdj(PXD005241_Max)
FLR_PXD012764_Max <- FLR_NotAdj(PXD012764_Max)
FLR_PXD019291_Max <- FLR_NotAdj(PXD019291_Max)


FLR_AllRice_pform_Max <-dplyr::bind_rows(FLR_PXD000923_Max,FLR_PXD002222_Max, FLR_PXD002756_Max,
                                   FLR_PXD004705_Max,FLR_PXD004939_Max, FLR_PXD005241_Max,
                                   FLR_PXD012764_Max,FLR_PXD019291_Max)

FLR_PXD000923_Max_01 <- FLR_PXD000923_Max[1:max(which(FLR_PXD000923_Max$FLR_Unadjusted<=0.01)),]
FLR_PXD002222_Max_01 <- FLR_PXD002222_Max[1:max(which(FLR_PXD002222_Max$FLR_Unadjusted<=0.01)),]
FLR_PXD002756_Max_01 <- FLR_PXD002756_Max[1:max(which(FLR_PXD002756_Max$FLR_Unadjusted<=0.01)),]
FLR_PXD004705_Max_01 <- FLR_PXD004705_Max[1:max(which(FLR_PXD004705_Max$FLR_Unadjusted<=0.01)),]
FLR_PXD004939_Max_01 <- FLR_PXD004939_Max[1:max(which(FLR_PXD004939_Max$FLR_Unadjusted<=0.01)),]
FLR_PXD005241_Max_01 <- FLR_PXD005241_Max[1:max(which(FLR_PXD005241_Max$FLR_Unadjusted<=0.01)),]
FLR_PXD012764_Max_01 <- FLR_PXD012764_Max[1:max(which(FLR_PXD012764_Max$FLR_Unadjusted<=0.01)),]
FLR_PXD019291_Max_01 <- FLR_PXD019291_Max[1:max(which(FLR_PXD019291_Max$FLR_Unadjusted<=0.01)),]

FLR_PXD000923_Max_05 <- FLR_PXD000923_Max[1:max(which(FLR_PXD000923_Max$FLR_Unadjusted<=0.05)),]
FLR_PXD002222_Max_05 <- FLR_PXD002222_Max[1:max(which(FLR_PXD002222_Max$FLR_Unadjusted<=0.05)),]
FLR_PXD002756_Max_05 <- FLR_PXD002756_Max[1:max(which(FLR_PXD002756_Max$FLR_Unadjusted<=0.05)),]
FLR_PXD004705_Max_05 <- FLR_PXD004705_Max[1:max(which(FLR_PXD004705_Max$FLR_Unadjusted<=0.05)),]
FLR_PXD004939_Max_05 <- FLR_PXD004939_Max[1:max(which(FLR_PXD004939_Max$FLR_Unadjusted<=0.05)),]
FLR_PXD005241_Max_05 <- FLR_PXD005241_Max[1:max(which(FLR_PXD005241_Max$FLR_Unadjusted<=0.05)),]
FLR_PXD012764_Max_05 <- FLR_PXD012764_Max[1:max(which(FLR_PXD012764_Max$FLR_Unadjusted<=0.05)),]
FLR_PXD019291_Max_05 <- FLR_PXD019291_Max[1:max(which(FLR_PXD019291_Max$FLR_Unadjusted<=0.05)),]


FLR_AllRice_pform_Max_R1<-dplyr::bind_rows(FLR_PXD000923_Max_01, FLR_PXD002222_Max_01, FLR_PXD002756_Max_01, FLR_PXD004705_Max_01,
                                       FLR_PXD004939_Max_01, FLR_PXD005241_Max_01, FLR_PXD012764_Max_01, FLR_PXD019291_Max_01)

FLR_AllRice_pform_Max_R5<-dplyr::bind_rows(FLR_PXD000923_Max_05, FLR_PXD002222_Max_05, FLR_PXD002756_Max_05, FLR_PXD004705_Max_05,
                                           FLR_PXD004939_Max_05, FLR_PXD005241_Max_05, FLR_PXD012764_Max_05, FLR_PXD019291_Max_05)


Table_FLR_AllRice_pform_Max_R5 <- FLR_AllRice_pform_Max_R5 %>%
  dplyr::select(dataset, Rnumber) %>% 
  group_by(dataset) %>% 
  summarise(NumPform_Max = max(Rnumber))

Table_FLR_AllRice_pform_Max_R1 <- FLR_AllRice_pform_Max_R1 %>%
  dplyr::select(dataset, Rnumber) %>% 
  group_by(dataset) %>% 
  summarise(NumPform_Max = max(Rnumber))

# Unadjusted data collapsed using MM method #

PXD000923_USI <- aggregate(Unadjusted_Score ~ All_USI, data = PXD000923, mean)
PXD002222_USI <- aggregate(Unadjusted_Score ~ All_USI, data = PXD002222, mean)
PXD002756_USI <- aggregate(Unadjusted_Score ~ All_USI, data = PXD002756, mean)
PXD004705_USI <- aggregate(Unadjusted_Score ~ All_USI, data = PXD004705, mean)
PXD004939_USI <- aggregate(Unadjusted_Score ~ All_USI, data = PXD004939, mean)
PXD005241_USI <- aggregate(Unadjusted_Score ~ All_USI, data = PXD005241, mean)
PXD012764_USI <- aggregate(Unadjusted_Score ~ All_USI, data = PXD012764, mean)
PXD019291_USI <- aggregate(Unadjusted_Score ~ All_USI, data = PXD019291, mean)

PXD000923_USI2 <- subset(PXD000923,select=-c(Unadjusted_Score))
PXD002222_USI2 <- subset(PXD002222,select=-c(Unadjusted_Score))
PXD002756_USI2 <- subset(PXD002756,select=-c(Unadjusted_Score))
PXD004705_USI2 <- subset(PXD004705,select=-c(Unadjusted_Score))
PXD004939_USI2 <- subset(PXD004939,select=-c(Unadjusted_Score))
PXD005241_USI2 <- subset(PXD005241,select=-c(Unadjusted_Score))
PXD012764_USI2 <- subset(PXD012764,select=-c(Unadjusted_Score))
PXD019291_USI2 <- subset(PXD019291,select=-c(Unadjusted_Score))

PXD000923_USI3 <- merge(PXD000923_USI,PXD000923_USI2, by="All_USI")
PXD002222_USI3 <- merge(PXD002222_USI,PXD002222_USI2, by="All_USI")
PXD002756_USI3 <- merge(PXD002756_USI,PXD002756_USI2, by="All_USI")
PXD004705_USI3 <- merge(PXD004705_USI,PXD004705_USI2, by="All_USI")
PXD004939_USI3 <- merge(PXD004939_USI,PXD004939_USI2, by="All_USI")
PXD005241_USI3 <- merge(PXD005241_USI,PXD005241_USI2, by="All_USI")
PXD012764_USI3 <- merge(PXD012764_USI,PXD012764_USI2, by="All_USI")
PXD019291_USI3 <- merge(PXD019291_USI,PXD019291_USI2, by="All_USI")

PXD000923_MM <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD000923_USI3, max)
PXD002222_MM <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD002222_USI3, max)
PXD002756_MM <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD002756_USI3, max)
PXD004705_MM <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD004705_USI3, max)
PXD004939_MM <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD004939_USI3, max)
PXD005241_MM <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD005241_USI3, max)
PXD012764_MM <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD012764_USI3, max)
PXD019291_MM <- aggregate(Unadjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = PXD019291_USI3, max)

FLR_PXD000923_MM <- FLR_NotAdj(PXD000923_MM)
FLR_PXD002222_MM <- FLR_NotAdj(PXD002222_MM)
FLR_PXD002756_MM <- FLR_NotAdj(PXD002756_MM)
FLR_PXD004705_MM <- FLR_NotAdj(PXD004705_MM)
FLR_PXD004939_MM <- FLR_NotAdj(PXD004939_MM)
FLR_PXD005241_MM <- FLR_NotAdj(PXD005241_MM)
FLR_PXD012764_MM <- FLR_NotAdj(PXD012764_MM)
FLR_PXD019291_MM <- FLR_NotAdj(PXD019291_MM)

FLR_AllRice_pform_MM <-dplyr::bind_rows(FLR_PXD000923_MM,FLR_PXD002222_MM, FLR_PXD002756_MM,
                                         FLR_PXD004705_MM,FLR_PXD004939_MM, FLR_PXD005241_MM,
                                         FLR_PXD012764_MM,FLR_PXD019291_MM)


FLR_PXD000923_MM_01 <- FLR_PXD000923_MM[1:max(which(FLR_PXD000923_MM$FLR_Unadjusted<=0.01)),]
FLR_PXD002222_MM_01 <- FLR_PXD002222_MM[1:max(which(FLR_PXD002222_MM$FLR_Unadjusted<=0.01)),]
FLR_PXD002756_MM_01 <- FLR_PXD002756_MM[1:max(which(FLR_PXD002756_MM$FLR_Unadjusted<=0.01)),]
FLR_PXD004705_MM_01 <- FLR_PXD004705_MM[1:max(which(FLR_PXD004705_MM$FLR_Unadjusted<=0.01)),]
FLR_PXD004939_MM_01 <- FLR_PXD004939_MM[1:max(which(FLR_PXD004939_MM$FLR_Unadjusted<=0.01)),]
FLR_PXD005241_MM_01 <- FLR_PXD005241_MM[1:max(which(FLR_PXD005241_MM$FLR_Unadjusted<=0.01)),]
FLR_PXD012764_MM_01 <- FLR_PXD012764_MM[1:max(which(FLR_PXD012764_MM$FLR_Unadjusted<=0.01)),]
FLR_PXD019291_MM_01 <- FLR_PXD019291_MM[1:max(which(FLR_PXD019291_MM$FLR_Unadjusted<=0.01)),]

FLR_PXD000923_MM_05 <- FLR_PXD000923_MM[1:max(which(FLR_PXD000923_MM$FLR_Unadjusted<=0.05)),]
FLR_PXD002222_MM_05 <- FLR_PXD002222_MM[1:max(which(FLR_PXD002222_MM$FLR_Unadjusted<=0.05)),]
FLR_PXD002756_MM_05 <- FLR_PXD002756_MM[1:max(which(FLR_PXD002756_MM$FLR_Unadjusted<=0.05)),]
FLR_PXD004705_MM_05 <- FLR_PXD004705_MM[1:max(which(FLR_PXD004705_MM$FLR_Unadjusted<=0.05)),]
FLR_PXD004939_MM_05 <- FLR_PXD004939_MM[1:max(which(FLR_PXD004939_MM$FLR_Unadjusted<=0.05)),]
FLR_PXD005241_MM_05 <- FLR_PXD005241_MM[1:max(which(FLR_PXD005241_MM$FLR_Unadjusted<=0.05)),]
FLR_PXD012764_MM_05 <- FLR_PXD012764_MM[1:max(which(FLR_PXD012764_MM$FLR_Unadjusted<=0.05)),]
FLR_PXD019291_MM_05 <- FLR_PXD019291_MM[1:max(which(FLR_PXD019291_MM$FLR_Unadjusted<=0.05)),]


FLR_AllRice_pform_MM_R1<-dplyr::bind_rows(FLR_PXD000923_MM_01, FLR_PXD002222_MM_01, FLR_PXD002756_MM_01, FLR_PXD004705_MM_01,
                                           FLR_PXD004939_MM_01, FLR_PXD005241_MM_01, FLR_PXD012764_MM_01, FLR_PXD019291_MM_01)

FLR_AllRice_pform_MM_R5<-dplyr::bind_rows(FLR_PXD000923_MM_05, FLR_PXD002222_MM_05, FLR_PXD002756_MM_05, FLR_PXD004705_MM_05,
                                           FLR_PXD004939_MM_05, FLR_PXD005241_MM_05, FLR_PXD012764_MM_05, FLR_PXD019291_MM_05)

Table_FLR_AllRice_pform_MM_R5 <- FLR_AllRice_pform_MM_R5 %>%
  dplyr::select(dataset, Rnumber) %>% 
  group_by(dataset) %>% 
  summarise(NumPform_Max = max(Rnumber))

Table_FLR_AllRice_pform_MM_R1 <- FLR_AllRice_pform_MM_R1 %>%
  dplyr::select(dataset, Rnumber) %>% 
  group_by(dataset) %>% 
  summarise(NumPform_Max = max(Rnumber))

# Checking numbers #


AllRice_PSM<-dplyr::bind_rows(PXD000923, PXD002222, PXD002756, PXD004705, PXD004939, PXD005241, PXD012764, PXD019291)

AllRice_PSM$Peptidoform <- paste0(AllRice_PSM$Peptide_mod,"_",AllRice_PSM$PTM.positions)


Unadjusted_Rice_Max <- aggregate(PTM_final_prob ~ dataset+Peptidoform, data = AllRice_PSM, max)

tab1(Unadjusted_Rice_Max$dataset)

# Binomial adjusted data collapsed using Max method #

library(plyr)

Bin_PXD000923 <- binAdjustPSM(PXD000923)
Bin_PXD002222 <- binAdjustPSM(PXD002222)
Bin_PXD002756 <- binAdjustPSM(PXD002756)
Bin_PXD004705 <- binAdjustPSM(PXD004705)
Bin_PXD004939 <- binAdjustPSM(PXD004939)
Bin_PXD005241 <- binAdjustPSM(PXD005241)
Bin_PXD012764 <- binAdjustPSM(PXD012764)
Bin_PXD019291 <- binAdjustPSM(PXD019291)

detach(package:plyr)

Bin_PXD000923_Max <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD000923, max)
Bin_PXD002222_Max <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD002222, max)
Bin_PXD002756_Max <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD002756, max)
Bin_PXD004705_Max <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD004705, max)
Bin_PXD004939_Max <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD004939, max)
Bin_PXD005241_Max <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD005241, max)
Bin_PXD012764_Max <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD012764, max)
Bin_PXD019291_Max <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD019291, max)


Bin_PXD000923_Max$Bin_Adjusted_Score <- Bin_PXD000923_Max$Bin_Adjusted_Score
Bin_PXD002222_Max$Bin_Adjusted_Score <- Bin_PXD002222_Max$Bin_Adjusted_Score
Bin_PXD002756_Max$Bin_Adjusted_Score <- Bin_PXD002756_Max$Bin_Adjusted_Score
Bin_PXD004705_Max$Bin_Adjusted_Score <- Bin_PXD004705_Max$Bin_Adjusted_Score
Bin_PXD004939_Max$Bin_Adjusted_Score <- Bin_PXD004939_Max$Bin_Adjusted_Score
Bin_PXD005241_Max$Bin_Adjusted_Score <- Bin_PXD005241_Max$Bin_Adjusted_Score
Bin_PXD012764_Max$Bin_Adjusted_Score <- Bin_PXD012764_Max$Bin_Adjusted_Score
Bin_PXD019291_Max$Bin_Adjusted_Score <- Bin_PXD019291_Max$Bin_Adjusted_Score



FLR_Bin_PXD000923_Max <- FLR_Adj(Bin_PXD000923_Max)
FLR_Bin_PXD002222_Max <- FLR_Adj(Bin_PXD002222_Max)
FLR_Bin_PXD002756_Max <- FLR_Adj(Bin_PXD002756_Max)
FLR_Bin_PXD004705_Max <- FLR_Adj(Bin_PXD004705_Max)
FLR_Bin_PXD004939_Max <- FLR_Adj(Bin_PXD004939_Max)
FLR_Bin_PXD005241_Max <- FLR_Adj(Bin_PXD005241_Max)
FLR_Bin_PXD012764_Max <- FLR_Adj(Bin_PXD012764_Max)
FLR_Bin_PXD019291_Max <- FLR_Adj(Bin_PXD019291_Max)


FLR_Bin_AllRice_pform_Max <-dplyr::bind_rows(FLR_Bin_PXD000923_Max,FLR_Bin_PXD002222_Max, FLR_Bin_PXD002756_Max,
                                         FLR_Bin_PXD004705_Max,FLR_Bin_PXD004939_Max, FLR_Bin_PXD005241_Max,
                                         FLR_Bin_PXD012764_Max,FLR_Bin_PXD019291_Max)



FLR_Bin_PXD000923_Max_01 <- FLR_Bin_PXD000923_Max[1:max(which(FLR_Bin_PXD000923_Max$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD002222_Max_01 <- FLR_Bin_PXD002222_Max[1:max(which(FLR_Bin_PXD002222_Max$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD002756_Max_01 <- FLR_Bin_PXD002756_Max[1:max(which(FLR_Bin_PXD002756_Max$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD004705_Max_01 <- FLR_Bin_PXD004705_Max[1:max(which(FLR_Bin_PXD004705_Max$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD004939_Max_01 <- FLR_Bin_PXD004939_Max[1:max(which(FLR_Bin_PXD004939_Max$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD005241_Max_01 <- FLR_Bin_PXD005241_Max[1:max(which(FLR_Bin_PXD005241_Max$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD012764_Max_01 <- FLR_Bin_PXD012764_Max[1:max(which(FLR_Bin_PXD012764_Max$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD019291_Max_01 <- FLR_Bin_PXD019291_Max[1:max(which(FLR_Bin_PXD019291_Max$FLR_Adj_Score<=0.01)),]

FLR_Bin_PXD000923_Max_05 <- FLR_Bin_PXD000923_Max[1:max(which(FLR_Bin_PXD000923_Max$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD002222_Max_05 <- FLR_Bin_PXD002222_Max[1:max(which(FLR_Bin_PXD002222_Max$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD002756_Max_05 <- FLR_Bin_PXD002756_Max[1:max(which(FLR_Bin_PXD002756_Max$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD004705_Max_05 <- FLR_Bin_PXD004705_Max[1:max(which(FLR_Bin_PXD004705_Max$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD004939_Max_05 <- FLR_Bin_PXD004939_Max[1:max(which(FLR_Bin_PXD004939_Max$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD005241_Max_05 <- FLR_Bin_PXD005241_Max[1:max(which(FLR_Bin_PXD005241_Max$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD012764_Max_05 <- FLR_Bin_PXD012764_Max[1:max(which(FLR_Bin_PXD012764_Max$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD019291_Max_05 <- FLR_Bin_PXD019291_Max[1:max(which(FLR_Bin_PXD019291_Max$FLR_Adj_Score<=0.05)),]


FLR_Bin_AllRice_pform_Max_R1<-dplyr::bind_rows(FLR_Bin_PXD000923_Max_01, FLR_Bin_PXD002222_Max_01, FLR_Bin_PXD002756_Max_01, FLR_Bin_PXD004705_Max_01,
                                           FLR_Bin_PXD004939_Max_01, FLR_Bin_PXD005241_Max_01, FLR_Bin_PXD012764_Max_01, FLR_Bin_PXD019291_Max_01)

FLR_Bin_AllRice_pform_Max_R5<-dplyr::bind_rows(FLR_Bin_PXD000923_Max_05, FLR_Bin_PXD002222_Max_05, FLR_Bin_PXD002756_Max_05, FLR_Bin_PXD004705_Max_05,
                                           FLR_Bin_PXD004939_Max_05, FLR_Bin_PXD005241_Max_05, FLR_Bin_PXD012764_Max_05, FLR_Bin_PXD019291_Max_05)



Table_FLR_Bin_AllRice_pform_Max_R5 <- FLR_Bin_AllRice_pform_Max_R5 %>%
  dplyr::select(dataset, Rnumber) %>% 
  group_by(dataset) %>% 
  summarize(NumPform_Max = max(Rnumber))

Table_FLR_Bin_AllRice_pform_Max_R1 <- FLR_Bin_AllRice_pform_Max_R1 %>%
  dplyr::select(dataset, Rnumber) %>% 
  group_by(dataset) %>% 
  summarize(NumPform_Max = max(Rnumber))

# Binomial adjusted data collapsed using MM method #

PXD000923_USI <- aggregate(PTM_final_prob ~ All_USI, data = PXD000923, mean)
PXD002222_USI <- aggregate(PTM_final_prob ~ All_USI, data = PXD002222, mean)
PXD002756_USI <- aggregate(PTM_final_prob ~ All_USI, data = PXD002756, mean)
PXD004705_USI <- aggregate(PTM_final_prob ~ All_USI, data = PXD004705, mean)
PXD004939_USI <- aggregate(PTM_final_prob ~ All_USI, data = PXD004939, mean)
PXD005241_USI <- aggregate(PTM_final_prob ~ All_USI, data = PXD005241, mean)
PXD012764_USI <- aggregate(PTM_final_prob ~ All_USI, data = PXD012764, mean)
PXD019291_USI <- aggregate(PTM_final_prob ~ All_USI, data = PXD019291, mean)

PXD000923_USI2 <- subset(PXD000923,select=-c(PTM_final_prob))
PXD002222_USI2 <- subset(PXD002222,select=-c(PTM_final_prob))
PXD002756_USI2 <- subset(PXD002756,select=-c(PTM_final_prob))
PXD004705_USI2 <- subset(PXD004705,select=-c(PTM_final_prob))
PXD004939_USI2 <- subset(PXD004939,select=-c(PTM_final_prob))
PXD005241_USI2 <- subset(PXD005241,select=-c(PTM_final_prob))
PXD012764_USI2 <- subset(PXD012764,select=-c(PTM_final_prob))
PXD019291_USI2 <- subset(PXD019291,select=-c(PTM_final_prob))

PXD000923_USI3 <- merge(PXD000923_USI,PXD000923_USI2, by="All_USI")
PXD002222_USI3 <- merge(PXD002222_USI,PXD002222_USI2, by="All_USI")
PXD002756_USI3 <- merge(PXD002756_USI,PXD002756_USI2, by="All_USI")
PXD004705_USI3 <- merge(PXD004705_USI,PXD004705_USI2, by="All_USI")
PXD004939_USI3 <- merge(PXD004939_USI,PXD004939_USI2, by="All_USI")
PXD005241_USI3 <- merge(PXD005241_USI,PXD005241_USI2, by="All_USI")
PXD012764_USI3 <- merge(PXD012764_USI,PXD012764_USI2, by="All_USI")
PXD019291_USI3 <- merge(PXD019291_USI,PXD019291_USI2, by="All_USI")

library(plyr)

Bin_PXD000923_USI <- binAdjustPSM(PXD000923_USI3)
Bin_PXD002222_USI <- binAdjustPSM(PXD002222_USI3)
Bin_PXD002756_USI <- binAdjustPSM(PXD002756_USI3)
Bin_PXD004705_USI <- binAdjustPSM(PXD004705_USI3)
Bin_PXD004939_USI <- binAdjustPSM(PXD004939_USI3)
Bin_PXD005241_USI <- binAdjustPSM(PXD005241_USI3)
Bin_PXD012764_USI <- binAdjustPSM(PXD012764_USI3)
Bin_PXD019291_USI <- binAdjustPSM(PXD019291_USI3)

detach(package:plyr)


Bin_PXD000923_MM <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD000923_USI, max)
Bin_PXD002222_MM <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD002222_USI, max)
Bin_PXD002756_MM <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD002756_USI, max)
Bin_PXD004705_MM <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD004705_USI, max)
Bin_PXD004939_MM <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD004939_USI, max)
Bin_PXD005241_MM <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD005241_USI, max)
Bin_PXD012764_MM <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD012764_USI, max)
Bin_PXD019291_MM <- aggregate(Bin_Adjusted_Score ~ Peptidoform+dataset+Peptide+Peptide_mod+PTM.positions, data = Bin_PXD019291_USI, max)


FLR_Bin_PXD000923_MM <- FLR_Adj(Bin_PXD000923_MM)
FLR_Bin_PXD002222_MM <- FLR_Adj(Bin_PXD002222_MM)
FLR_Bin_PXD002756_MM <- FLR_Adj(Bin_PXD002756_MM)
FLR_Bin_PXD004705_MM <- FLR_Adj(Bin_PXD004705_MM)
FLR_Bin_PXD004939_MM <- FLR_Adj(Bin_PXD004939_MM)
FLR_Bin_PXD005241_MM <- FLR_Adj(Bin_PXD005241_MM)
FLR_Bin_PXD012764_MM <- FLR_Adj(Bin_PXD012764_MM)
FLR_Bin_PXD019291_MM <- FLR_Adj(Bin_PXD019291_MM)

FLR_Bin_AllRice_pform_MM <-dplyr::bind_rows(FLR_Bin_PXD000923_MM,FLR_Bin_PXD002222_MM, FLR_Bin_PXD002756_MM,
                                        FLR_Bin_PXD004705_MM,FLR_Bin_PXD004939_MM, FLR_Bin_PXD005241_MM,
                                        FLR_Bin_PXD012764_MM,FLR_Bin_PXD019291_MM)


FLR_Bin_PXD000923_MM_01 <- FLR_Bin_PXD000923_MM[1:max(which(FLR_Bin_PXD000923_MM$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD002222_MM_01 <- FLR_Bin_PXD002222_MM[1:max(which(FLR_Bin_PXD002222_MM$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD002756_MM_01 <- FLR_Bin_PXD002756_MM[1:max(which(FLR_Bin_PXD002756_MM$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD004705_MM_01 <- FLR_Bin_PXD004705_MM[1:max(which(FLR_Bin_PXD004705_MM$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD004939_MM_01 <- FLR_Bin_PXD004939_MM[1:max(which(FLR_Bin_PXD004939_MM$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD005241_MM_01 <- FLR_Bin_PXD005241_MM[1:max(which(FLR_Bin_PXD005241_MM$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD012764_MM_01 <- FLR_Bin_PXD012764_MM[1:max(which(FLR_Bin_PXD012764_MM$FLR_Adj_Score<=0.01)),]
FLR_Bin_PXD019291_MM_01 <- FLR_Bin_PXD019291_MM[1:max(which(FLR_Bin_PXD019291_MM$FLR_Adj_Score<=0.01)),]

FLR_Bin_PXD000923_MM_05 <- FLR_Bin_PXD000923_MM[1:max(which(FLR_Bin_PXD000923_MM$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD002222_MM_05 <- FLR_Bin_PXD002222_MM[1:max(which(FLR_Bin_PXD002222_MM$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD002756_MM_05 <- FLR_Bin_PXD002756_MM[1:max(which(FLR_Bin_PXD002756_MM$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD004705_MM_05 <- FLR_Bin_PXD004705_MM[1:max(which(FLR_Bin_PXD004705_MM$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD004939_MM_05 <- FLR_Bin_PXD004939_MM[1:max(which(FLR_Bin_PXD004939_MM$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD005241_MM_05 <- FLR_Bin_PXD005241_MM[1:max(which(FLR_Bin_PXD005241_MM$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD012764_MM_05 <- FLR_Bin_PXD012764_MM[1:max(which(FLR_Bin_PXD012764_MM$FLR_Adj_Score<=0.05)),]
FLR_Bin_PXD019291_MM_05 <- FLR_Bin_PXD019291_MM[1:max(which(FLR_Bin_PXD019291_MM$FLR_Adj_Score<=0.05)),]


FLR_Bin_AllRice_pform_MM_R1<-dplyr::bind_rows(FLR_Bin_PXD000923_MM_01, FLR_Bin_PXD002222_MM_01, FLR_Bin_PXD002756_MM_01, FLR_Bin_PXD004705_MM_01,
                                               FLR_Bin_PXD004939_MM_01, FLR_Bin_PXD005241_MM_01, FLR_Bin_PXD012764_MM_01, FLR_Bin_PXD019291_MM_01)

FLR_Bin_AllRice_pform_MM_R5<-dplyr::bind_rows(FLR_Bin_PXD000923_MM_05, FLR_Bin_PXD002222_MM_05, FLR_Bin_PXD002756_MM_05, FLR_Bin_PXD004705_MM_05,
                                               FLR_Bin_PXD004939_MM_05, FLR_Bin_PXD005241_MM_05, FLR_Bin_PXD012764_MM_05, FLR_Bin_PXD019291_MM_05)

Table_FLR_Bin_AllRice_pform_MM_R5 <- FLR_Bin_AllRice_pform_MM_R5 %>%
  dplyr::select(dataset, Rnumber) %>% 
  group_by(dataset) %>% 
  summarise(NumPform_Max = max(Rnumber))

Table_FLR_Bin_AllRice_pform_MM_R1 <- FLR_Bin_AllRice_pform_MM_R1 %>%
  dplyr::select(dataset, Rnumber) %>% 
  group_by(dataset) %>% 
  summarise(NumPform_Max = max(Rnumber))

