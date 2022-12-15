

# Table 5

source('D:/Paper submission/R Functions/Binomial Scores Function Peptidoform level.R')
source('D:/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('D:/Paper submission/R Functions/Function frequency of site.R')
source('D:/Paper submission/R Functions/Traditional Analysis.R')
source('D:/Paper submission/R Functions/GSB_Function.R')

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
library(gmodels)


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


tab1(Adjusted_PXD000923$Amino)
tab1(Adjusted_PXD002222$Amino)
tab1(Adjusted_PXD002756$Amino)
tab1(Adjusted_PXD004705$Amino)
tab1(Adjusted_PXD004939$Amino)
tab1(Adjusted_PXD005241$Amino)
tab1(Adjusted_PXD012764$Amino)
tab1(Adjusted_PXD019291$Amino)

N_Prot_Adjusted_PXD000923 <- Adjusted_PXD000923[!duplicated(Adjusted_PXD000923[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD002222 <- Adjusted_PXD002222[!duplicated(Adjusted_PXD002222[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD002756 <- Adjusted_PXD002756[!duplicated(Adjusted_PXD002756[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD004705 <- Adjusted_PXD004705[!duplicated(Adjusted_PXD004705[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD004939 <- Adjusted_PXD004939[!duplicated(Adjusted_PXD004939[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD005241 <- Adjusted_PXD005241[!duplicated(Adjusted_PXD005241[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD012764 <- Adjusted_PXD012764[!duplicated(Adjusted_PXD012764[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD019291 <- Adjusted_PXD019291[!duplicated(Adjusted_PXD019291[c("PROTEIN_LOC")]),]


AllRice_Peptidoform<-dplyr::bind_rows(Adjusted_PXD000923, Adjusted_PXD002222, Adjusted_PXD002756, Adjusted_PXD004705,
                                      Adjusted_PXD004939, Adjusted_PXD005241, Adjusted_PXD012764, Adjusted_PXD019291)


FLR_PXD000923_Peptidoform_05 <- Adjusted_PXD000923[1:max(which(Adjusted_PXD000923$FLR_Adj_Score<=0.05)),]
FLR_PXD002222_Peptidoform_05 <- Adjusted_PXD002222[1:max(which(Adjusted_PXD002222$FLR_Adj_Score<=0.05)),]
FLR_PXD002756_Peptidoform_05 <- Adjusted_PXD002756[1:max(which(Adjusted_PXD002756$FLR_Adj_Score<=0.05)),]
FLR_PXD004705_Peptidoform_05 <- Adjusted_PXD004705[1:max(which(Adjusted_PXD004705$FLR_Adj_Score<=0.05)),]
FLR_PXD004939_Peptidoform_05 <- Adjusted_PXD004939[1:max(which(Adjusted_PXD004939$FLR_Adj_Score<=0.05)),]
FLR_PXD005241_Peptidoform_05 <- Adjusted_PXD005241[1:max(which(Adjusted_PXD005241$FLR_Adj_Score<=0.05)),]
FLR_PXD012764_Peptidoform_05 <- Adjusted_PXD012764[1:max(which(Adjusted_PXD012764$FLR_Adj_Score<=0.05)),]
FLR_PXD019291_Peptidoform_05 <- Adjusted_PXD019291[1:max(which(Adjusted_PXD019291$FLR_Adj_Score<=0.05)),]



FLR_PXD000923_Peptidoform_01 <- Adjusted_PXD000923[1:max(which(Adjusted_PXD000923$FLR_Adj_Score<=0.01)),]
FLR_PXD002222_Peptidoform_01 <- Adjusted_PXD002222[1:max(which(Adjusted_PXD002222$FLR_Adj_Score<=0.01)),]
FLR_PXD002756_Peptidoform_01 <- Adjusted_PXD002756[1:max(which(Adjusted_PXD002756$FLR_Adj_Score<=0.01)),]
FLR_PXD004705_Peptidoform_01 <- Adjusted_PXD004705[1:max(which(Adjusted_PXD004705$FLR_Adj_Score<=0.01)),]
FLR_PXD004939_Peptidoform_01 <- Adjusted_PXD004939[1:max(which(Adjusted_PXD004939$FLR_Adj_Score<=0.01)),]
FLR_PXD005241_Peptidoform_01 <- Adjusted_PXD005241[1:max(which(Adjusted_PXD005241$FLR_Adj_Score<=0.01)),]
FLR_PXD012764_Peptidoform_01 <- Adjusted_PXD012764[1:max(which(Adjusted_PXD012764$FLR_Adj_Score<=0.01)),]
FLR_PXD019291_Peptidoform_01 <- Adjusted_PXD019291[1:max(which(Adjusted_PXD019291$FLR_Adj_Score<=0.01)),]

FLR_Peptidoform_01<-dplyr::bind_rows(FLR_PXD000923_Peptidoform_01, FLR_PXD002222_Peptidoform_01, FLR_PXD002756_Peptidoform_01, FLR_PXD004705_Peptidoform_01,
                                      FLR_PXD004939_Peptidoform_01, FLR_PXD005241_Peptidoform_01, FLR_PXD012764_Peptidoform_01, FLR_PXD019291_Peptidoform_01)

FLR_Peptidoform_05<-dplyr::bind_rows(FLR_PXD000923_Peptidoform_05, FLR_PXD002222_Peptidoform_05, FLR_PXD002756_Peptidoform_05, FLR_PXD004705_Peptidoform_05,
                                     FLR_PXD004939_Peptidoform_05, FLR_PXD005241_Peptidoform_05, FLR_PXD012764_Peptidoform_05, FLR_PXD019291_Peptidoform_05)


N_Prot_Adjusted_PXD000923_01 <- FLR_PXD000923_Peptidoform_01[!duplicated(FLR_PXD000923_Peptidoform_01[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD002222_01 <- FLR_PXD002222_Peptidoform_01[!duplicated(FLR_PXD002222_Peptidoform_01[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD002756_01 <- FLR_PXD002756_Peptidoform_01[!duplicated(FLR_PXD002756_Peptidoform_01[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD004705_01 <- FLR_PXD004705_Peptidoform_01[!duplicated(FLR_PXD004705_Peptidoform_01[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD004939_01 <- FLR_PXD004939_Peptidoform_01[!duplicated(FLR_PXD004939_Peptidoform_01[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD005241_01 <- FLR_PXD005241_Peptidoform_01[!duplicated(FLR_PXD005241_Peptidoform_01[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD012764_01 <- FLR_PXD012764_Peptidoform_01[!duplicated(FLR_PXD012764_Peptidoform_01[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD019291_01 <- FLR_PXD019291_Peptidoform_01[!duplicated(FLR_PXD019291_Peptidoform_01[c("PROTEIN_LOC")]),]


N_Prot_Adjusted_PXD000923_05 <- FLR_PXD000923_Peptidoform_05[!duplicated(FLR_PXD000923_Peptidoform_05[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD002222_05 <- FLR_PXD002222_Peptidoform_05[!duplicated(FLR_PXD002222_Peptidoform_05[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD002756_05 <- FLR_PXD002756_Peptidoform_05[!duplicated(FLR_PXD002756_Peptidoform_05[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD004705_05 <- FLR_PXD004705_Peptidoform_05[!duplicated(FLR_PXD004705_Peptidoform_05[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD004939_05 <- FLR_PXD004939_Peptidoform_05[!duplicated(FLR_PXD004939_Peptidoform_05[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD005241_05 <- FLR_PXD005241_Peptidoform_05[!duplicated(FLR_PXD005241_Peptidoform_05[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD012764_05 <- FLR_PXD012764_Peptidoform_05[!duplicated(FLR_PXD012764_Peptidoform_05[c("PROTEIN_LOC")]),]
N_Prot_Adjusted_PXD019291_05 <- FLR_PXD019291_Peptidoform_05[!duplicated(FLR_PXD019291_Peptidoform_05[c("PROTEIN_LOC")]),]

AllRice_pform_Max<-dplyr::bind_rows(N_Prot_Adjusted_PXD000923, N_Prot_Adjusted_PXD002222, N_Prot_Adjusted_PXD002756, N_Prot_Adjusted_PXD004705,
                                       N_Prot_Adjusted_PXD004939, N_Prot_Adjusted_PXD005241, N_Prot_Adjusted_PXD012764, N_Prot_Adjusted_PXD019291)

AllRice_pform_Max_01<-dplyr::bind_rows(N_Prot_Adjusted_PXD000923_01, N_Prot_Adjusted_PXD002222_01, N_Prot_Adjusted_PXD002756_01, N_Prot_Adjusted_PXD004705_01,
                                       N_Prot_Adjusted_PXD004939_01, N_Prot_Adjusted_PXD005241_01, N_Prot_Adjusted_PXD012764_01, N_Prot_Adjusted_PXD019291_01)

AllRice_pform_Max_05<-dplyr::bind_rows(N_Prot_Adjusted_PXD000923_05, N_Prot_Adjusted_PXD002222_05, N_Prot_Adjusted_PXD002756_05, N_Prot_Adjusted_PXD004705_05,
                                       N_Prot_Adjusted_PXD004939_05, N_Prot_Adjusted_PXD005241_05, N_Prot_Adjusted_PXD012764_05, N_Prot_Adjusted_PXD019291_05)
tab1(AllRice_pform_Max$dataset)
tab1(AllRice_pform_Max_01$dataset)
tab1(AllRice_pform_Max_05$dataset)

Unique <- AllRice_Peptidoform[!duplicated(AllRice_Peptidoform[c("PROTEIN_LOC")]),]
Unique_01 <- AllRice_pform_Max_01[!duplicated(AllRice_pform_Max_01[c("PROTEIN_LOC")]),]
Unique_05 <- AllRice_pform_Max_05[!duplicated(AllRice_pform_Max_05[c("PROTEIN_LOC")]),]

AllRice_pform_Max_Final <- GSB_Function(FLR_Peptidoform_01,FLR_Peptidoform_05)
CrossTable(AllRice_pform_Max_Final$cat, AllRice_pform_Max_Final$Amino)

# Traditional analysis based on scores #
########################################


PXD000923_TA <- TA_Pform(PXD000923)
PXD002222_TA <- TA_Pform(PXD002222)
PXD002756_TA <- TA_Pform(PXD002756)
PXD004705_TA <- TA_Pform(PXD004705)
PXD004939_TA <- TA_Pform(PXD004939)
PXD005241_TA <- TA_Pform(PXD005241)
PXD012764_TA <- TA_Pform(PXD012764)
PXD019291_TA <- TA_Pform(PXD019291)


PXD000923_Peptidoform_99 <- PXD000923_TA[PXD000923_TA$PTM_final_prob>=0.99,]
PXD002222_Peptidoform_99 <- PXD002222_TA[PXD002222_TA$PTM_final_prob>=0.99,]
PXD002756_Peptidoform_99 <- PXD002756_TA[PXD002756_TA$PTM_final_prob>=0.99,]
PXD004705_Peptidoform_99 <- PXD004705_TA[PXD004705_TA$PTM_final_prob>=0.99,]
PXD004939_Peptidoform_99 <- PXD004939_TA[PXD004939_TA$PTM_final_prob>=0.99,]
PXD005241_Peptidoform_99 <- PXD005241_TA[PXD005241_TA$PTM_final_prob>=0.99,]
PXD012764_Peptidoform_99 <- PXD012764_TA[PXD012764_TA$PTM_final_prob>=0.99,]
PXD019291_Peptidoform_99 <- PXD019291_TA[PXD019291_TA$PTM_final_prob>=0.99,]


PXD000923_Peptidoform_95 <- PXD000923_TA[PXD000923_TA$PTM_final_prob>=0.95,]
PXD002222_Peptidoform_95 <- PXD002222_TA[PXD002222_TA$PTM_final_prob>=0.95,]
PXD002756_Peptidoform_95 <- PXD002756_TA[PXD002756_TA$PTM_final_prob>=0.95,]
PXD004705_Peptidoform_95 <- PXD004705_TA[PXD004705_TA$PTM_final_prob>=0.95,]
PXD004939_Peptidoform_95 <- PXD004939_TA[PXD004939_TA$PTM_final_prob>=0.95,]
PXD005241_Peptidoform_95 <- PXD005241_TA[PXD005241_TA$PTM_final_prob>=0.95,]
PXD012764_Peptidoform_95 <- PXD012764_TA[PXD012764_TA$PTM_final_prob>=0.95,]
PXD019291_Peptidoform_95 <- PXD019291_TA[PXD019291_TA$PTM_final_prob>=0.95,]



AllRice_Pform_TA<-dplyr::bind_rows(PXD000923_TA, PXD002222_TA,PXD002756_TA, PXD004705_TA,
                                     PXD004939_TA, PXD005241_TA, PXD012764_TA, PXD019291_TA)


N_Prot_PXD000923_TA_99 <- PXD000923_Peptidoform_99[!duplicated(PXD000923_Peptidoform_99[c("PROTEIN_LOC")]),]
N_Prot_PXD002222_TA_99 <- PXD002222_Peptidoform_99[!duplicated(PXD002222_Peptidoform_99[c("PROTEIN_LOC")]),]
N_Prot_PXD002756_TA_99 <- PXD002756_Peptidoform_99[!duplicated(PXD002756_Peptidoform_99[c("PROTEIN_LOC")]),]
N_Prot_PXD004705_TA_99 <- PXD004705_Peptidoform_99[!duplicated(PXD004705_Peptidoform_99[c("PROTEIN_LOC")]),]
N_Prot_PXD004939_TA_99 <- PXD004939_Peptidoform_99[!duplicated(PXD004939_Peptidoform_99[c("PROTEIN_LOC")]),]
N_Prot_PXD005241_TA_99 <- PXD005241_Peptidoform_99[!duplicated(PXD005241_Peptidoform_99[c("PROTEIN_LOC")]),]
N_Prot_PXD012764_TA_99 <- PXD012764_Peptidoform_99[!duplicated(PXD012764_Peptidoform_99[c("PROTEIN_LOC")]),]
N_Prot_PXD019291_TA_99 <- PXD019291_Peptidoform_99[!duplicated(PXD019291_Peptidoform_99[c("PROTEIN_LOC")]),]



N_Prot_PXD000923_TA_95 <- PXD000923_Peptidoform_95[!duplicated(PXD000923_Peptidoform_95[c("PROTEIN_LOC")]),]
N_Prot_PXD002222_TA_95 <- PXD002222_Peptidoform_95[!duplicated(PXD002222_Peptidoform_95[c("PROTEIN_LOC")]),]
N_Prot_PXD002756_TA_95 <- PXD002756_Peptidoform_95[!duplicated(PXD002756_Peptidoform_95[c("PROTEIN_LOC")]),]
N_Prot_PXD004705_TA_95 <- PXD004705_Peptidoform_95[!duplicated(PXD004705_Peptidoform_95[c("PROTEIN_LOC")]),]
N_Prot_PXD004939_TA_95 <- PXD004939_Peptidoform_95[!duplicated(PXD004939_Peptidoform_95[c("PROTEIN_LOC")]),]
N_Prot_PXD005241_TA_95 <- PXD005241_Peptidoform_95[!duplicated(PXD005241_Peptidoform_95[c("PROTEIN_LOC")]),]
N_Prot_PXD012764_TA_95 <- PXD012764_Peptidoform_95[!duplicated(PXD012764_Peptidoform_95[c("PROTEIN_LOC")]),]
N_Prot_PXD019291_TA_95 <- PXD019291_Peptidoform_95[!duplicated(PXD019291_Peptidoform_95[c("PROTEIN_LOC")]),]


AllRice_TA_99<-dplyr::bind_rows(N_Prot_PXD000923_TA_99, N_Prot_PXD002222_TA_99, N_Prot_PXD002756_TA_99, N_Prot_PXD004705_TA_99,
                                       N_Prot_PXD004939_TA_99, N_Prot_PXD005241_TA_99, N_Prot_PXD012764_TA_99, N_Prot_PXD019291_TA_99)

AllRice_TA_95<-dplyr::bind_rows(N_Prot_PXD000923_TA_95, N_Prot_PXD002222_TA_95, N_Prot_PXD002756_TA_95, N_Prot_PXD004705_TA_95,
                                N_Prot_PXD004939_TA_95, N_Prot_PXD005241_TA_95, N_Prot_PXD012764_TA_95, N_Prot_PXD019291_TA_95)



tab1(AllRice_Pform_TA$dataset)
tab1(AllRice_TA_99$dataset)
tab1(AllRice_TA_95$dataset)


Unique_AllRice_TA_99 <- AllRice_TA_99[!duplicated(AllRice_TA_99[c("PROTEIN_LOC")]),]
Unique_AllRice_TA_95 <- AllRice_TA_95[!duplicated(AllRice_TA_95[c("PROTEIN_LOC")]),]


