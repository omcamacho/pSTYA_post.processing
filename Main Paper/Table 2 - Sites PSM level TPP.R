
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
library(epiDisplay)
library(gmodels)


source('D:/Paper submission/R Functions/Binomial Scores Function PSM level.R')
source('D:/Paper submission/R Functions/Binomial Scores Function PSM level Syn7058.R')
source('D:/Paper submission/R Functions/FLR function Unadjusted.R')
source('D:/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('D:/Paper submission/R Functions/Function frequency of site.R')

# Arabidopsis #
###############


Ara_PXD008355 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD008355.csv')

Ara_PXD008355 <- dplyr::rename(Ara_PXD008355, PTM.positions = PTM_positions)

Ara_PXD008355 <- dplyr::rename(Ara_PXD008355, PTM.Score = PTM_score)

BinAdj_Ara_PSM<-binAdjustPSM(Ara_PXD008355)

Not_adj_Ara_PSM <- cbind.data.frame(BinAdj_Ara_PSM$Peptide_mod,BinAdj_Ara_PSM$Peptide,BinAdj_Ara_PSM$PTM.positions,BinAdj_Ara_PSM$PTM_final_prob)
names(Not_adj_Ara_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Unadjusted_Score")

Bin_adj_Ara_PSM <- cbind.data.frame(BinAdj_Ara_PSM$Peptide_mod,BinAdj_Ara_PSM$Peptide,BinAdj_Ara_PSM$PTM.positions,BinAdj_Ara_PSM$Bin_Adjusted_Score)
names(Bin_adj_Ara_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Bin_Adjusted_Score")

FLR_Not_Adj_Ara_PSM<-FLR_NotAdj(Not_adj_Ara_PSM)

FLR_Not_Adj_Ara_PSM_R10 <- FLR_Not_Adj_Ara_PSM[1:max(which(FLR_Not_Adj_Ara_PSM$FLR_Unadjusted<=0.1)),]
FLR_Not_Adj_Ara_PSM_R5 <- FLR_Not_Adj_Ara_PSM[1:max(which(FLR_Not_Adj_Ara_PSM$FLR_Unadjusted<=0.05)),]
FLR_Not_Adj_Ara_PSM_R2.5 <- FLR_Not_Adj_Ara_PSM[1:max(which(FLR_Not_Adj_Ara_PSM$FLR_Unadjusted<=0.025)),]
FLR_Not_Adj_Ara_PSM_R1 <- FLR_Not_Adj_Ara_PSM[1:max(which(FLR_Not_Adj_Ara_PSM$FLR_Unadjusted<=0.01)),]

FLR_Adj_Ara_PSM<-FLR_Adj(Bin_adj_Ara_PSM)

FLR_Adj_Ara_PSM_R10 <- FLR_Adj_Ara_PSM[1:max(which(FLR_Adj_Ara_PSM$FLR_Adj_Score<=0.1)),]
FLR_Adj_Ara_PSM_R5 <- FLR_Adj_Ara_PSM[1:max(which(FLR_Adj_Ara_PSM$FLR_Adj_Score<=0.05)),]
FLR_Adj_Ara_PSM_R2.5 <- FLR_Adj_Ara_PSM[1:max(which(FLR_Adj_Ara_PSM$FLR_Adj_Score<=0.025)),]
FLR_Adj_Ara_PSM_R1 <- FLR_Adj_Ara_PSM[1:max(which(FLR_Adj_Ara_PSM$FLR_Adj_Score<=0.01)),]


# Human #
#########


Hum_PXD000612 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD000612.csv')

Hum_PXD000612 <- dplyr::rename(Hum_PXD000612, PTM.positions = PTM_positions)

Hum_PXD000612 <- dplyr::rename(Hum_PXD000612, PTM.Score = PTM_score)

BinAdj_Hum_PSM<-binAdjustPSM(Hum_PXD000612)

Not_adj_Hum_PSM <- cbind.data.frame(BinAdj_Hum_PSM$Peptide_mod,BinAdj_Hum_PSM$Peptide,BinAdj_Hum_PSM$PTM.positions,BinAdj_Hum_PSM$PTM_final_prob)
names(Not_adj_Hum_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Unadjusted_Score")

Bin_adj_Hum_PSM <- cbind.data.frame(BinAdj_Hum_PSM$Peptide_mod,BinAdj_Hum_PSM$Peptide,BinAdj_Hum_PSM$PTM.positions,BinAdj_Hum_PSM$Bin_Adjusted_Score)
names(Bin_adj_Hum_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Bin_Adjusted_Score")

FLR_Not_Adj_Hum_PSM<-FLR_NotAdj(Not_adj_Hum_PSM)

FLR_Not_Adj_Hum_PSM_R10 <- FLR_Not_Adj_Hum_PSM[1:max(which(FLR_Not_Adj_Hum_PSM$FLR_Unadjusted<=0.1)),]
FLR_Not_Adj_Hum_PSM_R5 <- FLR_Not_Adj_Hum_PSM[1:max(which(FLR_Not_Adj_Hum_PSM$FLR_Unadjusted<=0.05)),]
FLR_Not_Adj_Hum_PSM_R2.5 <- FLR_Not_Adj_Hum_PSM[1:max(which(FLR_Not_Adj_Hum_PSM$FLR_Unadjusted<=0.025)),]
FLR_Not_Adj_Hum_PSM_R1 <- FLR_Not_Adj_Hum_PSM[1:max(which(FLR_Not_Adj_Hum_PSM$FLR_Unadjusted<=0.01)),]

FLR_Adj_Hum_PSM<-FLR_Adj(Bin_adj_Hum_PSM)

FLR_Adj_Hum_PSM_R10 <- FLR_Adj_Hum_PSM[1:max(which(FLR_Adj_Hum_PSM$FLR_Adj_Score<=0.1)),]
FLR_Adj_Hum_PSM_R5 <- FLR_Adj_Hum_PSM[1:max(which(FLR_Adj_Hum_PSM$FLR_Adj_Score<=0.05)),]
FLR_Adj_Hum_PSM_R2.5 <- FLR_Adj_Hum_PSM[1:max(which(FLR_Adj_Hum_PSM$FLR_Adj_Score<=0.025)),]
FLR_Adj_Hum_PSM_R1 <- FLR_Adj_Hum_PSM[1:max(which(FLR_Adj_Hum_PSM$FLR_Adj_Score<=0.01)),]


# Synthetic PXD000138 #
#######################


Syn_PXD000138 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD000138.csv')

Syn_PXD000138 <- dplyr::rename(Syn_PXD000138, PTM.positions = PTM_positions)

Syn_PXD000138 <- dplyr::rename(Syn_PXD000138, PTM.Score = PTM_score)

BinAdj_Syn138_PSM<-binAdjustPSM(Syn_PXD000138)

Not_adj_Syn138_PSM <- cbind.data.frame(BinAdj_Syn138_PSM$Peptide_mod,BinAdj_Syn138_PSM$Peptide,BinAdj_Syn138_PSM$PTM.positions,BinAdj_Syn138_PSM$PTM_final_prob, BinAdj_Syn138_PSM$Incorrect.count)
names(Not_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Unadjusted_Score", "Incorrect.count")

Bin_adj_Syn138_PSM <- cbind.data.frame(BinAdj_Syn138_PSM$Peptide_mod,BinAdj_Syn138_PSM$Peptide,BinAdj_Syn138_PSM$PTM.positions,BinAdj_Syn138_PSM$Bin_Adjusted_Score, BinAdj_Syn138_PSM$Incorrect.count)
names(Bin_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Bin_Adjusted_Score", "Incorrect.count")

FLR_Not_Adj_Syn138_PSM<-FLR_NotAdj(Not_adj_Syn138_PSM)

FLR_Not_Adj_Syn138_PSM_R10 <- FLR_Not_Adj_Syn138_PSM[1:max(which(FLR_Not_Adj_Syn138_PSM$FLR_Unadjusted<=0.1)),]
FLR_Not_Adj_Syn138_PSM_R5 <- FLR_Not_Adj_Syn138_PSM[1:max(which(FLR_Not_Adj_Syn138_PSM$FLR_Unadjusted<=0.05)),]
FLR_Not_Adj_Syn138_PSM_R2.5 <- FLR_Not_Adj_Syn138_PSM[1:max(which(FLR_Not_Adj_Syn138_PSM$FLR_Unadjusted<=0.025)),]
FLR_Not_Adj_Syn138_PSM_R1 <- FLR_Not_Adj_Syn138_PSM[1:max(which(FLR_Not_Adj_Syn138_PSM$FLR_Unadjusted<=0.01)),]

tab1(FLR_Not_Adj_Syn138_PSM_R10$Incorrect.count)
tab1(FLR_Not_Adj_Syn138_PSM_R5$Incorrect.count)
tab1(FLR_Not_Adj_Syn138_PSM_R2.5$Incorrect.count)
tab1(FLR_Not_Adj_Syn138_PSM_R1$Incorrect.count)

FLR_Adj_Syn138_PSM<-FLR_Adj(Bin_adj_Syn138_PSM)

FLR_Adj_Syn138_PSM_R10 <- FLR_Adj_Syn138_PSM[1:max(which(FLR_Adj_Syn138_PSM$FLR_Adj_Score<=0.1)),]
FLR_Adj_Syn138_PSM_R5 <- FLR_Adj_Syn138_PSM[1:max(which(FLR_Adj_Syn138_PSM$FLR_Adj_Score<=0.05)),]
FLR_Adj_Syn138_PSM_R2.5 <- FLR_Adj_Syn138_PSM[1:max(which(FLR_Adj_Syn138_PSM$FLR_Adj_Score<=0.025)),]
FLR_Adj_Syn138_PSM_R1 <- FLR_Not_Adj_Syn138_PSM[1:max(which(FLR_Adj_Syn138_PSM$FLR_Adj_Score<=0.01)),]

tab1(FLR_Adj_Syn138_PSM_R10$Incorrect.count)
tab1(FLR_Adj_Syn138_PSM_R5$Incorrect.count)
tab1(FLR_Adj_Syn138_PSM_R2.5$Incorrect.count)
tab1(FLR_Adj_Syn138_PSM_R1$Incorrect.count)


# Synthetic PXD007058 #
#######################


Syn_PXD007058 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD007058.csv')
Syn_PXD007058$Peptidoform <- paste0(Syn_PXD007058$Peptide_mod,"_",Syn_PXD007058$PTM_positions,"_",Syn_PXD007058$Pool)
Syn_PXD007058$All_Proteins <- Syn_PXD007058$Protein
Syn_PXD007058$Spectrum <- Syn_PXD007058$All_USI

BinAdj_Syn7058_PSM<-binAdjustSyn(Syn_PXD007058)

Not_adj_Syn7058_PSM <- cbind.data.frame(BinAdj_Syn7058_PSM$Peptide_mod,BinAdj_Syn7058_PSM$Peptide,BinAdj_Syn7058_PSM$PTM_positions,BinAdj_Syn7058_PSM$PTM_final_prob,BinAdj_Syn7058_PSM$Synthetic_full_match_False)
names(Not_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Unadjusted_Score", "Incorrect.count")

Bin_adj_Syn7058_PSM <- cbind.data.frame(BinAdj_Syn7058_PSM$Peptide_mod,BinAdj_Syn7058_PSM$Peptide,BinAdj_Syn7058_PSM$PTM_positions,BinAdj_Syn7058_PSM$Bin_Adjusted_Score,BinAdj_Syn7058_PSM$Synthetic_full_match_False)
names(Bin_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Bin_Adjusted_Score", "Incorrect.count")


FLR_Not_Adj_Syn7058_PSM<-FLR_NotAdj(Not_adj_Syn7058_PSM)

FLR_Not_Adj_Syn7058_PSM_R10 <- FLR_Not_Adj_Syn7058_PSM[1:max(which(FLR_Not_Adj_Syn7058_PSM$FLR_Unadjusted<=0.1)),]
FLR_Not_Adj_Syn7058_PSM_R5 <- FLR_Not_Adj_Syn7058_PSM[1:max(which(FLR_Not_Adj_Syn7058_PSM$FLR_Unadjusted<=0.05)),]
FLR_Not_Adj_Syn7058_PSM_R2.5 <- FLR_Not_Adj_Syn7058_PSM[1:max(which(FLR_Not_Adj_Syn7058_PSM$FLR_Unadjusted<=0.025)),]
FLR_Not_Adj_Syn7058_PSM_R1 <- FLR_Not_Adj_Syn7058_PSM[1:max(which(FLR_Not_Adj_Syn7058_PSM$FLR_Unadjusted<=0.01)),]

tab1(FLR_Not_Adj_Syn7058_PSM_R10$Incorrect.count)
tab1(FLR_Not_Adj_Syn7058_PSM_R5$Incorrect.count)
tab1(FLR_Not_Adj_Syn7058_PSM_R2.5$Incorrect.count)
tab1(FLR_Not_Adj_Syn7058_PSM_R1$Incorrect.count)

FLR_Adj_Syn7058_PSM<-FLR_Adj(Bin_adj_Syn7058_PSM)

FLR_Adj_Syn7058_PSM_R10 <- FLR_Adj_Syn7058_PSM[1:max(which(FLR_Adj_Syn7058_PSM$FLR_Adj_Score<=0.1)),]
FLR_Adj_Syn7058_PSM_R5 <- FLR_Adj_Syn7058_PSM[1:max(which(FLR_Adj_Syn7058_PSM$FLR_Adj_Score<=0.05)),]
FLR_Adj_Syn7058_PSM_R2.5 <- FLR_Adj_Syn7058_PSM[1:max(which(FLR_Adj_Syn7058_PSM$FLR_Adj_Score<=0.025)),]
FLR_Adj_Syn7058_PSM_R1 <- FLR_Adj_Syn7058_PSM[1:max(which(FLR_Adj_Syn7058_PSM$FLR_Adj_Score<=0.01)),]

tab1(FLR_Adj_Syn7058_PSM_R10$Incorrect.count)
tab1(FLR_Adj_Syn7058_PSM_R5$Incorrect.count)
tab1(FLR_Adj_Syn7058_PSM_R2.5$Incorrect.count)
tab1(FLR_Adj_Syn7058_PSM_R1$Incorrect.count)


FLR_Adj_Syn138_PSM_R10_WA <- FLR_Adj_Syn138_PSM_R10[FLR_Adj_Syn138_PSM_R10$Amino!="A",]
FLR_Adj_Syn138_PSM_R5_WA <- FLR_Adj_Syn138_PSM_R5[FLR_Adj_Syn138_PSM_R5$Amino!="A",]
FLR_Adj_Syn138_PSM_R2.5_WA <- FLR_Adj_Syn138_PSM_R2.5[FLR_Adj_Syn138_PSM_R2.5$Amino!="A",]
FLR_Adj_Syn138_PSM_R1_WA <- FLR_Adj_Syn138_PSM_R1[FLR_Adj_Syn138_PSM_R1$Amino!="A",]

CrossTable(FLR_Adj_Syn138_PSM_R10_WA$Incorrect.count,FLR_Adj_Syn138_PSM_R10_WA$Amino)
CrossTable(FLR_Adj_Syn138_PSM_R5_WA$Incorrect.count,FLR_Adj_Syn138_PSM_R5_WA$Amino)
CrossTable(FLR_Adj_Syn138_PSM_R2.5_WA$Incorrect.count,FLR_Adj_Syn138_PSM_R2.5_WA$Amino)
CrossTable(FLR_Adj_Syn138_PSM_R1_WA$Incorrect.count,FLR_Adj_Syn138_PSM_R1_WA$Amino)

FLR_Adj_Syn7058_PSM_R10_WA <- FLR_Adj_Syn7058_PSM_R10[FLR_Adj_Syn7058_PSM_R10$Amino!="A",]
FLR_Adj_Syn7058_PSM_R5_WA <- FLR_Adj_Syn7058_PSM_R5[FLR_Adj_Syn7058_PSM_R5$Amino!="A",]
FLR_Adj_Syn7058_PSM_R2.5_WA <- FLR_Adj_Syn7058_PSM_R2.5[FLR_Adj_Syn7058_PSM_R2.5$Amino!="A",]
FLR_Adj_Syn7058_PSM_R1_WA <- FLR_Adj_Syn7058_PSM_R1[FLR_Adj_Syn7058_PSM_R1$Amino!="A",]

CrossTable(FLR_Adj_Syn7058_PSM_R10_WA$Incorrect.count,FLR_Adj_Syn7058_PSM_R10_WA$Amino)
CrossTable(FLR_Adj_Syn7058_PSM_R5_WA$Incorrect.count,FLR_Adj_Syn7058_PSM_R5_WA$Amino)
CrossTable(FLR_Adj_Syn7058_PSM_R2.5_WA$Incorrect.count,FLR_Adj_Syn7058_PSM_R2.5_WA$Amino)
CrossTable(FLR_Adj_Syn7058_PSM_R1_WA$Incorrect.count,FLR_Adj_Syn7058_PSM_R1_WA$Amino)



FLR_Not_Adj_Syn138_PSM_R10_WA <- FLR_Not_Adj_Syn138_PSM_R10[FLR_Not_Adj_Syn138_PSM_R10$Amino!="A",]
FLR_Not_Adj_Syn138_PSM_R5_WA <- FLR_Not_Adj_Syn138_PSM_R5[FLR_Not_Adj_Syn138_PSM_R5$Amino!="A",]
FLR_Not_Adj_Syn138_PSM_R2.5_WA <- FLR_Not_Adj_Syn138_PSM_R2.5[FLR_Not_Adj_Syn138_PSM_R2.5$Amino!="A",]
FLR_Not_Adj_Syn138_PSM_R1_WA <- FLR_Not_Adj_Syn138_PSM_R1[FLR_Not_Adj_Syn138_PSM_R1$Amino!="A",]

CrossTable(FLR_Not_Adj_Syn138_PSM_R10_WA$Incorrect.count,FLR_Not_Adj_Syn138_PSM_R10_WA$Amino)
CrossTable(FLR_Not_Adj_Syn138_PSM_R5_WA$Incorrect.count,FLR_Not_Adj_Syn138_PSM_R5_WA$Amino)
CrossTable(FLR_Not_Adj_Syn138_PSM_R2.5_WA$Incorrect.count,FLR_Not_Adj_Syn138_PSM_R2.5_WA$Amino)
CrossTable(FLR_Not_Adj_Syn138_PSM_R1_WA$Incorrect.count,FLR_Not_Adj_Syn138_PSM_R1_WA$Amino)

FLR_Not_Adj_Syn7058_PSM_R10_WA <- FLR_Not_Adj_Syn7058_PSM_R10[FLR_Not_Adj_Syn7058_PSM_R10$Amino!="A",]
FLR_Not_Adj_Syn7058_PSM_R5_WA <- FLR_Not_Adj_Syn7058_PSM_R5[FLR_Not_Adj_Syn7058_PSM_R5$Amino!="A",]
FLR_Not_Adj_Syn7058_PSM_R2.5_WA <- FLR_Not_Adj_Syn7058_PSM_R2.5[FLR_Not_Adj_Syn7058_PSM_R2.5$Amino!="A",]
FLR_Not_Adj_Syn7058_PSM_R1_WA <- FLR_Not_Adj_Syn7058_PSM_R1[FLR_Not_Adj_Syn7058_PSM_R1$Amino!="A",]

CrossTable(FLR_Not_Adj_Syn7058_PSM_R10_WA$Incorrect.count,FLR_Not_Adj_Syn7058_PSM_R10_WA$Amino)
CrossTable(FLR_Not_Adj_Syn7058_PSM_R5_WA$Incorrect.count,FLR_Not_Adj_Syn7058_PSM_R5_WA$Amino)
CrossTable(FLR_Not_Adj_Syn7058_PSM_R2.5_WA$Incorrect.count,FLR_Not_Adj_Syn7058_PSM_R2.5_WA$Amino)
CrossTable(FLR_Not_Adj_Syn7058_PSM_R1_WA$Incorrect.count,FLR_Not_Adj_Syn7058_PSM_R1_WA$Amino)



