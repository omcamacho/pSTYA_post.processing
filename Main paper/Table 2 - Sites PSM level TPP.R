
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

source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level Arabidopsis.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level Syn7058.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Unadjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Prod Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Function frequency of site.R')


# Arabidopsis #
###############


Ara_PXD008355 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD008355.csv')
Ara_PXD008355$Peptidoform <- paste0(Ara_PXD008355$Peptide_mod,"_",Ara_PXD008355$PTM_positions)

BinAdj_Ara_PSM<-binAdjustAra(Ara_PXD008355)

#ProdAdj_Ara_PSM <- BinAdj_Ara_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjusted_Score = (1-prod(1-PTM_final_prob)))

Not_adj_Ara_PSM <- cbind.data.frame(BinAdj_Ara_PSM$Peptide_mod,BinAdj_Ara_PSM$Peptide,BinAdj_Ara_PSM$PTM_positions,BinAdj_Ara_PSM$PTM_final_prob)
names(Not_adj_Ara_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Bin_adj_Ara_PSM <- cbind.data.frame(BinAdj_Ara_PSM$Peptide_mod,BinAdj_Ara_PSM$Peptide,BinAdj_Ara_PSM$PTM_positions,BinAdj_Ara_PSM$NewScore3)
names(Bin_adj_Ara_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

#Prod_adj_Ara_PSM <- cbind.data.frame(ProdAdj_Ara_PSM$Peptide_mod,ProdAdj_Ara_PSM$Peptide,ProdAdj_Ara_PSM$PTM_positions,ProdAdj_Ara_PSM$Prod_Adjusted_Score)
#names(Prod_adj_Ara_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score")

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


#FLR_Prod_Adj_Ara_PSM<-FLR_Prod_Adj(Prod_adj_Ara_PSM)

#FLR_Prod_Adj_Ara_PSM_R10 <- FLR_Prod_Adj_Ara_PSM[1:max(which(FLR_Prod_Adj_Ara_PSM$FLR_Adj_Score<=0.1)),]
#FLR_Prod_Adj_Ara_PSM_R5 <- FLR_Prod_Adj_Ara_PSM[1:max(which(FLR_Prod_Adj_Ara_PSM$FLR_Adj_Score<=0.05)),]
#FLR_Prod_Adj_Ara_PSM_R2.5 <- FLR_Prod_Adj_Ara_PSM[1:max(which(FLR_Prod_Adj_Ara_PSM$FLR_Adj_Score<=0.025)),]
#FLR_Prod_Adj_Ara_PSM_R1 <- FLR_Prod_Adj_Ara_PSM[1:max(which(FLR_Prod_Adj_Ara_PSM$FLR_Adj_Score<=0.01)),]


# Human #
#########


Hum_PXD000612 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD000612.csv')
Hum_PXD000612$Peptidoform <- paste0(Hum_PXD000612$Peptide_mod,"_",Hum_PXD000612$PTM_positions)

BinAdj_Hum_PSM<-binAdjustPSM(Hum_PXD000612)

#ProdAdj_Hum_PSM <- BinAdj_Hum_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjusted_Score = (1-prod(1-PTM_final_prob)))

Not_adj_Hum_PSM <- cbind.data.frame(BinAdj_Hum_PSM$Peptide_mod,BinAdj_Hum_PSM$Peptide,BinAdj_Hum_PSM$PTM_positions,BinAdj_Hum_PSM$PTM_final_prob)
names(Not_adj_Hum_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Bin_adj_Hum_PSM <- cbind.data.frame(BinAdj_Hum_PSM$Peptide_mod,BinAdj_Hum_PSM$Peptide,BinAdj_Hum_PSM$PTM_positions,BinAdj_Hum_PSM$NewScore3)
names(Bin_adj_Hum_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

#Prod_adj_Hum_PSM <- cbind.data.frame(ProdAdj_Hum_PSM$Peptide_mod,ProdAdj_Hum_PSM$Peptide,ProdAdj_Hum_PSM$PTM_positions,ProdAdj_Hum_PSM$Prod_Adjusted_Score)
#names(Prod_adj_Hum_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score")

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


#FLR_Prod_Adj_Hum_PSM<-FLR_Prod_Adj(Prod_adj_Hum_PSM)

#FLR_Prod_Adj_Hum_PSM_R10 <- FLR_Prod_Adj_Hum_PSM[1:max(which(FLR_Prod_Adj_Hum_PSM$FLR_Adj_Score<=0.1)),]
#FLR_Prod_Adj_Hum_PSM_R5 <- FLR_Prod_Adj_Hum_PSM[1:max(which(FLR_Prod_Adj_Hum_PSM$FLR_Adj_Score<=0.05)),]
#FLR_Prod_Adj_Hum_PSM_R2.5 <- FLR_Prod_Adj_Hum_PSM[1:max(which(FLR_Prod_Adj_Hum_PSM$FLR_Adj_Score<=0.025)),]
#FLR_Prod_Adj_Hum_PSM_R1 <- FLR_Prod_Adj_Hum_PSM[1:max(which(FLR_Prod_Adj_Hum_PSM$FLR_Adj_Score<=0.01)),]



# Synthetic PXD000138 #
#######################


Syn_PXD000138 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD000138.csv')
Syn_PXD000138$Peptidoform <- paste0(Syn_PXD000138$Peptide_mod,"_",Syn_PXD000138$PTM_positions)

BinAdj_Syn138_PSM<-binAdjustSyn(Syn_PXD000138)

#ProdAdj_Syn138_PSM <- BinAdj_Syn138_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjusted_Score = (1-prod(1-PTM_final_prob)))

Not_adj_Syn138_PSM <- cbind.data.frame(BinAdj_Syn138_PSM$Peptide_mod,BinAdj_Syn138_PSM$Peptide,BinAdj_Syn138_PSM$PTM_positions,BinAdj_Syn138_PSM$Incorrect.count,BinAdj_Syn138_PSM$PTM_final_prob)
names(Not_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions", "Incorrect.count" ,"Unadjusted_Score")

Bin_adj_Syn138_PSM <- cbind.data.frame(BinAdj_Syn138_PSM$Peptide_mod,BinAdj_Syn138_PSM$Peptide,BinAdj_Syn138_PSM$PTM_positions,BinAdj_Syn138_PSM$Incorrect.count,BinAdj_Syn138_PSM$NewScore3)
names(Bin_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Incorrect.count" ,"Bin_Adjusted_Score")

#Prod_adj_Syn138_PSM <- cbind.data.frame(ProdAdj_Syn138_PSM$Peptide_mod,ProdAdj_Syn138_PSM$Peptide,ProdAdj_Syn138_PSM$PTM_positions,ProdAdj_Syn138_PSM$Incorrect.count,ProdAdj_Syn138_PSM$Prod_Adjusted_Score)
#names(Prod_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Incorrect.count" ,"Prod_Adjusted_Score")

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

#FLR_Prod_Adj_Syn138_PSM<-FLR_Prod_Adj(Prod_adj_Syn138_PSM)

#FLR_Prod_Adj_Syn138_PSM_R10 <- FLR_Prod_Adj_Syn138_PSM[1:max(which(FLR_Prod_Adj_Syn138_PSM$FLR_Adj_Score<=0.1)),]
#FLR_Prod_Adj_Syn138_PSM_R5 <- FLR_Prod_Adj_Syn138_PSM[1:max(which(FLR_Prod_Adj_Syn138_PSM$FLR_Adj_Score<=0.05)),]
#FLR_Prod_Adj_Syn138_PSM_R2.5 <- FLR_Prod_Adj_Syn138_PSM[1:max(which(FLR_Prod_Adj_Syn138_PSM$FLR_Adj_Score<=0.025)),]
#FLR_Prod_Adj_Syn138_PSM_R1 <- FLR_Prod_Adj_Syn138_PSM[1:max(which(FLR_Prod_Adj_Syn138_PSM$FLR_Adj_Score<=0.01)),]

#tab1(FLR_Prod_Adj_Syn138_PSM_R10$Incorrect.count)
#tab1(FLR_Prod_Adj_Syn138_PSM_R5$Incorrect.count)
#tab1(FLR_Prod_Adj_Syn138_PSM_R2.5$Incorrect.count)
#tab1(FLR_Prod_Adj_Syn138_PSM_R1$Incorrect.count)


# Synthetic PXD007058 #
#######################


Syn_PXD007058 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD007058.csv')
Syn_PXD007058$Peptidoform <- paste0(Syn_PXD007058$Peptide_mod,"_",Syn_PXD007058$PTM_positions)
Syn_PXD007058$All_Proteins <- Syn_PXD007058$Protein
Syn_PXD007058$Spectrum <- Syn_PXD007058$All_USI


BinAdj_Syn7058_PSM<-binAdjustSyn(Syn_PXD007058)

#ProdAdj_Syn7058_PSM <- BinAdj_Syn7058_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjusted_Score = (1-prod(1-PTM_final_prob)))

Not_adj_Syn7058_PSM <- cbind.data.frame(BinAdj_Syn7058_PSM$Peptide_mod,BinAdj_Syn7058_PSM$Peptide,BinAdj_Syn7058_PSM$PTM_positions,BinAdj_Syn7058_PSM$Synthetic_full_match_False, BinAdj_Syn7058_PSM$PTM_final_prob)
names(Not_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions", "Incorrect.count" ,"Unadjusted_Score")

Bin_adj_Syn7058_PSM <- cbind.data.frame(BinAdj_Syn7058_PSM$Peptide_mod,BinAdj_Syn7058_PSM$Peptide,BinAdj_Syn7058_PSM$PTM_positions,BinAdj_Syn7058_PSM$Synthetic_full_match_False,BinAdj_Syn7058_PSM$NewScore3)
names(Bin_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions", "Incorrect.count" ,"Bin_Adjusted_Score")

#Prod_adj_Syn7058_PSM <- cbind.data.frame(ProdAdj_Syn7058_PSM$Peptide_mod,ProdAdj_Syn7058_PSM$Peptide,ProdAdj_Syn7058_PSM$PTM_positions,ProdAdj_Syn7058_PSM$Synthetic_full_match_False,ProdAdj_Syn7058_PSM$Prod_Adjusted_Score)
#names(Prod_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions", "Incorrect.count" ,"Prod_Adjusted_Score")

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





#FLR_Prod_Adj_Syn7058_PSM<-FLR_Prod_Adj(Prod_adj_Syn7058_PSM)

#FLR_Prod_Adj_Syn7058_PSM_R10 <- FLR_Prod_Adj_Syn7058_PSM[1:max(which(FLR_Prod_Adj_Syn7058_PSM$FLR_Adj_Score<=0.1)),]
#FLR_Prod_Adj_Syn7058_PSM_R5 <- FLR_Prod_Adj_Syn7058_PSM[1:max(which(FLR_Prod_Adj_Syn7058_PSM$FLR_Adj_Score<=0.05)),]
#FLR_Prod_Adj_Syn7058_PSM_R2.5 <- FLR_Prod_Adj_Syn7058_PSM[1:max(which(FLR_Prod_Adj_Syn7058_PSM$FLR_Adj_Score<=0.025)),]
#FLR_Prod_Adj_Syn7058_PSM_R1 <- FLR_Prod_Adj_Syn7058_PSM[1:max(which(FLR_Prod_Adj_Syn7058_PSM$FLR_Adj_Score<=0.01)),]

#tab1(FLR_Prod_Adj_Syn7058_PSM_R10$Incorrect.count)
#tab1(FLR_Prod_Adj_Syn7058_PSM_R5$Incorrect.count)
#tab1(FLR_Prod_Adj_Syn7058_PSM_R2.5$Incorrect.count)
#tab1(FLR_Prod_Adj_Syn7058_PSM_R1$Incorrect.count)

