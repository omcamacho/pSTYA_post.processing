

# Table 7

source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function Peptidoform level.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/PeptidoformToPeptide_Num Phospo Complex Algorithm.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/PeptidoformToPeptide_Num Phospo Using Peptidoform Maximum.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/PeptidoformToPeptide_Num Phospo Using Peptidoform Mean.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Num protein sites and Unique peptides.R')


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

tab1(Adjusted_PXD000923$Amino)
tab1(Adjusted_PXD002222$Amino)
tab1(Adjusted_PXD002756$Amino)
tab1(Adjusted_PXD004705$Amino)
tab1(Adjusted_PXD004939$Amino)
tab1(Adjusted_PXD005241$Amino)
tab1(Adjusted_PXD012764$Amino)
tab1(Adjusted_PXD019291$Amino)

N_Prot_Adjusted_PXD000923<-Num_Proteins_Unique_Pep(Adjusted_PXD000923,Adjusted_PXD000923)
N_Prot_Adjusted_PXD002222<-Num_Proteins_Unique_Pep(Adjusted_PXD002222,Adjusted_PXD002222)
N_Prot_Adjusted_PXD002756<-Num_Proteins_Unique_Pep(Adjusted_PXD002756,Adjusted_PXD002756)
N_Prot_Adjusted_PXD004705<-Num_Proteins_Unique_Pep(Adjusted_PXD004705,Adjusted_PXD004705)
N_Prot_Adjusted_PXD004939<-Num_Proteins_Unique_Pep(Adjusted_PXD004939,Adjusted_PXD004939)
N_Prot_Adjusted_PXD005241<-Num_Proteins_Unique_Pep(Adjusted_PXD005241,Adjusted_PXD005241)
N_Prot_Adjusted_PXD012764<-Num_Proteins_Unique_Pep(Adjusted_PXD012764,Adjusted_PXD012764)
N_Prot_Adjusted_PXD019291<-Num_Proteins_Unique_Pep(Adjusted_PXD019291,Adjusted_PXD019291)

N_Prot_Adjusted_PXD000923[c(1,4)]
N_Prot_Adjusted_PXD002222[c(1,4)]
N_Prot_Adjusted_PXD002756[c(1,4)]
N_Prot_Adjusted_PXD004705[c(1,4)]
N_Prot_Adjusted_PXD004939[c(1,4)]
N_Prot_Adjusted_PXD005241[c(1,4)]
N_Prot_Adjusted_PXD012764[c(1,4)]
N_Prot_Adjusted_PXD019291[c(1,4)]

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


N_Prot_Adjusted_PXD000923_01_05<-Num_Proteins_Unique_Pep(FLR_PXD000923_Peptidoform_01,FLR_PXD000923_Peptidoform_05)
N_Prot_Adjusted_PXD002222_01_05<-Num_Proteins_Unique_Pep(FLR_PXD002222_Peptidoform_01,FLR_PXD002222_Peptidoform_05)
N_Prot_Adjusted_PXD002756_01_05<-Num_Proteins_Unique_Pep(FLR_PXD002756_Peptidoform_01,FLR_PXD002756_Peptidoform_05)
N_Prot_Adjusted_PXD004705_01_05<-Num_Proteins_Unique_Pep(FLR_PXD004705_Peptidoform_01,FLR_PXD004705_Peptidoform_05)
N_Prot_Adjusted_PXD004939_01_05<-Num_Proteins_Unique_Pep(FLR_PXD004939_Peptidoform_01,FLR_PXD004939_Peptidoform_05)
N_Prot_Adjusted_PXD005241_01_05<-Num_Proteins_Unique_Pep(FLR_PXD005241_Peptidoform_01,FLR_PXD005241_Peptidoform_05)
N_Prot_Adjusted_PXD012764_01_05<-Num_Proteins_Unique_Pep(FLR_PXD012764_Peptidoform_01,FLR_PXD012764_Peptidoform_05)
N_Prot_Adjusted_PXD019291_01_05<-Num_Proteins_Unique_Pep(FLR_PXD019291_Peptidoform_01,FLR_PXD019291_Peptidoform_05)

N_Prot_Adjusted_PXD000923_01_05[1:4]
N_Prot_Adjusted_PXD002222_01_05[1:4]
N_Prot_Adjusted_PXD002756_01_05[1:4]
N_Prot_Adjusted_PXD004705_01_05[1:4]
N_Prot_Adjusted_PXD004939_01_05[1:4]
N_Prot_Adjusted_PXD005241_01_05[1:4]
N_Prot_Adjusted_PXD012764_01_05[1:4]
N_Prot_Adjusted_PXD019291_01_05[1:4]



# Peptide level analysis #
#####################################################################

## CA ##
########

PepLevel_PXD000923_CA <- PeptidoformToPeptide_CA(Bin_Adjusted_PXD000923)
PepLevel_PXD002222_CA <- PeptidoformToPeptide_CA(Bin_Adjusted_PXD002222)
PepLevel_PXD002756_CA <- PeptidoformToPeptide_CA(Bin_Adjusted_PXD002756)
PepLevel_PXD004705_CA <- PeptidoformToPeptide_CA(Bin_Adjusted_PXD004705)
PepLevel_PXD004939_CA <- PeptidoformToPeptide_CA(Bin_Adjusted_PXD004939)
PepLevel_PXD005241_CA <- PeptidoformToPeptide_CA(Bin_Adjusted_PXD005241)
PepLevel_PXD012764_CA <- PeptidoformToPeptide_CA(Bin_Adjusted_PXD012764)
PepLevel_PXD019291_CA <- PeptidoformToPeptide_CA(Bin_Adjusted_PXD019291)

PepLevel_PXD000923_CA$Bin_Adjusted_Score <- PepLevel_PXD000923_CA$NewScore3
PepLevel_PXD002222_CA$Bin_Adjusted_Score <- PepLevel_PXD002222_CA$NewScore3
PepLevel_PXD002756_CA$Bin_Adjusted_Score <- PepLevel_PXD002756_CA$NewScore3
PepLevel_PXD004705_CA$Bin_Adjusted_Score <- PepLevel_PXD004705_CA$NewScore3
PepLevel_PXD004939_CA$Bin_Adjusted_Score <- PepLevel_PXD004939_CA$NewScore3
PepLevel_PXD005241_CA$Bin_Adjusted_Score <- PepLevel_PXD005241_CA$NewScore3
PepLevel_PXD012764_CA$Bin_Adjusted_Score <- PepLevel_PXD012764_CA$NewScore3
PepLevel_PXD019291_CA$Bin_Adjusted_Score <- PepLevel_PXD019291_CA$NewScore3


PepLevel_PXD000923_CA <- FLR_Adj(PepLevel_PXD000923_CA)
PepLevel_PXD002222_CA <- FLR_Adj(PepLevel_PXD002222_CA)
PepLevel_PXD002756_CA <- FLR_Adj(PepLevel_PXD002756_CA)
PepLevel_PXD004705_CA <- FLR_Adj(PepLevel_PXD004705_CA)
PepLevel_PXD004939_CA <- FLR_Adj(PepLevel_PXD004939_CA)
PepLevel_PXD005241_CA <- FLR_Adj(PepLevel_PXD005241_CA)
PepLevel_PXD012764_CA <- FLR_Adj(PepLevel_PXD012764_CA)
PepLevel_PXD019291_CA <- FLR_Adj(PepLevel_PXD019291_CA)

PepLevel_PXD000923_CA$dataset <- "PXD000923"
PepLevel_PXD002222_CA$dataset <- "PXD002222"
PepLevel_PXD002756_CA$dataset <- "PXD002756"
PepLevel_PXD004705_CA$dataset <- "PXD004705"
PepLevel_PXD004939_CA$dataset <- "PXD004939"
PepLevel_PXD005241_CA$dataset <- "PXD005241"
PepLevel_PXD012764_CA$dataset <- "PXD012764"
PepLevel_PXD019291_CA$dataset <- "PXD019291"

AllRice_Peptide_CA<-dplyr::bind_rows(PepLevel_PXD000923_CA, PepLevel_PXD002222_CA, PepLevel_PXD002756_CA, PepLevel_PXD004705_CA,
                                     PepLevel_PXD004939_CA, PepLevel_PXD005241_CA, PepLevel_PXD012764_CA, PepLevel_PXD019291_CA)


N_Prot_PepLevel_PXD000923_CA<-Num_Proteins_Unique_Pep(PepLevel_PXD000923_CA,PepLevel_PXD000923_CA)
N_Prot_PepLevel_PXD002222_CA<-Num_Proteins_Unique_Pep(PepLevel_PXD002222_CA,PepLevel_PXD002222_CA)
N_Prot_PepLevel_PXD002756_CA<-Num_Proteins_Unique_Pep(PepLevel_PXD002756_CA,PepLevel_PXD002756_CA)
N_Prot_PepLevel_PXD004705_CA<-Num_Proteins_Unique_Pep(PepLevel_PXD004705_CA,PepLevel_PXD004705_CA)
N_Prot_PepLevel_PXD004939_CA<-Num_Proteins_Unique_Pep(PepLevel_PXD004939_CA,PepLevel_PXD004939_CA)
N_Prot_PepLevel_PXD005241_CA<-Num_Proteins_Unique_Pep(PepLevel_PXD005241_CA,PepLevel_PXD005241_CA)
N_Prot_PepLevel_PXD012764_CA<-Num_Proteins_Unique_Pep(PepLevel_PXD012764_CA,PepLevel_PXD012764_CA)
N_Prot_PepLevel_PXD019291_CA<-Num_Proteins_Unique_Pep(PepLevel_PXD019291_CA,PepLevel_PXD019291_CA)

N_Prot_PepLevel_PXD000923_CA[c(1,4)]
N_Prot_PepLevel_PXD002222_CA[c(1,4)]
N_Prot_PepLevel_PXD002756_CA[c(1,4)]
N_Prot_PepLevel_PXD004705_CA[c(1,4)]
N_Prot_PepLevel_PXD004939_CA[c(1,4)]
N_Prot_PepLevel_PXD005241_CA[c(1,4)]
N_Prot_PepLevel_PXD012764_CA[c(1,4)]
N_Prot_PepLevel_PXD019291_CA[c(1,4)]

FLR_PepLevel_PXD000923_CA_05 <- PepLevel_PXD000923_CA[1:max(which(PepLevel_PXD000923_CA$FLR_Adj_Score<=0.05)),]
FLR_PepLevel_PXD002222_CA_05 <- PepLevel_PXD002222_CA[1:max(which(PepLevel_PXD002222_CA$FLR_Adj_Score<=0.05)),]
FLR_PepLevel_PXD002756_CA_05 <- PepLevel_PXD002756_CA[1:max(which(PepLevel_PXD002756_CA$FLR_Adj_Score<=0.05)),]
FLR_PepLevel_PXD004705_CA_05 <- PepLevel_PXD004705_CA[1:max(which(PepLevel_PXD004705_CA$FLR_Adj_Score<=0.05)),]
FLR_PepLevel_PXD004939_CA_05 <- PepLevel_PXD004939_CA[1:max(which(PepLevel_PXD004939_CA$FLR_Adj_Score<=0.05)),]
FLR_PepLevel_PXD005241_CA_05 <- PepLevel_PXD005241_CA[1:max(which(PepLevel_PXD005241_CA$FLR_Adj_Score<=0.05)),]
FLR_PepLevel_PXD012764_CA_05 <- PepLevel_PXD012764_CA[1:max(which(PepLevel_PXD012764_CA$FLR_Adj_Score<=0.05)),]
FLR_PepLevel_PXD019291_CA_05 <- PepLevel_PXD019291_CA[1:max(which(PepLevel_PXD019291_CA$FLR_Adj_Score<=0.05)),]

FLR_PepLevel_PXD000923_CA_01 <- PepLevel_PXD000923_CA[1:max(which(PepLevel_PXD000923_CA$FLR_Adj_Score<=0.01)),]
FLR_PepLevel_PXD002222_CA_01 <- PepLevel_PXD002222_CA[1:max(which(PepLevel_PXD002222_CA$FLR_Adj_Score<=0.01)),]
FLR_PepLevel_PXD002756_CA_01 <- PepLevel_PXD002756_CA[1:max(which(PepLevel_PXD002756_CA$FLR_Adj_Score<=0.01)),]
FLR_PepLevel_PXD004705_CA_01 <- PepLevel_PXD004705_CA[1:max(which(PepLevel_PXD004705_CA$FLR_Adj_Score<=0.01)),]
FLR_PepLevel_PXD004939_CA_01 <- PepLevel_PXD004939_CA[1:max(which(PepLevel_PXD004939_CA$FLR_Adj_Score<=0.01)),]
FLR_PepLevel_PXD005241_CA_01 <- PepLevel_PXD005241_CA[1:max(which(PepLevel_PXD005241_CA$FLR_Adj_Score<=0.01)),]
FLR_PepLevel_PXD012764_CA_01 <- PepLevel_PXD012764_CA[1:max(which(PepLevel_PXD012764_CA$FLR_Adj_Score<=0.01)),]
FLR_PepLevel_PXD019291_CA_01 <- PepLevel_PXD019291_CA[1:max(which(PepLevel_PXD019291_CA$FLR_Adj_Score<=0.01)),]


N_Prot_PepLevel_PXD000923_CA_01_05<-Num_Proteins_Unique_Pep(FLR_PepLevel_PXD000923_CA_01,FLR_PepLevel_PXD000923_CA_05)
N_Prot_PepLevel_PXD002222_CA_01_05<-Num_Proteins_Unique_Pep(FLR_PepLevel_PXD002222_CA_01,FLR_PepLevel_PXD002222_CA_05)
N_Prot_PepLevel_PXD002756_CA_01_05<-Num_Proteins_Unique_Pep(FLR_PepLevel_PXD002756_CA_01,FLR_PepLevel_PXD002756_CA_05)
N_Prot_PepLevel_PXD004705_CA_01_05<-Num_Proteins_Unique_Pep(FLR_PepLevel_PXD004705_CA_01,FLR_PepLevel_PXD004705_CA_05)
N_Prot_PepLevel_PXD004939_CA_01_05<-Num_Proteins_Unique_Pep(FLR_PepLevel_PXD004939_CA_01,FLR_PepLevel_PXD004939_CA_05)
N_Prot_PepLevel_PXD005241_CA_01_05<-Num_Proteins_Unique_Pep(FLR_PepLevel_PXD005241_CA_01,FLR_PepLevel_PXD005241_CA_05)
N_Prot_PepLevel_PXD012764_CA_01_05<-Num_Proteins_Unique_Pep(FLR_PepLevel_PXD012764_CA_01,FLR_PepLevel_PXD012764_CA_05)
N_Prot_PepLevel_PXD019291_CA_01_05<-Num_Proteins_Unique_Pep(FLR_PepLevel_PXD019291_CA_01,FLR_PepLevel_PXD019291_CA_05)

N_Prot_PepLevel_PXD000923_CA_01_05[1:4]
N_Prot_PepLevel_PXD002222_CA_01_05[1:4]
N_Prot_PepLevel_PXD002756_CA_01_05[1:4]
N_Prot_PepLevel_PXD004705_CA_01_05[1:4]
N_Prot_PepLevel_PXD004939_CA_01_05[1:4]
N_Prot_PepLevel_PXD005241_CA_01_05[1:4]
N_Prot_PepLevel_PXD012764_CA_01_05[1:4]
N_Prot_PepLevel_PXD019291_CA_01_05[1:4]

