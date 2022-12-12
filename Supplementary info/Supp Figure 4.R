
# Supp Figure 4

library("conflicted")
suppressPackageStartupMessages(library("tidyverse"))
conflict_prefer("filter", "dplyr")
conflict_prefer("summarize", "dplyr")
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
library(gmodels)
library(epiDisplay)


source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function Peptidoform level.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/PeptidoformToPeptide_Num Phospo Complex Algorithm.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/PeptidoformToPeptide_Num Phospo Using Peptidoform Maximum.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/PeptidoformToPeptide_Num Phospo Using Peptidoform Mean.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Function frequency of site.R')



################## Arabidopsis ###########################
##########################################################



Ara_PXD008355_pSTYA <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD008355.csv')


BinAdj_Ara_PXD008355_pSTYA<-binAdjustPform(Ara_PXD008355_pSTYA)

BinAdj_Ara_PXD008355_pSTYA$Bin_Adjusted_Score <- BinAdj_Ara_PXD008355_pSTYA$NewScore3

FLR_Ara_PXD008355_pSTYA <- FLR_Adj(BinAdj_Ara_PXD008355_pSTYA)


FLR_Ara_PXD008355_pSTYA_R5 <- FLR_Ara_PXD008355_pSTYA[1:max(which(FLR_Ara_PXD008355_pSTYA$FLR_Adj_Score<=0.05)),]
FLR_Ara_PXD008355_pSTYA_R1 <- FLR_Ara_PXD008355_pSTYA[1:max(which(FLR_Ara_PXD008355_pSTYA$FLR_Adj_Score<=0.01)),]


# Applying the Max, Mean and MM for collapsing at Peptide-sequence level

BinAdj_Ara_PXD008355_pSTYA$Amino <- str_sub(substr(BinAdj_Ara_PXD008355_pSTYA$Peptide,1,BinAdj_Ara_PXD008355_pSTYA$PTM_positions),-1)
BinAdj_Ara_PXD008355_pSTYA$Amino <- str_sub(substr(BinAdj_Ara_PXD008355_pSTYA$Peptide,1,BinAdj_Ara_PXD008355_pSTYA$PTM_positions),-1)

BinAdj_Ara_PXD008355_pSTYA_Pep_CA <- PeptidoformToPeptide_CA(BinAdj_Ara_PXD008355_pSTYA)
BinAdj_Ara_PXD008355_pSTYA_Pep_max <- PeptidoformToPeptide_max(BinAdj_Ara_PXD008355_pSTYA)
BinAdj_Ara_PXD008355_pSTYA_Pep_mean <- PeptidoformToPeptide_mean(BinAdj_Ara_PXD008355_pSTYA)

BinAdj_Ara_PXD008355_pSTYA_Pep_CA$Bin_Adjusted_Score <- BinAdj_Ara_PXD008355_pSTYA_Pep_CA$NewScore3
BinAdj_Ara_PXD008355_pSTYA_Pep_max$Bin_Adjusted_Score <- BinAdj_Ara_PXD008355_pSTYA_Pep_max$NewScore3
BinAdj_Ara_PXD008355_pSTYA_Pep_mean$Bin_Adjusted_Score <- BinAdj_Ara_PXD008355_pSTYA_Pep_mean$NewScore3


FLR_Ara_PXD008355_pSTYA_Pep_CA <- FLR_Adj(BinAdj_Ara_PXD008355_pSTYA_Pep_CA)
FLR_Ara_PXD008355_pSTYA_Pep_max <- FLR_Adj(BinAdj_Ara_PXD008355_pSTYA_Pep_max)
FLR_Ara_PXD008355_pSTYA_Pep_mean <- FLR_Adj(BinAdj_Ara_PXD008355_pSTYA_Pep_mean)


FLR_Ara_PXD008355_pSTYA_Pep_CA_R5 <- FLR_Ara_PXD008355_pSTYA_Pep_CA[1:max(which(FLR_Ara_PXD008355_pSTYA_Pep_CA$FLR_Adj_Score<=0.05)),]
FLR_Ara_PXD008355_pSTYA_Pep_max_R5 <- FLR_Ara_PXD008355_pSTYA_Pep_max[1:max(which(FLR_Ara_PXD008355_pSTYA_Pep_max$FLR_Adj_Score<=0.05)),]
FLR_Ara_PXD008355_pSTYA_Pep_mean_R5 <- FLR_Ara_PXD008355_pSTYA_Pep_mean[1:max(which(FLR_Ara_PXD008355_pSTYA_Pep_mean$FLR_Adj_Score<=0.05)),]
FLR_Ara_PXD008355_pSTYA_Pep_CA_R1 <- FLR_Ara_PXD008355_pSTYA_Pep_CA[1:max(which(FLR_Ara_PXD008355_pSTYA_Pep_CA$FLR_Adj_Score<=0.01)),]
FLR_Ara_PXD008355_pSTYA_Pep_max_R1 <- FLR_Ara_PXD008355_pSTYA_Pep_max[1:max(which(FLR_Ara_PXD008355_pSTYA_Pep_max$FLR_Adj_Score<=0.01)),]
FLR_Ara_PXD008355_pSTYA_Pep_mean_R1 <- FLR_Ara_PXD008355_pSTYA_Pep_mean[1:max(which(FLR_Ara_PXD008355_pSTYA_Pep_mean$FLR_Adj_Score<=0.01)),]


################################ Unique PSM Matches ###########################################

########### Peptidoform level ####################

Pro_Loc_Ara_PXD008355_pSTYA <- FLR_Ara_PXD008355_pSTYA[!duplicated(FLR_Ara_PXD008355_pSTYA$PROTEIN_LOC), ]
tab1(Pro_Loc_Ara_PXD008355_pSTYA$Amino, sort.group = "decreasing", cum.percent = TRUE)

Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1 <- cbind.data.frame(FLR_Ara_PXD008355_pSTYA_R1$PROTEIN_LOC, FLR_Ara_PXD008355_pSTYA_R1$Peptide, FLR_Ara_PXD008355_pSTYA_R1$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1 <- Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1[!duplicated(Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1$Protein_Location), ]
tab1(Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5 <- cbind.data.frame(FLR_Ara_PXD008355_pSTYA_R5$PROTEIN_LOC, FLR_Ara_PXD008355_pSTYA_R5$Peptide, FLR_Ara_PXD008355_pSTYA_R5$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5 <- Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5[!duplicated(Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5$Protein_Location), ]
tab1(Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5$Amino, sort.group = "decreasing", cum.percent = TRUE)



############# Peptide combined algorithm approach  ###################

Pro_Loc_Ara_PXD008355_pSTYA_Pep_CA <- BinAdj_Ara_PXD008355_pSTYA_Pep_CA[!duplicated(BinAdj_Ara_PXD008355_pSTYA_Pep_CA$PROTEIN_LOC), ]
tab1(Pro_Loc_Ara_PXD008355_pSTYA_Pep_CA$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1 <- cbind.data.frame(FLR_Ara_PXD008355_pSTYA_Pep_CA_R1$PROTEIN_LOC,FLR_Ara_PXD008355_pSTYA_Pep_CA_R1$Peptide, FLR_Ara_PXD008355_pSTYA_Pep_CA_R1$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1 <- Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1[!duplicated(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1$Protein_Location), ]
tab1(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5 <- cbind.data.frame(FLR_Ara_PXD008355_pSTYA_Pep_CA_R5$PROTEIN_LOC,FLR_Ara_PXD008355_pSTYA_Pep_CA_R5$Peptide , FLR_Ara_PXD008355_pSTYA_Pep_CA_R5$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5 <- Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5[!duplicated(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5$Protein_Location), ]
tab1(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5$Amino, sort.group = "decreasing", cum.percent = TRUE)



############# Peptide maximum approach  ###################

Pro_Loc_Ara_PXD008355_pSTYA_Pep_max <- BinAdj_Ara_PXD008355_pSTYA_Pep_max[!duplicated(BinAdj_Ara_PXD008355_pSTYA_Pep_max$PROTEIN_LOC), ]
tab1(Pro_Loc_Ara_PXD008355_pSTYA_Pep_max$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1 <- cbind.data.frame(FLR_Ara_PXD008355_pSTYA_Pep_max_R1$PROTEIN_LOC,FLR_Ara_PXD008355_pSTYA_Pep_max_R1$Peptide , FLR_Ara_PXD008355_pSTYA_Pep_max_R1$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1 <- Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1[!duplicated(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1$Protein_Location), ]
tab1(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5 <- cbind.data.frame(FLR_Ara_PXD008355_pSTYA_Pep_max_R5$PROTEIN_LOC,FLR_Ara_PXD008355_pSTYA_Pep_max_R5$Peptide , FLR_Ara_PXD008355_pSTYA_Pep_max_R5$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5 <- Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5[!duplicated(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5$Protein_Location), ]
tab1(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5$Amino, sort.group = "decreasing", cum.percent = TRUE)


############# Peptide Mean approach  ###################

Pro_Loc_Ara_PXD008355_pSTYA_Pep_mean <- BinAdj_Ara_PXD008355_pSTYA_Pep_mean[!duplicated(BinAdj_Ara_PXD008355_pSTYA_Pep_mean$PROTEIN_LOC), ]
tab1(Pro_Loc_Ara_PXD008355_pSTYA_Pep_mean$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1 <- cbind.data.frame(FLR_Ara_PXD008355_pSTYA_Pep_mean_R1$PROTEIN_LOC,FLR_Ara_PXD008355_pSTYA_Pep_mean_R1$Peptide , FLR_Ara_PXD008355_pSTYA_Pep_mean_R1$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1 <- Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1[!duplicated(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1$Protein_Location), ]
tab1(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5 <- cbind.data.frame(FLR_Ara_PXD008355_pSTYA_Pep_mean_R5$PROTEIN_LOC, FLR_Ara_PXD008355_pSTYA_Pep_mean_R5$Peptide, FLR_Ara_PXD008355_pSTYA_Pep_mean_R5$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5 <- Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5[!duplicated(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5$Protein_Location), ]
tab1(Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5$Amino, sort.group = "decreasing", cum.percent = TRUE)




################## Human ###########################
##########################################################


Hum_PXD000612_pSTYA<- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline//PD_PXD000612.csv')


BinAdj_Hum_PXD000612_pSTYA<-binAdjustPform(Hum_PXD000612_pSTYA)

BinAdj_Hum_PXD000612_pSTYA$Bin_Adjusted_Score <- BinAdj_Hum_PXD000612_pSTYA$NewScore3

FLR_Hum_PXD000612_pSTYA <- FLR_Adj(BinAdj_Hum_PXD000612_pSTYA)


FLR_Hum_PXD000612_pSTYA_R5 <- FLR_Hum_PXD000612_pSTYA[1:max(which(FLR_Hum_PXD000612_pSTYA$FLR_Adj_Score<=0.05)),]
FLR_Hum_PXD000612_pSTYA_R1 <- FLR_Hum_PXD000612_pSTYA[1:max(which(FLR_Hum_PXD000612_pSTYA$FLR_Adj_Score<=0.01)),]




# We are going to try to do the merge at the peptide-site level #
#################################################################

BinAdj_Hum_PXD000612_pSTYA$Amino <- str_sub(substr(BinAdj_Hum_PXD000612_pSTYA$Peptide,1,BinAdj_Hum_PXD000612_pSTYA$PTM_positions),-1)
BinAdj_Hum_PXD000612_pSTYA$Amino <- str_sub(substr(BinAdj_Hum_PXD000612_pSTYA$Peptide,1,BinAdj_Hum_PXD000612_pSTYA$PTM_positions),-1)

BinAdj_Hum_PXD000612_pSTYA_Pep_CA <- PeptidoformToPeptide_CA(BinAdj_Hum_PXD000612_pSTYA)
BinAdj_Hum_PXD000612_pSTYA_Pep_max <- PeptidoformToPeptide_max(BinAdj_Hum_PXD000612_pSTYA)
BinAdj_Hum_PXD000612_pSTYA_Pep_mean <- PeptidoformToPeptide_mean(BinAdj_Hum_PXD000612_pSTYA)

BinAdj_Hum_PXD000612_pSTYA_Pep_CA$Bin_Adjusted_Score <- BinAdj_Hum_PXD000612_pSTYA_Pep_CA$NewScore3
BinAdj_Hum_PXD000612_pSTYA_Pep_max$Bin_Adjusted_Score <- BinAdj_Hum_PXD000612_pSTYA_Pep_max$NewScore3
BinAdj_Hum_PXD000612_pSTYA_Pep_mean$Bin_Adjusted_Score <- BinAdj_Hum_PXD000612_pSTYA_Pep_mean$NewScore3

FLR_Hum_PXD000612_pSTYA_Pep_CA <- FLR_Adj(BinAdj_Hum_PXD000612_pSTYA_Pep_CA)
FLR_Hum_PXD000612_pSTYA_Pep_max <- FLR_Adj(BinAdj_Hum_PXD000612_pSTYA_Pep_max)
FLR_Hum_PXD000612_pSTYA_Pep_mean <- FLR_Adj(BinAdj_Hum_PXD000612_pSTYA_Pep_mean)




FLR_Hum_PXD000612_pSTYA_Pep_CA_R5 <- FLR_Hum_PXD000612_pSTYA_Pep_CA[1:max(which(FLR_Hum_PXD000612_pSTYA_Pep_CA$FLR_Adj_Score<=0.05)),]
FLR_Hum_PXD000612_pSTYA_Pep_max_R5 <- FLR_Hum_PXD000612_pSTYA_Pep_max[1:max(which(FLR_Hum_PXD000612_pSTYA_Pep_max$FLR_Adj_Score<=0.05)),]
FLR_Hum_PXD000612_pSTYA_Pep_mean_R5 <- FLR_Hum_PXD000612_pSTYA_Pep_mean[1:max(which(FLR_Hum_PXD000612_pSTYA_Pep_mean$FLR_Adj_Score<=0.05)),]
FLR_Hum_PXD000612_pSTYA_Pep_CA_R1 <- FLR_Hum_PXD000612_pSTYA_Pep_CA[1:max(which(FLR_Hum_PXD000612_pSTYA_Pep_CA$FLR_Adj_Score<=0.01)),]
FLR_Hum_PXD000612_pSTYA_Pep_max_R1 <- FLR_Hum_PXD000612_pSTYA_Pep_max[1:max(which(FLR_Hum_PXD000612_pSTYA_Pep_max$FLR_Adj_Score<=0.01)),]
FLR_Hum_PXD000612_pSTYA_Pep_mean_R1 <- FLR_Hum_PXD000612_pSTYA_Pep_mean[1:max(which(FLR_Hum_PXD000612_pSTYA_Pep_mean$FLR_Adj_Score<=0.01)),]



################################ Unique PSM Matches ###########################################

########### Peptidoform level ####################
Pro_Loc_Hum_PXD000612_pSTYA <- FLR_Hum_PXD000612_pSTYA[!duplicated(FLR_Hum_PXD000612_pSTYA$PROTEIN_LOC), ]
tab1(Pro_Loc_Hum_PXD000612_pSTYA$Amino, sort.group = "decreasing", cum.percent = TRUE)

Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1 <- cbind.data.frame(FLR_Hum_PXD000612_pSTYA_R1$PROTEIN_LOC, FLR_Hum_PXD000612_pSTYA_R1$Peptide, FLR_Hum_PXD000612_pSTYA_R1$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1 <- Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1[!duplicated(Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1$Protein_Location), ]
tab1(Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5 <- cbind.data.frame(FLR_Hum_PXD000612_pSTYA_R5$PROTEIN_LOC, FLR_Hum_PXD000612_pSTYA_R5$Peptide, FLR_Hum_PXD000612_pSTYA_R5$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5 <- Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5[!duplicated(Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5$Protein_Location), ]
tab1(Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5$Amino, sort.group = "decreasing", cum.percent = TRUE)


############# Peptide combined algorithm approach  ###################
Pro_Loc_Hum_PXD000612_pSTYA_CA <- BinAdj_Hum_PXD000612_pSTYA_Pep_CA[!duplicated(BinAdj_Hum_PXD000612_pSTYA_Pep_CA$PROTEIN_LOC), ]
tab1(Pro_Loc_Hum_PXD000612_pSTYA_CA$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1 <- cbind.data.frame(FLR_Hum_PXD000612_pSTYA_Pep_CA_R1$PROTEIN_LOC,FLR_Hum_PXD000612_pSTYA_Pep_CA_R1$Peptide, FLR_Hum_PXD000612_pSTYA_Pep_CA_R1$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1 <- Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1[!duplicated(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1$Protein_Location), ]
tab1(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5 <- cbind.data.frame(FLR_Hum_PXD000612_pSTYA_Pep_CA_R5$PROTEIN_LOC,FLR_Hum_PXD000612_pSTYA_Pep_CA_R5$Peptide , FLR_Hum_PXD000612_pSTYA_Pep_CA_R5$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5 <- Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5[!duplicated(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5$Protein_Location), ]
tab1(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5$Amino, sort.group = "decreasing", cum.percent = TRUE)



############# Peptide maximum approach  ###################
Pro_Loc_Hum_PXD000612_pSTYA_max <- BinAdj_Hum_PXD000612_pSTYA_Pep_max[!duplicated(BinAdj_Hum_PXD000612_pSTYA_Pep_max$PROTEIN_LOC), ]
tab1(Pro_Loc_Hum_PXD000612_pSTYA_max$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1 <- cbind.data.frame(FLR_Hum_PXD000612_pSTYA_Pep_max_R1$PROTEIN_LOC,FLR_Hum_PXD000612_pSTYA_Pep_max_R1$Peptide , FLR_Hum_PXD000612_pSTYA_Pep_max_R1$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1 <- Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1[!duplicated(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1$Protein_Location), ]
tab1(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5 <- cbind.data.frame(FLR_Hum_PXD000612_pSTYA_Pep_max_R5$PROTEIN_LOC,FLR_Hum_PXD000612_pSTYA_Pep_max_R5$Peptide , FLR_Hum_PXD000612_pSTYA_Pep_max_R5$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5 <- Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5[!duplicated(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5$Protein_Location), ]
tab1(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5$Amino, sort.group = "decreasing", cum.percent = TRUE)


############# Peptide Mean approach  ###################
Pro_Loc_Hum_PXD000612_pSTYA_mean <- BinAdj_Hum_PXD000612_pSTYA_Pep_mean[!duplicated(BinAdj_Hum_PXD000612_pSTYA_Pep_mean$PROTEIN_LOC), ]
tab1(Pro_Loc_Hum_PXD000612_pSTYA_mean$Amino, sort.group = "decreasing", cum.percent = TRUE)

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1 <- cbind.data.frame(FLR_Hum_PXD000612_pSTYA_Pep_mean_R1$PROTEIN_LOC,FLR_Hum_PXD000612_pSTYA_Pep_mean_R1$Peptide , FLR_Hum_PXD000612_pSTYA_Pep_mean_R1$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1 <- Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1[!duplicated(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1$Protein_Location), ]
tab1(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1$Amino, sort.group = "decreasing", cum.percent = TRUE)


Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5 <- cbind.data.frame(FLR_Hum_PXD000612_pSTYA_Pep_mean_R5$PROTEIN_LOC, FLR_Hum_PXD000612_pSTYA_Pep_mean_R5$Peptide, FLR_Hum_PXD000612_pSTYA_Pep_mean_R5$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5)<- c("Protein_Location","Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5 <- Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5[!duplicated(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5$Protein_Location), ]
tab1(Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5$Amino, sort.group = "decreasing", cum.percent = TRUE)


# Graph #
#########

Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1$method <- "Binomial_PformMax"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5$method <- "Binomial_PformMax"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1$method <- "Binomial_PeptideCA"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5$method <- "Binomial_PeptideCA"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1$method <- "Binomial_PeptideMax"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5$method <- "Binomial_PeptideMax"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1$method <- "Binomial_PeptideMM"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5$method <- "Binomial_PeptideMM"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1$Species<- "Arabidopsis"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5$Species <- "Arabidopsis"


Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1$FLR<- "FLR <= 1%"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5$FLR <- "FLR <= 5%"





Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1$method <- "Binomial_PformMax"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5$method <- "Binomial_PformMax"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1$method <- "Binomial_PeptideCA"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5$method <- "Binomial_PeptideCA"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1$method <- "Binomial_PeptideMax"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5$method <- "Binomial_PeptideMax"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1$method <- "Binomial_PeptideMM"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5$method <- "Binomial_PeptideMM"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1$Species<- "Arabidopsis"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5$Species <- "Arabidopsis"


Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5$FLR <- "FLR <= 5%"


GraphData_Ara <- do.call("rbind", list(Pro_Loc_FLR_Ara_PXD008355_pSTYA_R1, Pro_Loc_FLR_Ara_PXD008355_pSTYA_R5,
                                   Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R1, Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_CA_R5,
                                   Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R1, Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_max_R5,
                                   Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R1, Pro_Loc_FLR_Ara_PXD008355_pSTYA_Pep_mean_R5))

GraphData_Hum <- do.call("rbind", list(Pro_Loc_FLR_Hum_PXD000612_pSTYA_R1, Pro_Loc_FLR_Hum_PXD000612_pSTYA_R5,
                                       Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R1, Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_CA_R5,
                                       Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R1, Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_max_R5,
                                       Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R1, Pro_Loc_FLR_Hum_PXD000612_pSTYA_Pep_mean_R5))


ggplot(data = GraphData_Ara, aes(x = method, fill = FLR)) + ggtitle("Arabidopsis") +
  geom_bar(stat = "count", position = "dodge")  +
  coord_flip()
  
ggplot(data = GraphData_Hum, aes(x = method, fill = FLR)) + ggtitle("Human") +
    geom_bar(stat = "count", position = "dodge") +
  coord_flip()
