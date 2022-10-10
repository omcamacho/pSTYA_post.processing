
# Figure 6

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



################## Arabidopsis ###########################
##########################################################


Ara_PXD008355 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD008355.csv')


BinAdj_Ara_PXD008355<-binAdjustPform(Ara_PXD008355)
BinAdj_Ara_PXD008355$Amino <- str_sub(substr(BinAdj_Ara_PXD008355$Peptide,1,BinAdj_Ara_PXD008355$PTM_positions),-1)

BinAdj_Ara_PXD008355$Bin_Adjusted_Score<-BinAdj_Ara_PXD008355$NewScore3

FLR_Ara_PXD008355 <- FLR_Adj(BinAdj_Ara_PXD008355)


FLR_Ara_PXD008355_R5 <- FLR_Ara_PXD008355[1:max(which(FLR_Ara_PXD008355$FLR_Adj_Score<=0.05)),]
FLR_Ara_PXD008355_R1 <- FLR_Ara_PXD008355[1:max(which(FLR_Ara_PXD008355$FLR_Adj_Score<=0.01)),]


BinAdj_Ara_PXD008355_Pep_CA <- PeptidoformToPeptide_CA(BinAdj_Ara_PXD008355)
BinAdj_Ara_PXD008355_Pep_max <- PeptidoformToPeptide_max(BinAdj_Ara_PXD008355)
BinAdj_Ara_PXD008355_Pep_mean <- PeptidoformToPeptide_mean(BinAdj_Ara_PXD008355)

names(BinAdj_Ara_PXD008355_Pep_CA)[names(BinAdj_Ara_PXD008355_Pep_CA) == "NewScore3"]<-"Bin_Adjusted_Score"
names(BinAdj_Ara_PXD008355_Pep_max)[names(BinAdj_Ara_PXD008355_Pep_max) == "NewScore3"]<-"Bin_Adjusted_Score"
names(BinAdj_Ara_PXD008355_Pep_mean)[names(BinAdj_Ara_PXD008355_Pep_mean) == "NewScore3"]<-"Bin_Adjusted_Score"

FLR_Ara_PXD008355_Pep_CA <- FLR_Adj(BinAdj_Ara_PXD008355_Pep_CA)
FLR_Ara_PXD008355_Pep_max <- FLR_Adj(BinAdj_Ara_PXD008355_Pep_max)
FLR_Ara_PXD008355_Pep_mean <- FLR_Adj(BinAdj_Ara_PXD008355_Pep_mean)


FLR_Ara_PXD008355_Pep_CA_R5 <- FLR_Ara_PXD008355_Pep_CA[1:max(which(FLR_Ara_PXD008355_Pep_CA$FLR_Adj_Score<=0.05)),]
FLR_Ara_PXD008355_Pep_max_R5 <- FLR_Ara_PXD008355_Pep_max[1:max(which(FLR_Ara_PXD008355_Pep_max$FLR_Adj_Score<=0.05)),]
FLR_Ara_PXD008355_Pep_mean_R5 <- FLR_Ara_PXD008355_Pep_mean[1:max(which(FLR_Ara_PXD008355_Pep_mean$FLR_Adj_Score<=0.05)),]
FLR_Ara_PXD008355_Pep_CA_R1 <- FLR_Ara_PXD008355_Pep_CA[1:max(which(FLR_Ara_PXD008355_Pep_CA$FLR_Adj_Score<=0.01)),]
FLR_Ara_PXD008355_Pep_max_R1 <- FLR_Ara_PXD008355_Pep_max[1:max(which(FLR_Ara_PXD008355_Pep_max$FLR_Adj_Score<=0.01)),]
FLR_Ara_PXD008355_Pep_mean_R1 <- FLR_Ara_PXD008355_Pep_mean[1:max(which(FLR_Ara_PXD008355_Pep_mean$FLR_Adj_Score<=0.01)),]


################################ Counting unique protein sites ###########################################

########### Peptidoform level ####################

Pro_Loc_Ara_PXD008355 <- FLR_Ara_PXD008355[!duplicated(BinAdj_Ara_PXD008355$PROTEIN_LOC), ]


Pro_Loc_FLR_Ara_PXD008355_R1 <- cbind.data.frame(FLR_Ara_PXD008355_R1$PROTEIN_LOC, FLR_Ara_PXD008355_R1$PROTEIN_OCC, FLR_Ara_PXD008355_R1$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_R1)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_R1 <- Pro_Loc_FLR_Ara_PXD008355_R1[!duplicated(Pro_Loc_FLR_Ara_PXD008355_R1$Protein_Location), ]


Pro_Loc_FLR_Ara_PXD008355_R5 <- cbind.data.frame(FLR_Ara_PXD008355_R5$PROTEIN_LOC, FLR_Ara_PXD008355_R5$PROTEIN_OCC, FLR_Ara_PXD008355_R5$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_R5)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_R5 <- Pro_Loc_FLR_Ara_PXD008355_R5[!duplicated(Pro_Loc_FLR_Ara_PXD008355_R5$Protein_Location), ]


############# Peptide CA approach  ###################

Pro_Loc_Ara_PXD008355_Pep_CA <- BinAdj_Ara_PXD008355_Pep_CA[!duplicated(BinAdj_Ara_PXD008355_Pep_CA$PROTEIN_LOC), ]


Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1 <- cbind.data.frame(FLR_Ara_PXD008355_Pep_CA_R1$PROTEIN_LOC,FLR_Ara_PXD008355_Pep_CA_R1$PROTEIN_OCC, FLR_Ara_PXD008355_Pep_CA_R1$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1 <- Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1[!duplicated(Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1$Protein_Location), ]


Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5 <- cbind.data.frame(FLR_Ara_PXD008355_Pep_CA_R5$PROTEIN_LOC,FLR_Ara_PXD008355_Pep_CA_R5$PROTEIN_OCC , FLR_Ara_PXD008355_Pep_CA_R5$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5 <- Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5[!duplicated(Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5$Protein_Location), ]


############# Peptide maximum approach  ###################

Pro_Loc_Ara_PXD008355_Pep_max <- BinAdj_Ara_PXD008355_Pep_max[!duplicated(BinAdj_Ara_PXD008355_Pep_max$PROTEIN_LOC), ]


Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1 <- cbind.data.frame(FLR_Ara_PXD008355_Pep_max_R1$PROTEIN_LOC,FLR_Ara_PXD008355_Pep_max_R1$PROTEIN_OCC , FLR_Ara_PXD008355_Pep_max_R1$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1 <- Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1[!duplicated(Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1$Protein_Location), ]


Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5 <- cbind.data.frame(FLR_Ara_PXD008355_Pep_max_R5$PROTEIN_LOC,FLR_Ara_PXD008355_Pep_max_R5$PROTEIN_OCC , FLR_Ara_PXD008355_Pep_max_R5$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5 <- Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5[!duplicated(Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5$Protein_Location), ]


############# Peptide Mean approach  ###################

Pro_Loc_Ara_PXD008355_Pep_mean <- BinAdj_Ara_PXD008355_Pep_mean[!duplicated(BinAdj_Ara_PXD008355_Pep_mean$PROTEIN_LOC), ]

Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1 <- cbind.data.frame(FLR_Ara_PXD008355_Pep_mean_R1$PROTEIN_LOC,FLR_Ara_PXD008355_Pep_mean_R1$PROTEIN_OCC , FLR_Ara_PXD008355_Pep_mean_R1$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1 <- Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1[!duplicated(Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1$Protein_Location), ]


Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5 <- cbind.data.frame(FLR_Ara_PXD008355_Pep_mean_R5$PROTEIN_LOC, FLR_Ara_PXD008355_Pep_mean_R5$PROTEIN_OCC, FLR_Ara_PXD008355_Pep_mean_R5$Amino)
colnames(Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5 <- Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5[!duplicated(Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5$Protein_Location), ]



################## Human ###########################
##########################################################


Hum_PXD000612<- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD000612.csv')


BinAdj_Hum_PXD000612<-binAdjustPform(Hum_PXD000612)

BinAdj_Hum_PXD000612$Bin_Adjusted_Score<-BinAdj_Hum_PXD000612$NewScore3

FLR_Hum_PXD000612 <- FLR_Adj(BinAdj_Hum_PXD000612)


FLR_Hum_PXD000612_R5 <- FLR_Hum_PXD000612[1:max(which(FLR_Hum_PXD000612$FLR_Adj_Score<=0.05)),]
FLR_Hum_PXD000612_R1 <- FLR_Hum_PXD000612[1:max(which(FLR_Hum_PXD000612$FLR_Adj_Score<=0.01)),]



# We are going to try to do the merge at the peptide-site level #
#################################################################

BinAdj_Hum_PXD000612$Amino <- str_sub(substr(BinAdj_Hum_PXD000612$Peptide,1,BinAdj_Hum_PXD000612$PTM_positions),-1)


BinAdj_Hum_PXD000612_Pep_CA <- PeptidoformToPeptide_CA(BinAdj_Hum_PXD000612)
BinAdj_Hum_PXD000612_Pep_max <- PeptidoformToPeptide_max(BinAdj_Hum_PXD000612)
BinAdj_Hum_PXD000612_Pep_mean <- PeptidoformToPeptide_mean(BinAdj_Hum_PXD000612)

names(BinAdj_Hum_PXD000612_Pep_CA)[names(BinAdj_Hum_PXD000612_Pep_CA) == "NewScore3"]<-"Bin_Adjusted_Score"
names(BinAdj_Hum_PXD000612_Pep_max)[names(BinAdj_Hum_PXD000612_Pep_max) == "NewScore3"]<-"Bin_Adjusted_Score"
names(BinAdj_Hum_PXD000612_Pep_mean)[names(BinAdj_Hum_PXD000612_Pep_mean) == "NewScore3"]<-"Bin_Adjusted_Score"


FLR_Hum_PXD000612_Pep_CA <- FLR_Adj(BinAdj_Hum_PXD000612_Pep_CA)
FLR_Hum_PXD000612_Pep_max <- FLR_Adj(BinAdj_Hum_PXD000612_Pep_max)
FLR_Hum_PXD000612_Pep_mean <- FLR_Adj(BinAdj_Hum_PXD000612_Pep_mean)


FLR_Hum_PXD000612_Pep_CA_R5 <- FLR_Hum_PXD000612_Pep_CA[1:max(which(FLR_Hum_PXD000612_Pep_CA$FLR_Adj_Score<=0.05)),]
FLR_Hum_PXD000612_Pep_max_R5 <- FLR_Hum_PXD000612_Pep_max[1:max(which(FLR_Hum_PXD000612_Pep_max$FLR_Adj_Score<=0.05)),]
FLR_Hum_PXD000612_Pep_mean_R5 <- FLR_Hum_PXD000612_Pep_mean[1:max(which(FLR_Hum_PXD000612_Pep_mean$FLR_Adj_Score<=0.05)),]
FLR_Hum_PXD000612_Pep_CA_R1 <- FLR_Hum_PXD000612_Pep_CA[1:max(which(FLR_Hum_PXD000612_Pep_CA$FLR_Adj_Score<=0.01)),]
FLR_Hum_PXD000612_Pep_max_R1 <- FLR_Hum_PXD000612_Pep_max[1:max(which(FLR_Hum_PXD000612_Pep_max$FLR_Adj_Score<=0.01)),]
FLR_Hum_PXD000612_Pep_mean_R1 <- FLR_Hum_PXD000612_Pep_mean[1:max(which(FLR_Hum_PXD000612_Pep_mean$FLR_Adj_Score<=0.01)),]



################################ Counting frequency in protein location ###########################################

########### Peptidoform level ####################
Pro_Loc_Hum_PXD000612 <- FLR_Hum_PXD000612[!duplicated(FLR_Hum_PXD000612$PROTEIN_LOC), ]

Pro_Loc_FLR_Hum_PXD000612_R1 <- cbind.data.frame(FLR_Hum_PXD000612_R1$PROTEIN_LOC, FLR_Hum_PXD000612_R1$PROTEIN_OCC, FLR_Hum_PXD000612_R1$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_R1)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_R1 <- Pro_Loc_FLR_Hum_PXD000612_R1[!duplicated(Pro_Loc_FLR_Hum_PXD000612_R1$Protein_Location), ]


Pro_Loc_FLR_Hum_PXD000612_R5 <- cbind.data.frame(FLR_Hum_PXD000612_R5$PROTEIN_LOC, FLR_Hum_PXD000612_R5$PROTEIN_OCC, FLR_Hum_PXD000612_R5$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_R5)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_R5 <- Pro_Loc_FLR_Hum_PXD000612_R5[!duplicated(Pro_Loc_FLR_Hum_PXD000612_R5$Protein_Location), ]


############# Peptide CA approach  ###################

Pro_Loc_Hum_PXD000612_CA <- BinAdj_Hum_PXD000612_Pep_CA[!duplicated(BinAdj_Hum_PXD000612_Pep_CA$PROTEIN_LOC), ]


Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1 <- cbind.data.frame(FLR_Hum_PXD000612_Pep_CA_R1$PROTEIN_LOC,FLR_Hum_PXD000612_Pep_CA_R1$PROTEIN_OCC, FLR_Hum_PXD000612_Pep_CA_R1$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1 <- Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1[!duplicated(Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1$Protein_Location), ]


Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5 <- cbind.data.frame(FLR_Hum_PXD000612_Pep_CA_R5$PROTEIN_LOC,FLR_Hum_PXD000612_Pep_CA_R5$PROTEIN_OCC , FLR_Hum_PXD000612_Pep_CA_R5$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5 <- Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5[!duplicated(Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5$Protein_Location), ]



############# Peptide maximum approach  ###################

Pro_Loc_Hum_PXD000612_max <- BinAdj_Hum_PXD000612_Pep_max[!duplicated(BinAdj_Hum_PXD000612_Pep_max$PROTEIN_LOC), ]


Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1 <- cbind.data.frame(FLR_Hum_PXD000612_Pep_max_R1$PROTEIN_LOC,FLR_Hum_PXD000612_Pep_max_R1$PROTEIN_OCC , FLR_Hum_PXD000612_Pep_max_R1$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1 <- Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1[!duplicated(Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1$Protein_Location), ]


Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5 <- cbind.data.frame(FLR_Hum_PXD000612_Pep_max_R5$PROTEIN_LOC,FLR_Hum_PXD000612_Pep_max_R5$PROTEIN_OCC , FLR_Hum_PXD000612_Pep_max_R5$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5 <- Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5[!duplicated(Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5$Protein_Location), ]


############# Peptide Mean approach  ###################

Pro_Loc_Hum_PXD000612_mean <- BinAdj_Hum_PXD000612_Pep_mean[!duplicated(BinAdj_Hum_PXD000612_Pep_mean$PROTEIN_LOC), ]

Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1 <- cbind.data.frame(FLR_Hum_PXD000612_Pep_mean_R1$PROTEIN_LOC,FLR_Hum_PXD000612_Pep_mean_R1$PROTEIN_OCC , FLR_Hum_PXD000612_Pep_mean_R1$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1 <- Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1[!duplicated(Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1$Protein_Location), ]


Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5 <- cbind.data.frame(FLR_Hum_PXD000612_Pep_mean_R5$PROTEIN_LOC, FLR_Hum_PXD000612_Pep_mean_R5$PROTEIN_OCC, FLR_Hum_PXD000612_Pep_mean_R5$Amino)
colnames(Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5)<- c("Protein_Location","Unique_Peptide","Amino")
Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5 <- Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5[!duplicated(Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5$Protein_Location), ]

# Preparing labels for Plot #
#############################


Pro_Loc_FLR_Ara_PXD008355_R1$method <- "Binomial_PformMax"
Pro_Loc_FLR_Ara_PXD008355_R5$method <- "Binomial_PformMax"

Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1$method <- "Binomial_PeptideCA"
Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5$method <- "Binomial_PeptideCA"

Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1$method <- "Binomial_PeptideMax"
Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5$method <- "Binomial_PeptideMax"

Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1$method <- "Binomial_PeptideMM"
Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5$method <- "Binomial_PeptideMM"

Pro_Loc_FLR_Ara_PXD008355_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Ara_PXD008355_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1$Species<- "Arabidopsis"
Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5$Species <- "Arabidopsis"


Pro_Loc_FLR_Ara_PXD008355_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Ara_PXD008355_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1$FLR<- "FLR <= 1%"
Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5$FLR <- "FLR <= 5%"



Pro_Loc_FLR_Hum_PXD000612_R1$method <- "Binomial_PformMax"
Pro_Loc_FLR_Hum_PXD000612_R5$method <- "Binomial_PformMax"

Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1$method <- "Binomial_PeptideCA"
Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5$method <- "Binomial_PeptideCA"

Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1$method <- "Binomial_PeptideMax"
Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5$method <- "Binomial_PeptideMax"

Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1$method <- "Binomial_PeptideMM"
Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5$method <- "Binomial_PeptideMM"

Pro_Loc_FLR_Hum_PXD000612_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Hum_PXD000612_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1$Species<- "Arabidopsis"
Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5$Species <- "Arabidopsis"

Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1$Species <- "Arabidopsis"
Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5$Species <- "Arabidopsis"


Pro_Loc_FLR_Hum_PXD000612_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Hum_PXD000612_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5$FLR <- "FLR <= 5%"

Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1$FLR <- "FLR <= 1%"
Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5$FLR <- "FLR <= 5%"


GraphData_Ara <- do.call("rbind", list(Pro_Loc_FLR_Ara_PXD008355_R1, Pro_Loc_FLR_Ara_PXD008355_R5,
                                   Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R1, Pro_Loc_FLR_Ara_PXD008355_Pep_CA_R5,
                                   Pro_Loc_FLR_Ara_PXD008355_Pep_max_R1, Pro_Loc_FLR_Ara_PXD008355_Pep_max_R5,
                                   Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R1, Pro_Loc_FLR_Ara_PXD008355_Pep_mean_R5))

GraphData_Hum <- do.call("rbind", list(Pro_Loc_FLR_Hum_PXD000612_R1, Pro_Loc_FLR_Hum_PXD000612_R5,
                                       Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R1, Pro_Loc_FLR_Hum_PXD000612_Pep_CA_R5,
                                       Pro_Loc_FLR_Hum_PXD000612_Pep_max_R1, Pro_Loc_FLR_Hum_PXD000612_Pep_max_R5,
                                       Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R1, Pro_Loc_FLR_Hum_PXD000612_Pep_mean_R5))


# PLOTS #
#########


ggplot(data = GraphData_Ara, aes(x = method, fill = FLR)) + ggtitle("Arabidopsis") +
  geom_bar(stat = "count", position = "dodge")  +
  coord_flip()
  
ggplot(data = GraphData_Hum, aes(x = method, fill = FLR)) + ggtitle("Human") +
    geom_bar(stat = "count", position = "dodge") +
  coord_flip()
