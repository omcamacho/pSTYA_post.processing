
# Figure 3. Plotting FLR vs Ordered scores #
############################################


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

FLR_Adj_Ara_PSM<-FLR_Adj(Bin_adj_Ara_PSM)


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

FLR_Adj_Hum_PSM<-FLR_Adj(Bin_adj_Hum_PSM)


# Synthetic PXD000138 #
#######################


Syn_PXD000138 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD000138.csv')

Syn_PXD000138 <- dplyr::rename(Syn_PXD000138, PTM.positions = PTM_positions)

Syn_PXD000138 <- dplyr::rename(Syn_PXD000138, PTM.Score = PTM_score)

BinAdj_Syn138_PSM<-binAdjustPSM(Syn_PXD000138)

Not_adj_Syn138_PSM <- cbind.data.frame(BinAdj_Syn138_PSM$Peptide_mod,BinAdj_Syn138_PSM$Peptide,BinAdj_Syn138_PSM$PTM.positions,BinAdj_Syn138_PSM$PTM_final_prob)
names(Not_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Unadjusted_Score")

Bin_adj_Syn138_PSM <- cbind.data.frame(BinAdj_Syn138_PSM$Peptide_mod,BinAdj_Syn138_PSM$Peptide,BinAdj_Syn138_PSM$PTM.positions,BinAdj_Syn138_PSM$Bin_Adjusted_Score)
names(Bin_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Bin_Adjusted_Score")

FLR_Not_Adj_Syn138_PSM<-FLR_NotAdj(Not_adj_Syn138_PSM)

FLR_Adj_Syn138_PSM<-FLR_Adj(Bin_adj_Syn138_PSM)


# Synthetic PXD007058 #
#######################


Syn_PXD007058 <- read.csv(file = 'D:/Paper submission/Data TPP pipeline/PXD007058.csv')
Syn_PXD007058$Peptidoform <- paste0(Syn_PXD007058$Peptide_mod,"_",Syn_PXD007058$PTM_positions,"_",Syn_PXD007058$Pool)
Syn_PXD007058$All_Proteins <- Syn_PXD007058$Protein
Syn_PXD007058$Spectrum <- Syn_PXD007058$All_USI

BinAdj_Syn7058_PSM<-binAdjustSyn(Syn_PXD007058)

Not_adj_Syn7058_PSM <- cbind.data.frame(BinAdj_Syn7058_PSM$Peptide_mod,BinAdj_Syn7058_PSM$Peptide,BinAdj_Syn7058_PSM$PTM_positions,BinAdj_Syn7058_PSM$PTM_final_prob)
names(Not_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Unadjusted_Score")

Bin_adj_Syn7058_PSM <- cbind.data.frame(BinAdj_Syn7058_PSM$Peptide_mod,BinAdj_Syn7058_PSM$Peptide,BinAdj_Syn7058_PSM$PTM_positions,BinAdj_Syn7058_PSM$Bin_Adjusted_Score)
names(Bin_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM.positions" ,"Bin_Adjusted_Score")

FLR_Not_Adj_Syn7058_PSM<-FLR_NotAdj(Not_adj_Syn7058_PSM)

FLR_Adj_Syn7058_PSM<-FLR_Adj(Bin_adj_Syn7058_PSM)



# Graphs #
##########

library(tidyr)
library(ggplot2)
library(dplyr)
library(lattice)
library(devtools)


# Plotting Arabidopsis #


Graph.data.Ara <- data.frame(cbind(FLR_Adj_Ara_PSM$Rnumber,FLR_Not_Adj_Ara_PSM$FLR_Unadjusted,FLR_Adj_Ara_PSM$FLR_Adj_Score))

colnames(Graph.data.Ara) <- c("Order_of_Scores","Unadjusted","Binomial_Adjustment")

df_Graph_Ara <- gather(Graph.data.Ara,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Ara, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Ara_PXD008355") +
  theme(legend.text=element_text(size=12), axis.title = element_text(size = 12))

# Plotting Human #


Graph.data.Hum <- data.frame(cbind(FLR_Adj_Hum_PSM$Rnumber,FLR_Not_Adj_Hum_PSM$FLR_Unadjusted,FLR_Adj_Hum_PSM$FLR_Adj_Score))

colnames(Graph.data.Hum) <- c("Order_of_Scores","Unadjusted","Binomial_Adjustment")

df_Graph_Hum <- gather(Graph.data.Hum,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Hum, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Hum_PXD000612") +
  theme(legend.text=element_text(size=12), axis.title = element_text(size = 12))

# Plotting Synthetic PXD000138 #


Graph.data.Syn138 <- data.frame(cbind(FLR_Adj_Syn138_PSM$Rnumber,FLR_Not_Adj_Syn138_PSM$FLR_Unadjusted,FLR_Adj_Syn138_PSM$FLR_Adj_Score))

colnames(Graph.data.Syn138) <- c("Order_of_Scores","Unadjusted","Binomial_Adjustment")

df_Graph_Syn138 <- gather(Graph.data.Syn138,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Syn138, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Syn_PXD000138")+
  theme(legend.text=element_text(size=12), axis.title = element_text(size = 12))

# Plotting Synthetic PXD007058 #


Graph.data.Syn7058 <- data.frame(cbind(FLR_Adj_Syn7058_PSM$Rnumber,FLR_Not_Adj_Syn7058_PSM$FLR_Unadjusted,FLR_Adj_Syn7058_PSM$FLR_Adj_Score))

colnames(Graph.data.Syn7058) <- c("Order_of_Scores","Unadjusted","Binomial_Adjustment")

df_Graph_Syn7058 <- gather(Graph.data.Syn7058,Method ,FLR, -Order_of_Scores)

qplot(data=df_Graph_Syn7058, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Syn_PXD007058") +
theme(legend.text=element_text(size=12), axis.title = element_text(size = 12))
