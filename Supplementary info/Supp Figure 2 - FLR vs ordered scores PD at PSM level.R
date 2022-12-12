

# Supp Figure 2

# Synthetic dataset #
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

source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level Arabidopsis.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level Syn7058.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Unadjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Prod Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Function frequency of site.R')

# Arabidopsis #
###############


Ara_PXD008355 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD008355.csv')
Ara_PXD008355$Peptidoform <- paste0(Ara_PXD008355$Peptide_mod,"_",Ara_PXD008355$PTM_positions)

BinAdj_Ara_PSM<-binAdjustAra(Ara_PXD008355)

#ProdAdj_Ara_PSM <- BinAdj_Ara_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjasted_Score = (1-prod(1-PTM_final_prob)))

Not_adj_Ara_PSM <- cbind.data.frame(BinAdj_Ara_PSM$Peptide_mod,BinAdj_Ara_PSM$Peptide,BinAdj_Ara_PSM$PTM_positions,BinAdj_Ara_PSM$PTM_final_prob)
names(Not_adj_Ara_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Bin_adj_Ara_PSM <- cbind.data.frame(BinAdj_Ara_PSM$Peptide_mod,BinAdj_Ara_PSM$Peptide,BinAdj_Ara_PSM$PTM_positions,BinAdj_Ara_PSM$NewScore3)
names(Bin_adj_Ara_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

#Prod_adj_Ara_PSM <- cbind.data.frame(ProdAdj_Ara_PSM$Peptide_mod.x,ProdAdj_Ara_PSM$Peptide,ProdAdj_Ara_PSM$PTM_positions,ProdAdj_Ara_PSM$Prod_Adjasted_Score)
#names(Prod_adj_Ara_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score")

FLR_Not_Adj_Ara_PSM<-FLR_NotAdj(Not_adj_Ara_PSM)

FLR_Adj_Ara_PSM<-FLR_Adj(Bin_adj_Ara_PSM)

#FLR_Prod_Adj_Ara_PSM<-FLR_Prod_Adj(Prod_adj_Ara_PSM)



# Human #
#########


Hum_PXD000612 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD000612.csv')
Hum_PXD000612$Peptidoform <- paste0(Hum_PXD000612$Peptide_mod,"_",Hum_PXD000612$PTM_positions)

BinAdj_Hum_PSM<-binAdjustPSM(Hum_PXD000612)

#ProdAdj_Hum_PSM <- BinAdj_Hum_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjasted_Score = (1-prod(1-PTM_final_prob)))

Not_adj_Hum_PSM <- cbind.data.frame(BinAdj_Hum_PSM$Peptide_mod,BinAdj_Hum_PSM$Peptide,BinAdj_Hum_PSM$PTM_positions,BinAdj_Hum_PSM$PTM_final_prob)
names(Not_adj_Hum_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Bin_adj_Hum_PSM <- cbind.data.frame(BinAdj_Hum_PSM$Peptide_mod,BinAdj_Hum_PSM$Peptide,BinAdj_Hum_PSM$PTM_positions,BinAdj_Hum_PSM$NewScore3)
names(Bin_adj_Hum_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

#Prod_adj_Hum_PSM <- cbind.data.frame(ProdAdj_Hum_PSM$Peptide_mod.x,ProdAdj_Hum_PSM$Peptide,ProdAdj_Hum_PSM$PTM_positions,ProdAdj_Hum_PSM$Prod_Adjasted_Score)
#names(Prod_adj_Hum_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score")

FLR_Not_Adj_Hum_PSM<-FLR_NotAdj(Not_adj_Hum_PSM)

FLR_Adj_Hum_PSM<-FLR_Adj(Bin_adj_Hum_PSM)

#FLR_Prod_Adj_Hum_PSM<-FLR_Prod_Adj(Prod_adj_Hum_PSM)



# Synthetic PXD000138 #
#######################


Syn_PXD000138 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD000138.csv')
Syn_PXD000138$Peptidoform <- paste0(Syn_PXD000138$Peptide_mod,"_",Syn_PXD000138$PTM_positions)

BinAdj_Syn138_PSM<-binAdjustSyn(Syn_PXD000138)

#ProdAdj_Syn138_PSM <- BinAdj_Syn138_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjasted_Score = (1-prod(1-PTM_final_prob)))

Not_adj_Syn138_PSM <- cbind.data.frame(BinAdj_Syn138_PSM$Peptide_mod,BinAdj_Syn138_PSM$Peptide,BinAdj_Syn138_PSM$PTM_positions,BinAdj_Syn138_PSM$PTM_final_prob)
names(Not_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Bin_adj_Syn138_PSM <- cbind.data.frame(BinAdj_Syn138_PSM$Peptide_mod,BinAdj_Syn138_PSM$Peptide,BinAdj_Syn138_PSM$PTM_positions,BinAdj_Syn138_PSM$NewScore3)
names(Bin_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

#Prod_adj_Syn138_PSM <- cbind.data.frame(ProdAdj_Syn138_PSM$Peptide_mod.x,ProdAdj_Syn138_PSM$Peptide,ProdAdj_Syn138_PSM$PTM_positions,ProdAdj_Syn138_PSM$Prod_Adjasted_Score)
#names(Prod_adj_Syn138_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score")

FLR_Not_Adj_Syn138_PSM<-FLR_NotAdj(Not_adj_Syn138_PSM)

FLR_Adj_Syn138_PSM<-FLR_Adj(Bin_adj_Syn138_PSM)

#FLR_Prod_Adj_Syn138_PSM<-FLR_Prod_Adj(Prod_adj_Syn138_PSM)



# Synthetic PXD007058 #
#######################


Syn_PXD007058 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD007058.csv')
Syn_PXD007058$Peptidoform <- paste0(Syn_PXD007058$Peptide_mod,"_",Syn_PXD007058$PTM_positions)
Syn_PXD007058$PTM_final_prob <- Syn_PXD007058$PTM_score * Syn_PXD007058$Score
Syn_PXD007058$All_Proteins <- Syn_PXD007058$Protein
Syn_PXD007058$Spectrum <- Syn_PXD007058$All_USI
BinAdj_Syn7058_PSM<-binAdjustSyn(Syn_PXD007058)

#ProdAdj_Syn7058_PSM <- BinAdj_Syn7058_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjasted_Score = (1-prod(1-PTM_final_prob)))

Not_adj_Syn7058_PSM <- cbind.data.frame(BinAdj_Syn7058_PSM$Peptide_mod,BinAdj_Syn7058_PSM$Peptide,BinAdj_Syn7058_PSM$PTM_positions,BinAdj_Syn7058_PSM$PTM_final_prob)
names(Not_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Bin_adj_Syn7058_PSM <- cbind.data.frame(BinAdj_Syn7058_PSM$Peptide_mod,BinAdj_Syn7058_PSM$Peptide,BinAdj_Syn7058_PSM$PTM_positions,BinAdj_Syn7058_PSM$NewScore3)
names(Bin_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

#Prod_adj_Syn7058_PSM <- cbind.data.frame(ProdAdj_Syn7058_PSM$Peptide_mod.x,ProdAdj_Syn7058_PSM$Peptide,ProdAdj_Syn7058_PSM$PTM_positions,ProdAdj_Syn7058_PSM$Prod_Adjasted_Score)
#names(Prod_adj_Syn7058_PSM)<-c("Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score")

FLR_Not_Adj_Syn7058_PSM<-FLR_NotAdj(Not_adj_Syn7058_PSM)

FLR_Adj_Syn7058_PSM<-FLR_Adj(Bin_adj_Syn7058_PSM)

#FLR_Prod_Adj_Syn7058_PSM<-FLR_Prod_Adj(Prod_adj_Syn7058_PSM)



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



qplot(data=df_Graph_Ara, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Ara_PXD008355")


# Plotting Human #


Graph.data.Hum <- data.frame(cbind(FLR_Adj_Hum_PSM$Rnumber,FLR_Not_Adj_Hum_PSM$FLR_Unadjusted,FLR_Adj_Hum_PSM$FLR_Adj_Score))

colnames(Graph.data.Hum) <- c("Order_of_Scores","Unadjusted","Binomial_Adjustment")

df_Graph_Hum <- gather(Graph.data.Hum,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Hum, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Hum_PXD000612")


# Plotting Synthetic PXD000138 #


Graph.data.Syn138 <- data.frame(cbind(FLR_Adj_Syn138_PSM$Rnumber,FLR_Not_Adj_Syn138_PSM$FLR_Unadjusted,FLR_Adj_Syn138_PSM$FLR_Adj_Score))

colnames(Graph.data.Syn138) <- c("Order_of_Scores","Unadjusted","Binomial_Adjustment")

df_Graph_Syn138 <- gather(Graph.data.Syn138,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Syn138, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Syn_PXD000138")


# Plotting Synthetic PXD007058 #


Graph.data.Syn7058 <- data.frame(cbind(FLR_Adj_Syn7058_PSM$Rnumber,FLR_Not_Adj_Syn7058_PSM$FLR_Unadjusted,FLR_Adj_Syn7058_PSM$FLR_Adj_Score))

colnames(Graph.data.Syn7058) <- c("Order_of_Scores","Unadjusted","Binomial_Adjustment")

df_Graph_Syn7058 <- gather(Graph.data.Syn7058,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Syn7058, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Syn_PXD007058")