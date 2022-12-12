
library("conflicted")
suppressPackageStartupMessages(library("tidyverse"))
conflict_prefer("filter", "dplyr")
#> [conflicted] Will prefer dplyr::filter over any other package

library("kableExtra")
library(plyr)
library(dplyr)
library(stringr)
library(useful)
library("data.table")
conflict_prefer("mutate", "dplyr")
library(epiDisplay)

source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level Arabidopsis.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function PSM level Syn7058.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Unadjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Bin Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/FLR function Prod Adjusted.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Function frequency of site.R')

Ara_PXD008355 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD008355.csv')
Hum_PXD000612 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD000612.csv')
Syn_PXD000138 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD000138.csv')
Syn_PXD007058 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data PD pipeline/PD_PXD007058.csv')



# Arabidopsis #
###############


Ara_PXD008355$Peptidoform <- paste0(Ara_PXD008355$Peptide_mod,"_",Ara_PXD008355$PTM_positions)

BinAdj_Ara_PSM<-binAdjustAra(Ara_PXD008355)

ProdAdj_Ara_PSM <- BinAdj_Ara_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjusted_Score = (1-prod(1-PTM_final_prob)))

ProdAdj_Ara_Pform <- ProdAdj_Ara_PSM[!duplicated(ProdAdj_Ara_PSM$Peptidoform), ]
# Max 

BinAdj_Ara_Pform_Max <- aggregate(NewScore3 ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Ara_PSM, max)

total_BinAdj_Ara_Pform_Max <- merge(BinAdj_Ara_Pform_Max,Ara_PXD008355,by="Peptidoform")

BinAdj_Ara_Pform_Max <- total_BinAdj_Ara_Pform_Max %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean

BinAdj_Ara_Pform_Mean <- aggregate(NewScore3 ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Ara_PSM, mean)

total_BinAdj_Ara_Pform_Mean <- merge(BinAdj_Ara_Pform_Mean,Ara_PXD008355,by="Peptidoform")

BinAdj_Ara_Pform_Mean <- total_BinAdj_Ara_Pform_Mean %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean of the PSM first

Ara_PSM_Mean <- aggregate(PTM_final_prob ~ All_USI , data = Ara_PXD008355, mean)

names(Ara_PSM_Mean) <- c("All_USI", "PTM_final_prob")

Ara_PXD008355_USI <- subset(Ara_PXD008355,select=-c(PTM_final_prob))

Spectrum_Ara_Pform_Mean <- merge(Ara_PXD008355_USI,Ara_PSM_Mean,by="All_USI")

Bin_Ara_PSM_Mean <- binAdjustAra(Spectrum_Ara_Pform_Mean)

BinAdj_Ara_Pform_MM <- aggregate(NewScore3 ~ Peptide_mod + Peptidoform, data = Bin_Ara_PSM_Mean, max)

names(BinAdj_Ara_Pform_MM) <- c("Peptide_mod", "Peptidoform" ,"BinMeanScore")

Total_BinAdj_Ara_Pform_MM <- merge(BinAdj_Ara_Pform_MM,Bin_Ara_PSM_Mean,by="Peptidoform")

BinAdj_Ara_Pform_MM <- Total_BinAdj_Ara_Pform_MM %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

BinAdj_Ara_Pform_MM$FinalScore3 <- NULL


# THE SAME FOR UNADJUSTED DATA #
################################

# Max

UnAdj_Ara_Pform_Max <- aggregate(PTM_final_prob ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Ara_PSM, max)

total_UnAdj_Ara_Pform_Max <- merge(UnAdj_Ara_Pform_Max,Ara_PXD008355,by="Peptidoform")

UnAdj_Ara_Pform_Max <- total_UnAdj_Ara_Pform_Max %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean

UnAdj_Ara_Pform_Mean <- aggregate(PTM_final_prob ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Ara_PSM, mean)

total_UnAdj_Ara_Pform_Mean <- merge(UnAdj_Ara_Pform_Mean,Ara_PXD008355,by="Peptidoform")

UnAdj_Ara_Pform_Mean <- total_UnAdj_Ara_Pform_Mean %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean of the PSM first

UnAdj_Ara_PSM_Mean <- aggregate(PTM_final_prob ~ All_USI , data = BinAdj_Ara_PSM, mean)

names(UnAdj_Ara_PSM_Mean) <- c("All_USI", "BinMeanScore")

Spectrum_UnAdj_Ara_Pform_Mean <- merge(UnAdj_Ara_PSM_Mean,BinAdj_Ara_PSM,by="All_USI")

UnAdj_Ara_Pform_MM <- aggregate(BinMeanScore ~ Peptide_mod + Peptidoform, data = Spectrum_UnAdj_Ara_Pform_Mean, max)

names(UnAdj_Ara_Pform_MM) <- c("Peptide_mod", "Peptidoform" ,"BinMeanScore")

Total_UnAdj_Ara_Pform_MM <- merge(UnAdj_Ara_Pform_MM,BinAdj_Ara_PSM,by="Peptidoform")

UnAdj_Ara_Pform_MM <- Total_UnAdj_Ara_Pform_MM %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

UnAdj_Ara_Pform_MM$FinalScore3 <- NULL

# EXtracting info for FLR calculation #


Not_adj_Ara_Pform_Max <- cbind.data.frame(UnAdj_Ara_Pform_Max$Peptidoform, UnAdj_Ara_Pform_Max$Peptide_mod,UnAdj_Ara_Pform_Max$Peptide,UnAdj_Ara_Pform_Max$PTM_positions,UnAdj_Ara_Pform_Max$PTM_final_prob.x)
names(Not_adj_Ara_Pform_Max)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Not_adj_Ara_Pform_Mean <- cbind.data.frame(UnAdj_Ara_Pform_Mean$Peptidoform, UnAdj_Ara_Pform_Mean$Peptide_mod,UnAdj_Ara_Pform_Mean$Peptide,UnAdj_Ara_Pform_Mean$PTM_positions,UnAdj_Ara_Pform_Mean$PTM_final_prob.x)
names(Not_adj_Ara_Pform_Mean)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Not_adj_Ara_Pform_MM <- cbind.data.frame(UnAdj_Ara_Pform_MM$Peptidoform, UnAdj_Ara_Pform_MM$Peptide_mod.x,UnAdj_Ara_Pform_MM$Peptide,UnAdj_Ara_Pform_MM$PTM_positions,UnAdj_Ara_Pform_MM$BinMeanScore)
names(Not_adj_Ara_Pform_MM)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")


Prod_adj_Ara_Pform <- cbind.data.frame(ProdAdj_Ara_Pform$Peptidoform, ProdAdj_Ara_Pform$Peptide_mod,ProdAdj_Ara_Pform$Peptide,ProdAdj_Ara_Pform$PTM_positions,ProdAdj_Ara_Pform$Prod_Adjusted_Score)
names(Prod_adj_Ara_Pform)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score")


Bin_adj_Ara_Pform_Max <- cbind.data.frame(BinAdj_Ara_Pform_Max$Peptidoform, BinAdj_Ara_Pform_Max$Peptide_mod,BinAdj_Ara_Pform_Max$Peptide,BinAdj_Ara_Pform_Max$PTM_positions,BinAdj_Ara_Pform_Max$NewScore3)
names(Bin_adj_Ara_Pform_Max)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

Bin_adj_Ara_Pform_Mean <- cbind.data.frame(BinAdj_Ara_Pform_Mean$Peptidoform, BinAdj_Ara_Pform_Mean$Peptide_mod,BinAdj_Ara_Pform_Mean$Peptide,BinAdj_Ara_Pform_Mean$PTM_positions,BinAdj_Ara_Pform_Mean$NewScore3)
names(Bin_adj_Ara_Pform_Mean)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

Bin_adj_Ara_Pform_MM <- cbind.data.frame(BinAdj_Ara_Pform_MM$Peptidoform, BinAdj_Ara_Pform_MM$Peptide_mod.x,BinAdj_Ara_Pform_MM$Peptide,BinAdj_Ara_Pform_MM$PTM_positions,BinAdj_Ara_Pform_MM$BinMeanScore)
names(Bin_adj_Ara_Pform_MM)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")



FLR_Not_adj_Ara_Pform_Max<-FLR_NotAdj(Not_adj_Ara_Pform_Max)

FLR_Not_adj_Ara_Pform_Mean<-FLR_NotAdj(Not_adj_Ara_Pform_Mean)

FLR_Not_adj_Ara_Pform_MM<-FLR_NotAdj(Not_adj_Ara_Pform_MM)

FLR_Prod_adj_Ara_Pform<-FLR_Prod_Adj(Prod_adj_Ara_Pform)

FLR_Bin_adj_Ara_Pform_Max<-FLR_Adj(Bin_adj_Ara_Pform_Max)

FLR_Bin_adj_Ara_Pform_Mean<-FLR_Adj(Bin_adj_Ara_Pform_Mean)

FLR_Bin_adj_Ara_Pform_MM<-FLR_Adj(Bin_adj_Ara_Pform_MM)


# Plotting Arabidopsis #


Graph.data.Ara <- data.frame(cbind(FLR_Not_adj_Ara_Pform_Max$Rnumber,FLR_Not_adj_Ara_Pform_Max$FLR_Unadjusted,FLR_Not_adj_Ara_Pform_Mean$FLR_Unadjusted,FLR_Not_adj_Ara_Pform_MM$FLR_Unadjusted,
                                   FLR_Bin_adj_Ara_Pform_Max$FLR_Adj_Score, FLR_Bin_adj_Ara_Pform_Mean$FLR_Adj_Score,FLR_Bin_adj_Ara_Pform_MM$FLR_Adj_Score, FLR_Prod_adj_Ara_Pform$FLR_Adj_Score))

colnames(Graph.data.Ara) <- c("Order_of_Scores","Unadjusted_PformMax","Unadjusted_PformMean ", "Unadjusted_PformMM", 
                              "Binomial_PformMax", "Binomial_PformMean", "Binomial_PformMM", "Product_Pform")

df_Graph_Ara <- gather(Graph.data.Ara,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Ara, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Ara_PXD008355")



# Human #
#########


Hum_PXD000612$Peptidoform <- paste0(Hum_PXD000612$Peptide_mod,"_",Hum_PXD000612$PTM_positions)

BinAdj_Hum_PSM<-binAdjustPSM(Hum_PXD000612)

ProdAdj_Hum_PSM <- BinAdj_Hum_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjusted_Score = (1-prod(1-PTM_final_prob)))

ProdAdj_Hum_Pform <- ProdAdj_Hum_PSM[!duplicated(ProdAdj_Hum_PSM$Peptidoform), ]


#  Max

BinAdj_Hum_Pform_Max <- aggregate(NewScore3 ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Hum_PSM, max)

total_BinAdj_Hum_Pform_Max <- merge(BinAdj_Hum_Pform_Max,Hum_PXD000612,by="Peptidoform")

BinAdj_Hum_Pform_Max <- total_BinAdj_Hum_Pform_Max %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean

BinAdj_Hum_Pform_Mean <- aggregate(NewScore3 ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Hum_PSM, mean)

total_BinAdj_Hum_Pform_Mean <- merge(BinAdj_Hum_Pform_Mean,Hum_PXD000612,by="Peptidoform")

BinAdj_Hum_Pform_Mean <- total_BinAdj_Hum_Pform_Mean %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean of the PSM first

Hum_PSM_Mean <- aggregate(PTM_final_prob ~ Spectrum , data = Hum_PXD000612, mean)

names(Hum_PSM_Mean) <- c("Spectrum", "PTM_final_prob")

Hum_PXD000612_USI <- subset(Hum_PXD000612,select=-c(PTM_final_prob))

Spectrum_Hum_Pform_Mean <- merge(Hum_PXD000612_USI,Hum_PSM_Mean,by="Spectrum")

Bin_Hum_PSM_Mean <- binAdjustPSM(Spectrum_Hum_Pform_Mean)

BinAdj_Hum_Pform_MM <- aggregate(NewScore3 ~ Peptide_mod + Peptidoform, data = Bin_Hum_PSM_Mean, max)

names(BinAdj_Hum_Pform_MM) <- c("Peptide_mod", "Peptidoform" ,"BinMeanScore")

Total_BinAdj_Hum_Pform_MM <- merge(BinAdj_Hum_Pform_MM,Bin_Hum_PSM_Mean,by="Peptidoform")

BinAdj_Hum_Pform_MM <- Total_BinAdj_Hum_Pform_MM %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

BinAdj_Hum_Pform_MM$FinalScore3 <- NULL


# THE SAME FOR UNADJUSTED DATA #
################################

# Max

UnAdj_Hum_Pform_Max <- aggregate(PTM_final_prob ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Hum_PSM, max)

total_UnAdj_Hum_Pform_Max <- merge(UnAdj_Hum_Pform_Max,Hum_PXD000612,by="Peptidoform")

UnAdj_Hum_Pform_Max <- total_UnAdj_Hum_Pform_Max %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean

UnAdj_Hum_Pform_Mean <- aggregate(PTM_final_prob ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Hum_PSM, mean)

total_UnAdj_Hum_Pform_Mean <- merge(UnAdj_Hum_Pform_Mean,Hum_PXD000612,by="Peptidoform")

UnAdj_Hum_Pform_Mean <- total_UnAdj_Hum_Pform_Mean %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean of the PSM first

UnAdj_Hum_PSM_Mean <- aggregate(PTM_final_prob ~ Spectrum , data = BinAdj_Hum_PSM, mean)

names(UnAdj_Hum_PSM_Mean) <- c("Spectrum", "BinMeanScore")

Spectrum_UnAdj_Hum_Pform_Mean <- merge(UnAdj_Hum_PSM_Mean,BinAdj_Hum_PSM,by="Spectrum")

UnAdj_Hum_Pform_MM <- aggregate(BinMeanScore ~ Peptide_mod + Peptidoform, data = Spectrum_UnAdj_Hum_Pform_Mean, max)

names(UnAdj_Hum_Pform_MM) <- c("Peptide_mod", "Peptidoform" ,"BinMeanScore")

Total_UnAdj_Hum_Pform_MM <- merge(UnAdj_Hum_Pform_MM,BinAdj_Hum_PSM,by="Peptidoform")

UnAdj_Hum_Pform_MM <- Total_UnAdj_Hum_Pform_MM %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

UnAdj_Hum_Pform_MM$FinalScore3 <- NULL

# EXtracting info for FLR calculation #


Not_adj_Hum_Pform_Max <- cbind.data.frame(UnAdj_Hum_Pform_Max$Peptidoform, UnAdj_Hum_Pform_Max$Peptide_mod,UnAdj_Hum_Pform_Max$Peptide,UnAdj_Hum_Pform_Max$PTM_positions,UnAdj_Hum_Pform_Max$PTM_final_prob.x)
names(Not_adj_Hum_Pform_Max)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Not_adj_Hum_Pform_Mean <- cbind.data.frame(UnAdj_Hum_Pform_Mean$Peptidoform, UnAdj_Hum_Pform_Mean$Peptide_mod,UnAdj_Hum_Pform_Mean$Peptide,UnAdj_Hum_Pform_Mean$PTM_positions,UnAdj_Hum_Pform_Mean$PTM_final_prob.x)
names(Not_adj_Hum_Pform_Mean)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Not_adj_Hum_Pform_MM <- cbind.data.frame(UnAdj_Hum_Pform_MM$Peptidoform, UnAdj_Hum_Pform_MM$Peptide_mod.x,UnAdj_Hum_Pform_MM$Peptide,UnAdj_Hum_Pform_MM$PTM_positions,UnAdj_Hum_Pform_MM$BinMeanScore)
names(Not_adj_Hum_Pform_MM)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score")

Prod_adj_Hum_Pform <- cbind.data.frame(ProdAdj_Hum_Pform$Peptidoform, ProdAdj_Hum_Pform$Peptide_mod,ProdAdj_Hum_Pform$Peptide,ProdAdj_Hum_Pform$PTM_positions,ProdAdj_Hum_Pform$Prod_Adjusted_Score)
names(Prod_adj_Hum_Pform)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score")

Bin_adj_Hum_Pform_Max <- cbind.data.frame(BinAdj_Hum_Pform_Max$Peptidoform, BinAdj_Hum_Pform_Max$Peptide_mod,BinAdj_Hum_Pform_Max$Peptide,BinAdj_Hum_Pform_Max$PTM_positions,BinAdj_Hum_Pform_Max$NewScore3)
names(Bin_adj_Hum_Pform_Max)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

Bin_adj_Hum_Pform_Mean <- cbind.data.frame(BinAdj_Hum_Pform_Mean$Peptidoform, BinAdj_Hum_Pform_Mean$Peptide_mod,BinAdj_Hum_Pform_Mean$Peptide,BinAdj_Hum_Pform_Mean$PTM_positions,BinAdj_Hum_Pform_Mean$NewScore3)
names(Bin_adj_Hum_Pform_Mean)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")

Bin_adj_Hum_Pform_MM <- cbind.data.frame(BinAdj_Hum_Pform_MM$Peptidoform, BinAdj_Hum_Pform_MM$Peptide_mod.x,BinAdj_Hum_Pform_MM$Peptide,BinAdj_Hum_Pform_MM$PTM_positions,BinAdj_Hum_Pform_MM$BinMeanScore)
names(Bin_adj_Hum_Pform_MM)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score")



FLR_Not_adj_Hum_Pform_Max<-FLR_NotAdj(Not_adj_Hum_Pform_Max)

FLR_Not_adj_Hum_Pform_Mean<-FLR_NotAdj(Not_adj_Hum_Pform_Mean)

FLR_Not_adj_Hum_Pform_MM<-FLR_NotAdj(Not_adj_Hum_Pform_MM)

FLR_Prod_adj_Hum_Pform<-FLR_Prod_Adj(Prod_adj_Hum_Pform)

FLR_Bin_adj_Hum_Pform_Max<-FLR_Adj(Bin_adj_Hum_Pform_Max)

FLR_Bin_adj_Hum_Pform_Mean<-FLR_Adj(Bin_adj_Hum_Pform_Mean)

FLR_Bin_adj_Hum_Pform_MM<-FLR_Adj(Bin_adj_Hum_Pform_MM)


# Plotting Human #


Graph.data.Hum <- data.frame(cbind(FLR_Not_adj_Hum_Pform_Max$Rnumber,FLR_Not_adj_Hum_Pform_Max$FLR_Unadjusted,FLR_Not_adj_Hum_Pform_Mean$FLR_Unadjusted,FLR_Not_adj_Hum_Pform_MM$FLR_Unadjusted,
                                   FLR_Bin_adj_Hum_Pform_Max$FLR_Adj_Score, FLR_Bin_adj_Hum_Pform_Mean$FLR_Adj_Score,FLR_Bin_adj_Hum_Pform_MM$FLR_Adj_Score,FLR_Prod_adj_Hum_Pform$FLR_Adj_Score))

colnames(Graph.data.Hum) <- c("Order_of_Scores","Unadjusted_PformMax","Unadjusted_PformMean ", "Unadjusted_PformMM", 
                              "Binomial_PformMax", "Binomial_PformMean", "Binomial_PformMM", "Product_Pform")

df_Graph_Hum <- gather(Graph.data.Hum,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Hum, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Hum_PXD000612")



# Synthetic PXD000138 #
#######################



Syn_PXD000138$Peptidoform <- paste0(Syn_PXD000138$Peptide_mod,"_",Syn_PXD000138$PTM_positions)

BinAdj_Syn138_PSM<-binAdjustSyn(Syn_PXD000138)

ProdAdj_Syn138_PSM <- BinAdj_Syn138_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjusted_Score = (1-prod(1-PTM_final_prob)))

ProdAdj_Syn138_Pform <- ProdAdj_Syn138_PSM[!duplicated(ProdAdj_Syn138_PSM$Peptidoform), ]


# Max

BinAdj_Syn138_Pform_Max <- aggregate(NewScore3 ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Syn138_PSM, max)

total_BinAdj_Syn138_Pform_Max <- merge(BinAdj_Syn138_Pform_Max,Syn_PXD000138,by="Peptidoform")

BinAdj_Syn138_Pform_Max <- total_BinAdj_Syn138_Pform_Max %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# mean

BinAdj_Syn138_Pform_Mean <- aggregate(NewScore3 ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Syn138_PSM, mean)

total_BinAdj_Syn138_Pform_Mean <- merge(BinAdj_Syn138_Pform_Mean,Syn_PXD000138,by="Peptidoform")

BinAdj_Syn138_Pform_Mean <- total_BinAdj_Syn138_Pform_Mean %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean of the PSM first

UnAdj_Syn138_PSM_Mean <- aggregate(PTM_final_prob ~ Spectrum , data = Syn_PXD000138, mean)

names(UnAdj_Syn138_PSM_Mean) <- c("Spectrum", "PTM_final_prob")

Syn_PXD000138_USI2 <- dplyr::select(Syn_PXD000138, -c(PTM_final_prob))

Spectrum_UnAdj_Syn138_Pform_Mean <- merge(UnAdj_Syn138_PSM_Mean,Syn_PXD000138_USI2,by="Spectrum")

Bin_Syn138_USI <- binAdjustSyn(Spectrum_UnAdj_Syn138_Pform_Mean)


BinAdj_Syn138_Pform_1 <- aggregate(NewScore3 ~ Peptide_mod + Peptidoform , data = Bin_Syn138_USI, max)

names(BinAdj_Syn138_Pform_1) <- c("Peptide_mod", "Peptidoform" ,"BinMeanScore")

Total_BinAdj_Syn138_Pform_MM <- merge(BinAdj_Syn138_Pform_1,BinAdj_Syn138_PSM,by="Peptidoform")

BinAdj_Syn138_Pform_MM <- Total_BinAdj_Syn138_Pform_MM %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

BinAdj_Syn138_Pform_MM$FinalScore3 <- NULL


# THE SAME FOR UNADJUSTED DATA #
################################

# Max

UnAdj_Syn138_Pform_Max <- aggregate(PTM_final_prob ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Syn138_PSM, max)

total_UnAdj_Syn138_Pform_Max <- merge(UnAdj_Syn138_Pform_Max,Syn_PXD000138,by="Peptidoform")

UnAdj_Syn138_Pform_Max <- total_UnAdj_Syn138_Pform_Max %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean

UnAdj_Syn138_Pform_Mean <- aggregate(PTM_final_prob ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Syn138_PSM, mean)

total_UnAdj_Syn138_Pform_Mean <- merge(UnAdj_Syn138_Pform_Mean,Syn_PXD000138,by="Peptidoform")

UnAdj_Syn138_Pform_Mean <- total_UnAdj_Syn138_Pform_Mean %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean of the PSM first


UnAdj_Syn138_Pform_1 <- aggregate(PTM_final_prob ~ Peptide_mod + Peptidoform, data = Spectrum_UnAdj_Syn138_Pform_Mean, max)

names(UnAdj_Syn138_Pform_1) <- c("Peptide_mod", "Peptidoform" ,"BinMeanScore")

Total_UnAdj_Syn138_Pform_MM <- merge(UnAdj_Syn138_Pform_1,BinAdj_Syn138_PSM,by="Peptidoform")

UnAdj_Syn138_Pform_MM <- Total_UnAdj_Syn138_Pform_MM %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

UnAdj_Syn138_Pform_MM$FinalScore3 <- NULL

# EXtracting info for FLR calculation #


Not_adj_Syn138_Pform_Max <- cbind.data.frame(UnAdj_Syn138_Pform_Max$Peptidoform, UnAdj_Syn138_Pform_Max$Peptide_mod,UnAdj_Syn138_Pform_Max$Peptide,UnAdj_Syn138_Pform_Max$PTM_positions,UnAdj_Syn138_Pform_Max$PTM_final_prob.x, UnAdj_Syn138_Pform_Max$Incorrect.count)
names(Not_adj_Syn138_Pform_Max)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score", "Incorrect.count")

Not_adj_Syn138_Pform_Mean <- cbind.data.frame(UnAdj_Syn138_Pform_Mean$Peptidoform, UnAdj_Syn138_Pform_Mean$Peptide_mod,UnAdj_Syn138_Pform_Mean$Peptide,UnAdj_Syn138_Pform_Mean$PTM_positions,UnAdj_Syn138_Pform_Mean$PTM_final_prob.x, UnAdj_Syn138_Pform_Mean$Incorrect.count)
names(Not_adj_Syn138_Pform_Mean)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score", "Incorrect.count")

Not_adj_Syn138_Pform_MM <- cbind.data.frame(UnAdj_Syn138_Pform_MM$Peptidoform, UnAdj_Syn138_Pform_MM$Peptide_mod.x,UnAdj_Syn138_Pform_MM$Peptide,UnAdj_Syn138_Pform_MM$PTM_positions,UnAdj_Syn138_Pform_MM$BinMeanScore, UnAdj_Syn138_Pform_MM$Incorrect.count)
names(Not_adj_Syn138_Pform_MM)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score", "Incorrect.count")

Prod_adj_Syn138_Pform <- cbind.data.frame(ProdAdj_Syn138_Pform$Peptidoform, ProdAdj_Syn138_Pform$Peptide_mod,ProdAdj_Syn138_Pform$Peptide,ProdAdj_Syn138_Pform$PTM_positions,ProdAdj_Syn138_Pform$Prod_Adjusted_Score, ProdAdj_Syn138_Pform$Incorrect.count)
names(Prod_adj_Syn138_Pform)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score", "Incorrect.count")

Bin_adj_Syn138_Pform_Max <- cbind.data.frame(BinAdj_Syn138_Pform_Max$Peptidoform, BinAdj_Syn138_Pform_Max$Peptide_mod,BinAdj_Syn138_Pform_Max$Peptide,BinAdj_Syn138_Pform_Max$PTM_positions,BinAdj_Syn138_Pform_Max$NewScore3, BinAdj_Syn138_Pform_Max$Incorrect.count)
names(Bin_adj_Syn138_Pform_Max)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score", "Incorrect.count")

Bin_adj_Syn138_Pform_Mean <- cbind.data.frame(BinAdj_Syn138_Pform_Mean$Peptidoform, BinAdj_Syn138_Pform_Mean$Peptide_mod,BinAdj_Syn138_Pform_Mean$Peptide,BinAdj_Syn138_Pform_Mean$PTM_positions,BinAdj_Syn138_Pform_Mean$NewScore3, BinAdj_Syn138_Pform_Mean$Incorrect.count)
names(Bin_adj_Syn138_Pform_Mean)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score", "Incorrect.count")

Bin_adj_Syn138_Pform_MM <- cbind.data.frame(BinAdj_Syn138_Pform_MM$Peptidoform, BinAdj_Syn138_Pform_MM$Peptide_mod.x,BinAdj_Syn138_Pform_MM$Peptide,BinAdj_Syn138_Pform_MM$PTM_positions,BinAdj_Syn138_Pform_MM$BinMeanScore, BinAdj_Syn138_Pform_MM$Incorrect.count)
names(Bin_adj_Syn138_Pform_MM)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score", "Incorrect.count")



FLR_Not_adj_Syn138_Pform_Max<-FLR_NotAdj(Not_adj_Syn138_Pform_Max)

FLR_Not_adj_Syn138_Pform_Mean<-FLR_NotAdj(Not_adj_Syn138_Pform_Mean)

FLR_Not_adj_Syn138_Pform_MM<-FLR_NotAdj(Not_adj_Syn138_Pform_MM)

FLR_Prod_adj_Syn138_Pform<-FLR_Prod_Adj(Prod_adj_Syn138_Pform)

FLR_Bin_adj_Syn138_Pform_Max<-FLR_Adj(Bin_adj_Syn138_Pform_Max)

FLR_Bin_adj_Syn138_Pform_Mean<-FLR_Adj(Bin_adj_Syn138_Pform_Mean)

FLR_Bin_adj_Syn138_Pform_MM<-FLR_Adj(Bin_adj_Syn138_Pform_MM)



# Plotting Synthetic PXD000138 #


Graph.data.Syn138 <- data.frame(cbind(FLR_Not_adj_Syn138_Pform_Max$Rnumber,FLR_Not_adj_Syn138_Pform_Max$FLR_Unadjusted,FLR_Not_adj_Syn138_Pform_Mean$FLR_Unadjusted,FLR_Not_adj_Syn138_Pform_MM$FLR_Unadjusted,
                                      FLR_Bin_adj_Syn138_Pform_Max$FLR_Adj_Score, FLR_Bin_adj_Syn138_Pform_Mean$FLR_Adj_Score,FLR_Bin_adj_Syn138_Pform_MM$FLR_Adj_Score, FLR_Prod_adj_Syn138_Pform$FLR_Adj_Score))

colnames(Graph.data.Syn138) <- c("Order_of_Scores","Unadjusted_PformMax","Unadjusted_PformMean ", "Unadjusted_PformMM", 
                                 "Binomial_PformMax", "Binomial_PformMean", "Binomial_PformMM", "Product_Pform")

df_Graph_Syn138 <- gather(Graph.data.Syn138,Method ,FLR, -Order_of_Scores)



qplot(data=df_Graph_Syn138, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Syn_PXD000138")




# Synthetic PXD007058 #
#######################


Syn_PXD007058$Peptidoform <- paste0(Syn_PXD007058$Peptide_mod,"_",Syn_PXD007058$PTM_positions)

Syn_PXD007058$PTM_final_prob <- Syn_PXD007058$PTM_score*Syn_PXD007058$Score
Syn_PXD007058$All_Proteins <- Syn_PXD007058$Protein
Syn_PXD007058$Spectrum <- Syn_PXD007058$All_USI

BinAdj_Syn7058_PSM<-binAdjustSyn(Syn_PXD007058)

ProdAdj_Syn7058_PSM <- BinAdj_Syn7058_PSM %>% group_by(PROTEIN_LOC) %>% mutate(Prod_Adjusted_Score = (1-prod(1-PTM_final_prob)))

ProdAdj_Syn7058_Pform <- ProdAdj_Syn7058_PSM[!duplicated(ProdAdj_Syn7058_PSM$Peptidoform), ]


# Max

BinAdj_Syn7058_Pform_Max <- aggregate(NewScore3 ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Syn7058_PSM, max)

total_BinAdj_Syn7058_Pform_Max <- merge(BinAdj_Syn7058_Pform_Max,Syn_PXD007058,by="Peptidoform")

BinAdj_Syn7058_Pform_Max <- total_BinAdj_Syn7058_Pform_Max %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean

BinAdj_Syn7058_Pform_Mean <- aggregate(NewScore3 ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Syn7058_PSM, mean)

total_BinAdj_Syn7058_Pform_Mean <- merge(BinAdj_Syn7058_Pform_Mean,Syn_PXD007058,by="Peptidoform")

BinAdj_Syn7058_Pform_Mean <- total_BinAdj_Syn7058_Pform_Mean %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean of the PSM first

UnAdj_Syn7058_PSM_Mean <- aggregate(PTM_final_prob ~ All_USI , data = Syn_PXD007058, mean)

names(UnAdj_Syn7058_PSM_Mean) <- c("All_USI", "PTM_final_prob")

Syn_PXD007058_USI2 <- dplyr::select(Syn_PXD007058, -c(PTM_final_prob))

Spectrum_UnAdj_Syn7058_Pform_Mean <- merge(UnAdj_Syn7058_PSM_Mean,Syn_PXD007058_USI2,by="All_USI")

Bin_Syn7058_USI <- binAdjustSyn(Spectrum_UnAdj_Syn7058_Pform_Mean)


BinAdj_Syn7058_Pform_1 <- aggregate(NewScore3 ~ Peptide_mod + Peptidoform , data = Bin_Syn7058_USI, max)

names(BinAdj_Syn7058_Pform_1) <- c("Peptide_mod", "Peptidoform" ,"BinMeanScore")

Total_BinAdj_Syn7058_Pform_MM <- merge(BinAdj_Syn7058_Pform_1,BinAdj_Syn7058_PSM,by="Peptidoform")

BinAdj_Syn7058_Pform_MM <- Total_BinAdj_Syn7058_Pform_MM %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

BinAdj_Syn7058_Pform_MM$FinalScore3 <- NULL

# THE SAME FOR UNADJUSTED DATA #
################################

# Max

UnAdj_Syn7058_Pform_Max <- aggregate(PTM_final_prob ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Syn7058_PSM, max)

total_UnAdj_Syn7058_Pform_Max <- merge(UnAdj_Syn7058_Pform_Max,Syn_PXD007058,by="Peptidoform")

UnAdj_Syn7058_Pform_Max <- total_UnAdj_Syn7058_Pform_Max %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean

UnAdj_Syn7058_Pform_Mean <- aggregate(PTM_final_prob ~ Peptidoform+PROTEIN_LOC, data = BinAdj_Syn7058_PSM, mean)

total_UnAdj_Syn7058_Pform_Mean <- merge(UnAdj_Syn7058_Pform_Mean,Syn_PXD007058,by="Peptidoform")

UnAdj_Syn7058_Pform_Mean <- total_UnAdj_Syn7058_Pform_Mean %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

# Mean of the PSM first

UnAdj_Syn7058_Pform_1 <- aggregate(PTM_final_prob ~ Peptide_mod + Peptidoform, data = Spectrum_UnAdj_Syn7058_Pform_Mean, max)

names(UnAdj_Syn7058_Pform_1) <- c("Peptide_mod", "Peptidoform" ,"BinMeanScore")

Total_UnAdj_Syn7058_Pform_MM <- merge(UnAdj_Syn7058_Pform_1,BinAdj_Syn7058_PSM,by="Peptidoform")

UnAdj_Syn7058_Pform_MM <- Total_UnAdj_Syn7058_Pform_MM %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()

UnAdj_Syn7058_Pform_MM$FinalScore3 <- NULL

# EXtracting info for FLR calculation #



Not_adj_Syn7058_Pform_Max <- cbind.data.frame(UnAdj_Syn7058_Pform_Max$Peptidoform, UnAdj_Syn7058_Pform_Max$Peptide_mod,UnAdj_Syn7058_Pform_Max$Peptide,UnAdj_Syn7058_Pform_Max$PTM_positions,UnAdj_Syn7058_Pform_Max$PTM_final_prob.x, UnAdj_Syn7058_Pform_Max$Synthetic_full_match)
names(Not_adj_Syn7058_Pform_Max)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score", "Incorrect.count")

Not_adj_Syn7058_Pform_Mean <- cbind.data.frame(UnAdj_Syn7058_Pform_Mean$Peptidoform, UnAdj_Syn7058_Pform_Mean$Peptide_mod,UnAdj_Syn7058_Pform_Mean$Peptide,UnAdj_Syn7058_Pform_Mean$PTM_positions,UnAdj_Syn7058_Pform_Mean$PTM_final_prob.x, UnAdj_Syn7058_Pform_Mean$Synthetic_full_match)
names(Not_adj_Syn7058_Pform_Mean)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score", "Incorrect.count")

Not_adj_Syn7058_Pform_MM <- cbind.data.frame(UnAdj_Syn7058_Pform_MM$Peptidoform, UnAdj_Syn7058_Pform_MM$Peptide_mod.x,UnAdj_Syn7058_Pform_MM$Peptide,UnAdj_Syn7058_Pform_MM$PTM_positions,UnAdj_Syn7058_Pform_MM$BinMeanScore, UnAdj_Syn7058_Pform_MM$Synthetic_full_match)
names(Not_adj_Syn7058_Pform_MM)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Unadjusted_Score", "Incorrect.count")

Prod_adj_Syn7058_Pform <- cbind.data.frame(ProdAdj_Syn7058_Pform$Peptidoform, ProdAdj_Syn7058_Pform$Peptide_mod,ProdAdj_Syn7058_Pform$Peptide,ProdAdj_Syn7058_Pform$PTM_positions,ProdAdj_Syn7058_Pform$Prod_Adjusted_Score, ProdAdj_Syn7058_Pform$Synthetic_full_match)
names(Prod_adj_Syn7058_Pform)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Prod_Adjusted_Score", "Incorrect.count")

Bin_adj_Syn7058_Pform_Max <- cbind.data.frame(BinAdj_Syn7058_Pform_Max$Peptidoform, BinAdj_Syn7058_Pform_Max$Peptide_mod,BinAdj_Syn7058_Pform_Max$Peptide,BinAdj_Syn7058_Pform_Max$PTM_positions,BinAdj_Syn7058_Pform_Max$NewScore3, BinAdj_Syn7058_Pform_Max$Synthetic_full_match)
names(Bin_adj_Syn7058_Pform_Max)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score", "Incorrect.count")

Bin_adj_Syn7058_Pform_Mean <- cbind.data.frame(BinAdj_Syn7058_Pform_Mean$Peptidoform, BinAdj_Syn7058_Pform_Mean$Peptide_mod,BinAdj_Syn7058_Pform_Mean$Peptide,BinAdj_Syn7058_Pform_Mean$PTM_positions,BinAdj_Syn7058_Pform_Mean$NewScore3, BinAdj_Syn7058_Pform_Mean$Synthetic_full_match)
names(Bin_adj_Syn7058_Pform_Mean)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score", "Incorrect.count")

Bin_adj_Syn7058_Pform_MM <- cbind.data.frame(BinAdj_Syn7058_Pform_MM$Peptidoform, BinAdj_Syn7058_Pform_MM$Peptide_mod.x,BinAdj_Syn7058_Pform_MM$Peptide,BinAdj_Syn7058_Pform_MM$PTM_positions,BinAdj_Syn7058_Pform_MM$BinMeanScore, BinAdj_Syn7058_Pform_MM$Synthetic_full_match)
names(Bin_adj_Syn7058_Pform_MM)<-c("Peptidoform","Peptide_mod", "Peptide", "PTM_positions" ,"Bin_Adjusted_Score", "Incorrect.count")


FLR_Not_adj_Syn7058_Pform_Max<-FLR_NotAdj(Not_adj_Syn7058_Pform_Max)

FLR_Not_adj_Syn7058_Pform_Mean<-FLR_NotAdj(Not_adj_Syn7058_Pform_Mean)

FLR_Not_adj_Syn7058_Pform_MM<-FLR_NotAdj(Not_adj_Syn7058_Pform_MM)

FLR_Prod_adj_Syn7058_Pform<-FLR_Prod_Adj(Prod_adj_Syn7058_Pform)

FLR_Bin_adj_Syn7058_Pform_Max<-FLR_Adj(Bin_adj_Syn7058_Pform_Max)

FLR_Bin_adj_Syn7058_Pform_Mean<-FLR_Adj(Bin_adj_Syn7058_Pform_Mean)

FLR_Bin_adj_Syn7058_Pform_MM<-FLR_Adj(Bin_adj_Syn7058_Pform_MM)



# Plotting Synthetic PXD007058 #


Graph.data.Syn7058 <- data.frame(cbind(FLR_Not_adj_Syn7058_Pform_Max$Rnumber,FLR_Not_adj_Syn7058_Pform_Max$FLR_Unadjusted,FLR_Not_adj_Syn7058_Pform_Mean$FLR_Unadjusted,FLR_Not_adj_Syn7058_Pform_MM$FLR_Unadjusted,
                                       FLR_Bin_adj_Syn7058_Pform_Max$FLR_Adj_Score, FLR_Bin_adj_Syn7058_Pform_Mean$FLR_Adj_Score,FLR_Bin_adj_Syn7058_Pform_MM$FLR_Adj_Score, FLR_Prod_adj_Syn7058_Pform$FLR_Adj_Score))

colnames(Graph.data.Syn7058) <- c("Order_of_Scores","Unadjusted_PformMax","Unadjusted_PformMean ", "Unadjusted_PformMM", 
                                  "Binomial_PformMax", "Binomial_PformMean", "Binomial_PformMM", "Product_Pform")


df_Graph_Syn7058 <- gather(Graph.data.Syn7058,Method ,FLR, -Order_of_Scores)

qplot(data=df_Graph_Syn7058, x=Order_of_Scores, y=FLR, geom='line', color=Method, main="Syn_PXD007058")

