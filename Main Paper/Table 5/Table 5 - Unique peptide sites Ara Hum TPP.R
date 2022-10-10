

# Table 5

library("conflicted")
suppressPackageStartupMessages(library("tidyverse"))
conflict_prefer("filter", "dplyr")
#> [conflicted] Will prefer dplyr::filter over any other package
library("broom")
library("kableExtra")
library(plyr)
library(dplyr)
library(stringr)
library(useful)
library("data.table")
conflict_prefer("mutate", "dplyr")
conflict_prefer("rename", "dplyr")
library(epiDisplay)
library(ggplot2)
library(janitor)



source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/Binomial Scores Function Peptidoform level.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/PeptidoformToPeptide_Num Phospo Complex Algorithm.R')
source('C:/Users/OMCamacho/Desktop/Paper submission/R Functions/PeptidoformToPeptide_Num Phospo Complex Algorithm Only 1 pep peptide.R')

Ara_PXD008355 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD008355.csv')
Hum_PXD000612 <- read.csv(file = 'C:/Users/OMCamacho/Desktop/Paper submission/Data TPP pipeline/PXD000612.csv')

Ara_PXD008355$dataset <-"Ara_PXD008355"
Hum_PXD000612$dataset <-"Hum_PXD000612"


# Peptidoform level #
#####################

PepformAdj_Ara_PXD008355 <- binAdjustPform(Ara_PXD008355)
PepformAdj_Hum_PXD000612 <- binAdjustPform(Hum_PXD000612)

UniquePROT_Sites_Pform_Ara<- PepformAdj_Ara_PXD008355[!duplicated(PepformAdj_Ara_PXD008355$PROTEIN_LOC), ]

UniquePROT_Sites_Pform_Ara <- rename(UniquePROT_Sites_Pform_Ara,Bin_Adjusted_Score=NewScore3)

UniquePROT_Sites_Pform_Ara$num_Pho <- as.character(str_count(UniquePROT_Sites_Pform_Ara$Peptidoform, "Phospho"))

UniquePROT_Sites_Pform_Hum<- PepformAdj_Hum_PXD000612[!duplicated(PepformAdj_Hum_PXD000612$PROTEIN_LOC), ]

UniquePROT_Sites_Pform_Hum <- rename(UniquePROT_Sites_Pform_Hum,Bin_Adjusted_Score=NewScore3)

UniquePROT_Sites_Pform_Hum$num_Pho <- as.character(str_count(UniquePROT_Sites_Pform_Hum$Peptidoform, "Phospho"))


Table_Pform_Ara <- 
  UniquePROT_Sites_Pform_Ara %>%
  dplyr::count(num_Pho) %>%
  mutate(prop = prop.table(n))


Table_Pform_Hum <- 
  UniquePROT_Sites_Pform_Hum %>%
  dplyr::count(num_Pho) %>%
  mutate(prop = prop.table(n))


#  Peptide level - only one match per peptide - No grouping #
#############################################################

PeptideAdj_Ara_PXD008355_G1 <- PeptidoformToPeptide_CA_1(PepformAdj_Ara_PXD008355)
PeptideAdj_Hum_PXD000612_G1 <- PeptidoformToPeptide_CA_1(PepformAdj_Hum_PXD000612)


UniquePROT_Sites_Ptide_Ara_G1<- PeptideAdj_Ara_PXD008355_G1[!duplicated(PeptideAdj_Ara_PXD008355_G1$PROTEIN_LOC), ]

UniquePROT_Sites_Ptide_Ara_G1$num_Pho <- as.character(str_count(UniquePROT_Sites_Ptide_Ara_G1$Peptidoform, "Phospho"))

OnlyNumPho_Ara_G1 <- as.data.frame(UniquePROT_Sites_Ptide_Ara_G1$num_Pho)

names(OnlyNumPho_Ara_G1)[1] <- "num_Pho"

UniquePROT_Sites_Ptide_Hum_G1<- PeptideAdj_Hum_PXD000612_G1[!duplicated(PeptideAdj_Hum_PXD000612_G1$PROTEIN_LOC), ]

UniquePROT_Sites_Ptide_Hum_G1$num_Pho <- as.character(str_count(UniquePROT_Sites_Ptide_Hum_G1$Peptidoform, "Phospho"))

OnlyNumPho_Hum_G1 <- as.data.frame(UniquePROT_Sites_Ptide_Hum_G1$num_Pho)

names(OnlyNumPho_Hum_G1)[1] <- "num_Pho"

Table_Ptide_Ara_G1 <- 
  OnlyNumPho_Ara_G1 %>%
  dplyr::count(num_Pho) %>%
  mutate(prop = prop.table(n))

Table_Ptide_Hum_G1 <- 
  OnlyNumPho_Hum_G1 %>%
  dplyr::count(num_Pho) %>%
  mutate(prop = prop.table(n))


#  Peptide level - one match per number of phospho #
####################################################

PeptideAdj_Ara_PXD008355 <- PeptidoformToPeptide_CA(PepformAdj_Ara_PXD008355)
PeptideAdj_Hum_PXD000612 <- PeptidoformToPeptide_CA(PepformAdj_Hum_PXD000612)


UniquePROT_Sites_Ptide_Ara<- PeptideAdj_Ara_PXD008355[!duplicated(PeptideAdj_Ara_PXD008355$PROTEIN_LOC), ]

UniquePROT_Sites_Ptide_Ara$num_Pho <- as.character(str_count(UniquePROT_Sites_Ptide_Ara$Peptidoform, "Phospho"))

OnlyNumPho_Ara <- as.data.frame(UniquePROT_Sites_Ptide_Ara$num_Pho)

names(OnlyNumPho_Ara)[1] <- "num_Pho"

UniquePROT_Sites_Ptide_Hum<- PeptideAdj_Hum_PXD000612[!duplicated(PeptideAdj_Hum_PXD000612$PROTEIN_LOC), ]

UniquePROT_Sites_Ptide_Hum$num_Pho <- as.character(str_count(UniquePROT_Sites_Ptide_Hum$Peptidoform, "Phospho"))

OnlyNumPho_Hum <- as.data.frame(UniquePROT_Sites_Ptide_Hum$num_Pho)

names(OnlyNumPho_Hum)[1] <- "num_Pho"

Table_Ptide_Ara <- 
  OnlyNumPho_Ara %>%
  dplyr::count(num_Pho) %>%
  mutate(prop = prop.table(n))

Table_Ptide_Hum <- 
  OnlyNumPho_Hum%>%
  dplyr::count(num_Pho) %>%
  mutate(prop = prop.table(n))
