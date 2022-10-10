

# Gold, Silver, Bronze function #

GSB_Function <- function(dataf){
  
  AllRice_pform_Max <- dataf
  
  AllRice_pform_Max_01 <- AllRice_pform_Max[AllRice_pform_Max$New_FLR_PEP<0.01,]
  
  AllRice_pform_Max_01_ND <- AllRice_pform_Max_01 [!duplicated(AllRice_pform_Max_01[c("PROTEIN_LOC","dataset")]),]
  
  AllRice_pform_Max_01_ND <-transform(AllRice_pform_Max_01_ND,Frequency.Observed = ave(seq(nrow(AllRice_pform_Max_01_ND)),PROTEIN_LOC,FUN=length))
  
  hist(AllRice_pform_Max_01_ND$Frequency.Observed, main="Frequency of observations")
  
  AllRice_pform_Max_Gold <- AllRice_pform_Max_01_ND[AllRice_pform_Max_01_ND$Frequency.Observed>1,]
  
  AllRice_pform_Max_Gold_UPS <- AllRice_pform_Max_Gold [!duplicated(AllRice_pform_Max_Gold[c("PROTEIN_LOC")]),]
  
  AllRice_pform_Max_Gold_UPS$cat<-"Gold"
  
  # Silver New_FLR_PEP<0.01 observed in 1 data set only and New_FLR_PEP<0.0 5 in 2 or more data sets#
  
  AllRice_pform_Max_01_Silver <- AllRice_pform_Max_01_ND[AllRice_pform_Max_01_ND$Frequency.Observed<2,]
  
  hist(AllRice_pform_Max_01_Silver$Frequency.Observed, main="Frequency of observations")
  
  AllRice_pform_Max_05 <- AllRice_pform_Max[AllRice_pform_Max$New_FLR_PEP>=0.01 & AllRice_pform_Max$New_FLR_PEP<0.05,]
  
  AllRice_pform_Max_05_ND <- AllRice_pform_Max_05 [!duplicated(AllRice_pform_Max_05[c("PROTEIN_LOC","dataset")]),]
  
  AllRice_pform_Max_05_ND <-transform(AllRice_pform_Max_05_ND,Frequency.Observed = ave(seq(nrow(AllRice_pform_Max_05_ND)),PROTEIN_LOC,FUN=length))
  
  hist(AllRice_pform_Max_05_ND$Frequency.Observed, main="Frequency of observations")
  
  AllRice_pform_Max_05_Silver <- AllRice_pform_Max_05_ND[AllRice_pform_Max_05_ND$Frequency.Observed>1,]
  
  AllRice_pform_Max_Silver <- dplyr::bind_rows(AllRice_pform_Max_01_Silver, AllRice_pform_Max_05_Silver)
  
  AllRice_pform_Max_Silver_U <- subset(AllRice_pform_Max_Silver, !(PROTEIN_LOC %in% AllRice_pform_Max_Gold_UPS$PROTEIN_LOC))
  
  AllRice_pform_Max_Silver_UPS <- AllRice_pform_Max_Silver_U [!duplicated(AllRice_pform_Max_Silver_U[c("PROTEIN_LOC")]),]
  
  AllRice_pform_Max_Silver_UPS$cat<-"Silver"
  
  # Bronze
  
  AllRice_pform_Max_Bronze <- AllRice_pform_Max_05_ND[AllRice_pform_Max_05_ND$Frequency.Observed<2,]
  
  AllRice_pform_Max_GS <- dplyr::bind_rows(AllRice_pform_Max_Gold_UPS, AllRice_pform_Max_Silver_UPS)
  
  AllRice_pform_Max_Bronze_U <- subset(AllRice_pform_Max_Bronze, !(PROTEIN_LOC %in% AllRice_pform_Max_GS$PROTEIN_LOC))
  
  AllRice_pform_Max_Bronze_UPS <- AllRice_pform_Max_Bronze_U [!duplicated(AllRice_pform_Max_Bronze_U[c("PROTEIN_LOC")]),]
  
  AllRice_pform_Max_Bronze_UPS$cat<-"Bronze"
  
  AllRice_pform_Max_pooled <- dplyr::bind_rows(AllRice_pform_Max_Gold_UPS, AllRice_pform_Max_Silver_UPS, AllRice_pform_Max_Bronze_UPS)
  
  AllRice_pform_Max_Final <- AllRice_pform_Max_pooled [!duplicated(AllRice_pform_Max_pooled[c("PROTEIN_LOC")]),]
  
  return(AllRice_pform_Max_Final)
}