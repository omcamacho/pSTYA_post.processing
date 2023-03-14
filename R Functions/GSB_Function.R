

# Gold, Silver, Bronze function #

GSB_Function <- function(dataf_01,data_05){
  
  
  dataf_01_ND <- dataf_01 [!duplicated(dataf_01[c("PROTEIN_LOC","dataset")]),]
  
  dataf_01_ND <-transform(dataf_01_ND,Frequency.Observed = ave(seq(nrow(dataf_01_ND)),PROTEIN_LOC,FUN=length))
  
  hist(dataf_01_ND$Frequency.Observed, main="Frequency of observations")
  
  Pform_Gold <- dataf_01_ND[dataf_01_ND$Frequency.Observed>1,]
  
  pform_Gold_UPS <- Pform_Gold [!duplicated(Pform_Gold[c("PROTEIN_LOC")]),]
  
  pform_Gold_UPS$cat<-"Gold"
  
  # Silver New_FLR_PEP<0.01 observed in 1 data set only and New_FLR_PEP<0.0 5 in 2 or more data sets#
  
  pform_Silver <- dataf_01_ND[dataf_01_ND$Frequency.Observed<2,]
  
  pform_Silver_U <- subset(pform_Silver, !(PROTEIN_LOC %in% pform_Gold_UPS$PROTEIN_LOC))
  
  pform_Silver_UPS <- pform_Silver_U [!duplicated(pform_Silver_U[c("PROTEIN_LOC")]),]
  
  pform_Silver_UPS$cat<-"Silver"
  
  # Bronze
  
  pform_Bronze <- data_05
  
  pform_GS <- dplyr::bind_rows(pform_Gold_UPS, pform_Silver_UPS)
  
  pform_Bronze_U <- subset(pform_Bronze, !(PROTEIN_LOC %in% pform_GS$PROTEIN_LOC))
  
  pform_Bronze_UPS <- pform_Bronze_U [!duplicated(pform_Bronze_U[c("PROTEIN_LOC")]),]
  
  pform_Bronze_UPS$cat<-"Bronze"
  
  pform_pooled <- dplyr::bind_rows(pform_Gold_UPS, pform_Silver_UPS, pform_Bronze_UPS)
  
  pform_Final <- pform_pooled [!duplicated(pform_pooled[c("PROTEIN_LOC")]),]
  
  return(pform_Final)
}