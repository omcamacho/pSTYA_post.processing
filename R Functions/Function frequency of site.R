library(dplyr)

FFreq <- function(dataf){

  dataf <- dataf[order(dataf$All_Proteins),]
  
matches <-list()
matches2 <-list()
N_Count=c()
for (protein in unique(dataf$All_Proteins)) {
  
  
  DataProtein <- dataf %>% filter(All_Proteins==protein)
  df <-NULL
  matches <-NULL
  sumdata <-NULL
  matches2 <-NULL
  DataProteinPosition <- cbind.data.frame(DataProtein$Spectrum, DataProtein$Peptide_mod, DataProtein$All_Proteins, DataProtein$PROTEIN_POS_NUM)
  colnames(DataProteinPosition) <- c("Spectrum", "Peptide_mod", "All_Proteins", "PROTEIN_POS_NUM")
  DataProteinLocation <- cbind.data.frame(DataProtein$Spectrum, DataProtein$PROTEIN_beg, DataProtein$PROTEIN_end)
  colnames(DataProteinLocation) <- c("Spectrum", "PROTEIN_beg", "PROTEIN_end")
  
  DataProteinLocationND <- DataProteinLocation[!duplicated(DataProteinLocation$Spectrum),]
  
  for (i in 1:length(DataProteinPosition$Spectrum)) {
    value <-NULL
   
    for (j in 1:length(DataProteinLocationND$Spectrum)) {
    
      value[j] <- dplyr::between(DataProteinPosition$PROTEIN_POS_NUM[i],DataProteinLocationND$PROTEIN_beg[j], DataProteinLocationND$PROTEIN_end[j])
      
      
    }
    matches[[i]] <-value
    df <- do.call("rbind",matches) #combine all vectors into a matrix
    sumdata<-rowSums(df)
   
  }
  
 
  N_Count = c(N_Count, sumdata)
  
}
counts <- as.data.frame(N_Count)

finalSet <- cbind.data.frame(dataf,counts)

return(finalSet)

}