

# This function takes the Binomial adjusted scores and bring them to peptide-site level by using the maximum (mean(peptidoform)) to choose 
# peptide-site level match # PLEASE NOTE THAT CAN BE DRAWS IN MEAN SCORING BETWEEN DIFFERENT PEPTIDOFORMS AND HENCE MORE THAN 1 SOLUTION FOR 1 SITE

# From Pepidoform site to Peptide site #
########################################


PeptidoformToPeptide_mean <- function(dframe){
  
  
 KeepCol <- c("NewScore3","Peptide_mod", "Peptide", "Peptidoform", "PTM_positions", "PTM_Protein_Positions", "PTM_pos", "PROTEIN_POS_NUM", "PROTEIN_LOC", "PROTEIN_OCC", "Amino")
  
  dframe2 <- dframe[,names(dframe) %in% KeepCol]
  
  dframe2$Peptide_Pos <- paste0( dframe2$Peptide, "_", dframe2$PTM_positions)
  
  dframe2$Peptide_num_Pho <- paste0( dframe2$Peptide, "_",str_count(dframe2$Peptidoform, "Phosphorylation"))
  
  dframe2 <- dframe2 %>% 
    group_by(Peptide_mod) %>% 
    mutate(MeanScore = mean(NewScore3))

Adjusted_PepMean <- dframe2 %>% group_by(Peptide_num_Pho) %>% top_n(1, MeanScore)

Adjusted_PepMean <- Adjusted_PepMean[!duplicated(Adjusted_PepMean), ]

Adjusted_PepMean$MeanScore <- NULL

return(Adjusted_PepMean)
}