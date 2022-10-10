

# This function takes the Binomial adjusted scores and bring them to peptide-site level by using the maximum (mean(peptidoform)) to choose 
# peptide-site level match # PLEASE NOTE THAT CAN BE DRAWS IN MEAN SCORING BETWEEN DIFFERENT PEPTIDOFORMS AND HENCE MORE THAN 1 SOLUTION FOR 1 SITE

# From Pepidoform site to Peptide site #
########################################


PeptidoformToPeptide_max <- function(dframe){
  
  
  dframe$Peptide_num_Pho <- paste0( dframe$Peptide, "_",str_count(dframe$Peptidoform, "Phosphorylation"))
  
  Adjusted_PepMax <- dframe %>% group_by(Peptide_num_Pho) %>% top_n(1, NewScore3)
  
  Adjusted1 <- as.data.frame(cbind( Adjusted_PepMax$Peptide_mod))
  
  
  colnames(Adjusted1)<-c("Peptide_mod")
  
  
 KeepCol <- c("NewScore3","Peptide_mod", "Peptide", "Peptidoform", "PTM_positions", "PTM_Protein_Positions", "PTM_pos", "PROTEIN_POS_NUM", "PROTEIN_LOC", "PROTEIN_OCC", "Amino")
  
  dframe2 <- dframe[,names(dframe) %in% KeepCol]
  
  dframe2$Peptide_Pos <- paste0( dframe2$Peptide, "_", dframe2$PTM_positions)
  
  dframe2$Peptide_num_Pho <- paste0( dframe2$Peptide, "_",str_count(dframe2$Peptidoform, "Phosphorylation"))
  
  datset2 <- merge(x=Adjusted1, y=dframe2, by="Peptide_mod", all.x = TRUE)

  Adjusted_PepMax <- Adjusted_PepMax[!duplicated(Adjusted_PepMax), ]


return(datset2)
}