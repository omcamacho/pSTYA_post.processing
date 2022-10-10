
# This programme calculates at protein level and number of unique peptides #
############################################################################

Num_Proteins_Unique_Pep <- function(dframeR1,dframeR5){
 
  dframeR1$Amino <- str_sub(substr(dframeR1$Peptide,1,dframeR1$PTM_positions),-1)  
NumProt_R1 <- cbind.data.frame(dframeR1$PROTEIN_LOC, dframeR1$PROTEIN_OCC, dframeR1$Amino)
colnames(NumProt_R1)<- c("Protein_Location","Unique_Peptide","Amino")
NumProt_R1 <- NumProt_R1[!duplicated(NumProt_R1), ]
UniquePep_R1 <- NumProt_R1
UniquePep_R1$Protein_Location <-NULL
UniquePep_R1$Amino <-NULL
UniquePep_R1 <- as.data.frame(UniquePep_R1[!duplicated(UniquePep_R1), ])
colnames(UniquePep_R1)<- c("Unique_Peptide")

prot_R1 <- length(NumProt_R1$Protein_Location)
pep_R1 <- length(UniquePep_R1[,1])

dframeR5$Amino <- str_sub(substr(dframeR5$Peptide,1,dframeR5$PTM_positions),-1) 
NumProt_R5 <- cbind.data.frame(dframeR5$PROTEIN_LOC, dframeR5$PROTEIN_OCC, dframeR5$Amino)
colnames(NumProt_R5)<- c("Protein_Location","Unique_Peptide","Amino")
NumProt_R5 <- NumProt_R5[!duplicated(NumProt_R5), ]
UniquePep_R5 <- NumProt_R5
UniquePep_R5$Protein_Location <-NULL
UniquePep_R5$Amino <-NULL
UniquePep_R5 <- as.data.frame(UniquePep_R5[!duplicated(UniquePep_R5), ])
colnames(UniquePep_R5)<- c("Unique_Peptide")

prot_R5 <- length(NumProt_R5$Protein_Location)
pep_R5 <- length(UniquePep_R5[,1])

metrics <- c(prot_R1, prot_R5, pep_R1, pep_R5)

results <- as.list(c(metrics, NumProt_R1, UniquePep_R1, NumProt_R5, UniquePep_R5))

return(results)
}