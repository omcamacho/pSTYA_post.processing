
# Calculating FLR #
###################

FLR_Prod_Adj <- function(dataf) {

# Count total STYs and As #
  ReducedTable2 <-  dataf
ReducedTable2$num_PTM_Y <- str_count(ReducedTable2$Peptide, "Y")
ReducedTable2$num_PTM_S <- str_count(ReducedTable2$Peptide, "S")
ReducedTable2$num_PTM_T <- str_count(ReducedTable2$Peptide, "T")
ReducedTable2$num_PTM_A <- str_count(ReducedTable2$Peptide, "A")

ratio <- sum(ReducedTable2$num_PTM_Y+ReducedTable2$num_PTM_S+ReducedTable2$num_PTM_T)/sum(ReducedTable2$num_PTM_A)


ReducedTable2$Amino <- str_sub(substr(ReducedTable2$Peptide,1,ReducedTable2$PTM_positions),-1)


ReducedTable2$value <- ifelse( ReducedTable2$Amino == "A", 1, 0)

ReducedTable2$PepwithA<- as.character(gregexpr(pattern ='A\\[Phospho',ReducedTable2$Peptide_mod))

ReducedTable2$AsN <-ifelse(ReducedTable2$PepwithA=='-1',0,1)



# Order by Score 3#

ReducedTable2 <- ReducedTable2[order(ReducedTable2$Peptide),]

ReducedTable2sorted3 <- ReducedTable2[order(-ReducedTable2$Prod_Adjusted_Score),]

ReducedTable2sorted3$cums <- cumsum(ReducedTable2sorted3$value)

ReducedTable2sorted3$cums2 <- cumsum(ReducedTable2sorted3$AsN)

ReducedTable2sorted3$Rnumber <- seq.int(nrow(ReducedTable2sorted3))

ReducedTable2sorted3 <- transform(ReducedTable2sorted3, FLR_Adj_Score = ((ratio*cums2)+cums2)/Rnumber)


return(ReducedTable2sorted3)
}

