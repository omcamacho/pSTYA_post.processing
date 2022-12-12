
binAdjustSyn <- function(dataf){
  PTMscore <- dataf
  
  
  PTMscore$PTM_pos<-paste0(as.character(PTMscore$Peptide),"_",as.character(PTMscore$PTM_positions))
  
  PTMscore$Single_Protein<-word(PTMscore$Protein,1,1,":")
  PTMscore$PRO_pos<-word(as.character(PTMscore$PTM_Protein_Positions),1,1,":")
  # PTMscore$PRO_pos<-str_remove(PTMscore$PRO_pos, "0;")
  PTMscore$Peptide_mod2 <-str_remove(PTMscore$Peptide_mod, "\\-")
  
  PTMscore$PTM_Beginning <- gregexpr(pattern ='\\[',PTMscore$Peptide_mod2)
  PTMscore$PTM_End <- gregexpr(pattern ='\\]',PTMscore$Peptide_mod2)
  
  PTMscore$PTM_end2 <- lapply(1:nrow(PTMscore), function(i) unlist(PTMscore$PTM_End[i])-1)
  PTMscore$PTM_beg2 <- lapply(1:nrow(PTMscore), function(i) unlist(PTMscore$PTM_Beginning[i]) - 1)
  PTMscore$PTM_length <- lapply(1:nrow(PTMscore), function(i)  unlist(PTMscore$PTM_End[i]) -  unlist(PTMscore$PTM_Beginning[i]))
  
  X1<-ldply(PTMscore$PTM_length,function(s){t(data.frame(unlist(s)))})
  X1$col1 <- X1[1]
  X1$col2 <- X1[1]+X1[2]+1
  X1$col3 <- X1[1]+X1[2]+X1[3]+2
  X1$col4 <- X1[1]+X1[2]+X1[3]+X1[4]+3
  X1$col5 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+4
  Y1<-ldply(PTMscore$PTM_end2,function(s){t(data.frame(unlist(s)))})
  Y1$col1 <- Y1[1]
  Y1$col2 <- Y1[2]
  Y1$col3 <- Y1[3]
  Y1$col4 <- Y1[4]
  Y1$col5 <- Y1[5]
  
  X1$Vector<-lapply(1:nrow(X1), function(i) c(is.element(Y1$col1[i,]-X1$col1[i,],PTMscore$PTM_positions[i]),is.element(Y1$col2[i,]-X1$col2[i,],
                                                                                                                       PTMscore$PTM_positions[i])
                                              ,is.element(Y1$col3[i,]-X1$col3[i,],PTMscore$PTM_positions[i]),is.element(Y1$col4[i,]-X1$col4[i,],PTMscore$PTM_positions[i]),
                                              is.element(Y1$col5[i,]-X1$col5[i,],PTMscore$PTM_positions[i]))) 
  
  PTMscore$PRO_pos_list <- lapply(1:nrow(X1), function(i) as.list(el(strsplit(PTMscore$PRO_pos[i], ";"))))
  
  Z1<-ldply(PTMscore$PRO_pos_list,function(s){t(data.frame(unlist(s)))})
  
  Z1$Vector <- lapply(1:nrow(X1), function(i) c(Z1[i,1],Z1[i,2],Z1[i,3],Z1[i,4],Z1[i,5]))
  
  PTMscore$PROTEIN_POS<- PTMscore$PRO_pos
  
  
  
  PTMscore$PROTEIN_POS_NUM <-as.numeric(PTMscore$PROTEIN_POS)
  
  
  PTMscoreSorted <- PTMscore[order((-PTMscore$PROTEIN_POS_NUM)),]
  PTMscoreSorted <- PTMscoreSorted[order(PTMscoreSorted$Single_Protein),]
  
  PTMscoreSorted$PTM_Beginning <- NULL
  PTMscoreSorted$PTM_End <- NULL
  
  PTMscoreSorted$PTM_Lth<-nchar(PTMscoreSorted$Peptide)
  PTMscoreSorted$PTM_right<-PTMscoreSorted$PTM_Lth-PTMscoreSorted$PTM_positions
  
  PTMscoreSorted$PROTEIN_end<-PTMscoreSorted$PTM_Protein_Positions+PTMscoreSorted$PTM_right
  PTMscoreSorted$PROTEIN_beg<-PTMscoreSorted$PROTEIN_end-PTMscoreSorted$PTM_Lth+1
  
  
  PTMscoreSorted$PROTEIN_LOC <-paste0(PTMscoreSorted$Single_Protein,"_",PTMscoreSorted$PROTEIN_POS_NUM)

  
  # Attempt to use binomial distribution to account for number of matches found during search #
  
  PTM2<- PTMscoreSorted %>% group_by(PROTEIN_LOC) %>% dplyr::mutate(count_S = n())
  

  # Number of different possibilities to be observed N #
  
  PTM3 <- FFreq(PTM2)
  
  PTM3$SpectA<- as.character(gregexpr(pattern ='A\\[Phospho',PTM3$Peptide_mod))
  
  PTM3$Num_As <-ifelse( PTM3$SpectA=='-1',0,1)
  
  PTM2Unique <- PTM3 %>% dplyr::arrange(All_USI) %>% distinct(All_USI, .keep_all = TRUE) %>% ungroup()
  
  
  PTM3$probChance3 <-sum(PTM2Unique$Num_As)/length(PTM2Unique$All_USI)
  
  
  # Calculation of binomial prob #
  
  
  PTM3$bioprob3 <- dbinom(PTM3$count_S, PTM3$N_Count, PTM3$probChance3)
  
  # Combine as independent p1*(1-P2) #
  
  PTM3$PTM_final_prob<-PTM3$PTM_score*PTM3$Score
  
  PTM3$NewScore3 <- (1-PTM3$bioprob3)*PTM3$PTM_final_prob
  

return(PTM3)
}