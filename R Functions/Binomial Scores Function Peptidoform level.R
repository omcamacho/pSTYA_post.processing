
binAdjustPform <- function(dataf){
  PTMscore <- dataf
  PTMscore$Peptidoform <- paste0(PTMscore$Peptide_mod,"_",PTMscore$PTM.positions)
  
  
  PTMscore$PTM_pos<-paste0(as.character(PTMscore$Peptide),"_",as.character(PTMscore$PTM.positions))
  
  PTMscore$Single_Protein<-word(PTMscore$All_Proteins,1,1,":")
  PTMscore$PRO_pos<-word(as.character(PTMscore$All_PTM_protein_positions),1,1,":")
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
  X1$col6 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+5
  X1$col7 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+6
  X1$col8 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+7
  X1$col9 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+8
  X1$col10 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+X1[10]+9
  X1$col11 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+X1[10]+X1[11]+10
  X1$col12 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+X1[10]+X1[11]+X1[12]+11
  X1$col13 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+X1[10]+X1[11]+X1[12]+X1[13]+12
  X1$col14 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+X1[10]+X1[11]+X1[12]+X1[13]+X1[14]+13
  X1$col15 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+X1[10]+X1[11]+X1[12]+X1[13]+X1[14]+X1[15]+14
  X1$col16 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+X1[10]+X1[11]+X1[12]+X1[13]+X1[14]+X1[15]+X1[16]+15
  X1$col17 <- X1[1]+X1[2]+X1[3]+X1[4]+X1[5]+X1[6]+X1[7]+X1[8]+X1[9]+X1[10]+X1[11]+X1[12]+X1[13]+X1[14]+X1[15]+X1[16]++X1[17]+16
  Y1<-ldply(PTMscore$PTM_end2,function(s){t(data.frame(unlist(s)))})
  Y1$col1 <- Y1[1]
  Y1$col2 <- Y1[2]
  Y1$col3 <- Y1[3]
  Y1$col4 <- Y1[4]
  Y1$col5 <- Y1[5]
  Y1$col6 <- Y1[6]
  Y1$col7 <- Y1[7]
  Y1$col8 <- Y1[8]
  Y1$col9 <- Y1[9]
  Y1$col10 <- Y1[10]
  Y1$col11 <- Y1[11]
  Y1$col12 <- Y1[12]
  Y1$col13 <- Y1[13]
  Y1$col14 <- Y1[14]
  Y1$col15 <- Y1[15]
  Y1$col16 <- Y1[16]
  Y1$col17 <- Y1[17]
  
  
  X1$Vector<-lapply(1:nrow(X1), function(i) c(is.element(Y1$col1[i,]-X1$col1[i,],PTMscore$PTM.positions[i]),is.element(Y1$col2[i,]-X1$col2[i,],
                                                                                                                       PTMscore$PTM.positions[i])
                                              ,is.element(Y1$col3[i,]-X1$col3[i,],PTMscore$PTM.positions[i]),is.element(Y1$col4[i,]-X1$col4[i,],PTMscore$PTM.positions[i]),
                                              is.element(Y1$col5[i,]-X1$col5[i,],PTMscore$PTM.positions[i]),is.element(Y1$col6[i,]-X1$col6[i,],PTMscore$PTM.positions[i])
                                              ,is.element(Y1$col7[i,]-X1$col7[i,],PTMscore$PTM.positions[i]),is.element(Y1$col8[i,]-X1$col8[i,],PTMscore$PTM.positions[i]),
                                              is.element(Y1$col9[i,]-X1$col9[i,],PTMscore$PTM.positions[i]),is.element(Y1$col10[i,]-X1$col10[i,],
                                                                                                                       PTMscore$PTM.positions[i])
                                              ,is.element(Y1$col11[i,]-X1$col11[i,],PTMscore$PTM.positions[i]),is.element(Y1$col12[i,]-X1$col12[i,],PTMscore$PTM.positions[i]),
                                              is.element(Y1$col13[i,]-X1$col13[i,],PTMscore$PTM.positions[i]),is.element(Y1$col14[i,]-X1$col14[i,],PTMscore$PTM.positions[i])
                                              ,is.element(Y1$col15[i,]-X1$col15[i,],PTMscore$PTM.positions[i]),is.element(Y1$col16[i,]-X1$col16[i,],PTMscore$PTM.positions[i])
                                              ,is.element(Y1$col17[i,]-X1$col17[i,],PTMscore$PTM.positions[i]))) 
  
  PTMscore$PRO_pos_list <- lapply(1:nrow(X1), function(i) as.list(el(strsplit(PTMscore$PRO_pos[i], ";"))))
  
  Z1<-ldply(PTMscore$PRO_pos_list,function(s){t(data.frame(unlist(s)))})
  
  Z1$Vector <- lapply(1:nrow(X1), function(i) c(Z1[i,1],Z1[i,2],Z1[i,3],Z1[i,4],Z1[i,5],Z1[i,6],Z1[i,7],Z1[i,8],Z1[i,9],Z1[i,10],Z1[i,11],Z1[i,12],Z1[i,13],Z1[i,14],Z1[i,15],Z1[i,16],Z1[i,17]))
  
  
  PTMscore$PROTEIN_POS<-lapply(1:nrow(X1), function(i) Z1$Vector[[i]][X1$Vector[[i]]])
  
  
  
  PTMscore$PROTEIN_POS_NUM <-as.numeric(PTMscore$PROTEIN_POS)
  
  
  PTMscoreSorted <- PTMscore[order((-PTMscore$PROTEIN_POS_NUM)),]
  PTMscoreSorted <- PTMscoreSorted[order(PTMscoreSorted$Single_Protein),]
  
  PTMscoreSorted$PTM_Beginning <- NULL
  PTMscoreSorted$PTM_End <- NULL
  
  PTMscoreSorted$PTM_Lth<-nchar(PTMscoreSorted$Peptide)
  PTMscoreSorted$PTM_right<-PTMscoreSorted$PTM_Lth-PTMscoreSorted$PTM.positions
  
  PTMscoreSorted$PROTEIN_end<-PTMscoreSorted$PROTEIN_POS_NUM+PTMscoreSorted$PTM_right
  PTMscoreSorted$PROTEIN_beg<-PTMscoreSorted$PROTEIN_end-PTMscoreSorted$PTM_Lth+1
  
  PTMscoreSorted$PROTEIN_LOC <-paste0(PTMscoreSorted$Single_Protein,"_",PTMscoreSorted$PROTEIN_POS_NUM)
  
  
  # Attempt to use binomial distribution to account for number of matches found during search #
  
  PTM2<- PTMscoreSorted %>% group_by(PROTEIN_LOC) %>% dplyr::mutate(count_S = n())
  
  
  
  # Number of different possibilities to be observed N #
  
  PTM3 <- FFreq(PTM2)
  
  PTM3$SpectA<- as.character(gregexpr(pattern ='A\\[Phospho',PTM3$Peptide_mod))
  
  PTM3$Num_As <-ifelse( PTM3$SpectA=='-1',0,1)
  
  PTM2Unique <- PTM3 %>% dplyr::arrange(Spectrum) %>% distinct(Spectrum, .keep_all = TRUE) %>% ungroup()
  
  
  
  PTM3$probChance3 <-sum(PTM2Unique$Num_As)/length(PTM2Unique$Spectrum)
  
  
  # Calculation of binomial prob #
  
  
  PTM3$bioprob3 <- dbinom(PTM3$count_S, PTM3$N_Count, PTM3$probChance3)
  
  # Combine as independent p1*(1-P2) #
  
  
  PTM3$Bin_Adjusted_Score <- (1-PTM3$bioprob3)*PTM3$PTM_final_prob
  
  
  BinomialSet3 <- aggregate(Bin_Adjusted_Score ~ Peptidoform+PROTEIN_LOC, data = PTM3, max)
  
  
  total <- merge(PTMscore,BinomialSet3,by="Peptidoform")
  
  
  ReducedTable2 <- total %>% dplyr::arrange(Peptidoform) %>% distinct(Peptidoform, .keep_all = TRUE) %>% ungroup()
  

return(ReducedTable2)
}