
# Binomial adjustment function for Arabidopsis data set TPP pipeline #
######################################################################

binAdjustAra <- function(dataf){
  PTMscore <- dataf 
  PTMscore$Peptidoform <- paste0(PTMscore$Peptide_mod,"_",PTMscore$PTM_positions)
  PTMscore$Amino <- str_sub(substr(PTMscore$Peptide,1,PTMscore$PTM_positions),-1)
  head(PTMscore)
  
  
  
  PTMscore$PTM_pos<-paste0(as.character(PTMscore$Peptide),"_",as.character(PTMscore$PTM_positions))
  
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
  
  PTMscore$PROTEIN_POS<-lapply(1:nrow(X1), function(i) Z1$Vector[[i]][X1$Vector[[i]]])
  
  
  
  PTMscore$PROTEIN_POS_NUM <-as.numeric(PTMscore$PROTEIN_POS)
  
  
  PTMscoreSorted <- PTMscore[order((-PTMscore$PROTEIN_POS_NUM)),]
  PTMscoreSorted <- PTMscoreSorted[order(PTMscoreSorted$Single_Protein),]
  
  PTMscoreSorted$PTM_Beginning <- NULL
  PTMscoreSorted$PTM_End <- NULL
  
  BEG<-ldply(PTMscoreSorted$PTM_beg2,function(s){t(data.frame(unlist(s)))})
  PTMscoreSorted$PTM_Beginning<-as.vector(BEG[,1])
  PTMscoreSorted["PTM_Beginning"][PTMscoreSorted["PTM_Beginning"] == 0] <- 1
  PTMscoreSorted$PTM_Lth<-nchar(PTMscoreSorted$Peptide)
  
  PTMscoreSorted$PROTEIN_beg<-PTMscoreSorted$PROTEIN_POS_NUM-PTMscoreSorted$PTM_Beginning+1
  PTMscoreSorted$PROTEIN_end<-PTMscoreSorted$PROTEIN_beg+PTMscoreSorted$PTM_Lth-1
  
  PTMscoreSorted <- 
    PTMscoreSorted %>%
    group_by(Single_Protein) %>%
    mutate(lag.Beg = dplyr::lag(PROTEIN_beg, n = 1, default = NA))
  PTMscoreSorted <- 
    PTMscoreSorted %>%
    group_by(Single_Protein) %>%
    mutate(lag.End = dplyr::lag(PROTEIN_end, n = 1, default = NA))
  
  PTMscoreSorted<-
    PTMscoreSorted %>% 
    mutate(lag.Beg = coalesce(lag.Beg,PROTEIN_beg))
  PTMscoreSorted<-
    PTMscoreSorted %>% 
    mutate(lag.End = coalesce(lag.End,PROTEIN_end))
  
  
  
  PTMscoreSorted$CAT <- ifelse(((PTMscoreSorted$PROTEIN_POS_NUM >= PTMscoreSorted$lag.Beg)&(PTMscoreSorted$PROTEIN_POS_NUM <= PTMscoreSorted$lag.End))|(PTMscoreSorted$PROTEIN_POS_NUM<=2), "Other", "Start")
  
  
  
  PTMscoreSorted<-PTMscoreSorted %>%
    group_by(Single_Protein) %>%
    mutate(grp = cumsum(CAT == "Start")+1)
  
  PTMscoreSorted$PROTEIN_LOC <-paste0(PTMscoreSorted$Single_Protein,"_",PTMscoreSorted$PROTEIN_POS_NUM)
  PTMscoreSorted$PROTEIN_OCC <-paste0(PTMscoreSorted$Single_Protein,"_",PTMscoreSorted$grp)
  
  
  # Attempt to use binomial distribution to account for number of matches found during search #
  
  PTM2<- PTMscoreSorted %>% group_by(PROTEIN_LOC) %>% dplyr::mutate(count_S = n())
  
  
  
  # Number of different possibilities to be observed N #
  
  PTM2Unique <- cbind.data.frame(PTM2$Spectrum,PTM2$PROTEIN_OCC,PTM2$Peptide_mod)
  
  colnames(PTM2Unique)<-c("Spectrum","PROTEIN_OCC","Peptide_mod")
  
  PTM2Unique <- PTM2Unique %>% dplyr::arrange(Spectrum) %>% distinct(Spectrum, .keep_all = TRUE) %>% ungroup()
  
  
  PTM2Unique<- PTM2Unique %>% group_by(PROTEIN_OCC) %>% dplyr::mutate(count_ALL = n())
  
  PTM2Unique$SpectA<- as.character(gregexpr(pattern ='A\\[Phospho',PTM2Unique$Peptide_mod))
  
  PTM2Unique$Num_As <-ifelse( PTM2Unique$SpectA=='-1',0,1)
  
  PTM2Unique$PROTEIN_OCC <- NULL
  
  PTM3 <- merge(PTM2Unique,PTM2, by="Spectrum")
  
  
  
  PTM3$probChance3 <-sum(PTM2Unique$Num_As)/length(PTM2Unique$Spectrum)
  
  
  # Calculation of binomial prob #
  
  maxCount<-max(PTM3$count_ALL)
  
  PTM3$bioprob3 <- dbinom(PTM3$count_S, PTM3$count_ALL, PTM3$probChance3)
  
  # Combine as independent p1*(1-P2) #
  
  
  PTM3$NewScore3 <- (1-PTM3$bioprob3)*PTM3$PTM_final_prob
  
  # Multiplying to the PMS Score #
  
  
  
  

  
return(PTM3)
}