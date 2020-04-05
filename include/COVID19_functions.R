# ------------------------------------------------------------------------------------------------------------
# Use to unpivot some data of assumption
# Take a .csv file path of a "square" matrix
# Return the same but unpivoted data
# 'Colnames' should be a 3-sized vector

unpiv_triangl_csv <- function(file, col, header=T, sep='\t', dec=',', encoding='utf-8'){
  up <- read.csv(file, header=header, sep=sep, dec=dec, encoding=encoding)
  up <- melt(up, id.vars = 1, na.rm = T)
  colnames(up) <- col
  up[,2] <- gsub('[A-Z]', '', up[,2])
  up[,2] <- gsub('\\.', '-', up[,2])
  up[,2] <- gsub('80[-\\.]*', '80-++', up[,2])
  i <- up[,col[1]]!=up[,col[2]]
  low <- data.frame(up[i,col[2]], up[i, col[1]], up[i, col[3]])
  colnames(low) <- col
  return(rbind(up, low))
}

# ------------------------------------------------------------------------------------------------------------
# Addition of transmission table
#l <- H[c('household', 'activity')]
# ADD_TRANS_TAB(H[c('household', 'activity', 'household')])
ADD_TRANS_TAB <- function(l, oprt='+'){
  if(length(l)<1){
    return(l[[1]])
  } else {
    col <- colnames(l[[1]])
    stopifnot(identical(col,colnames(l[[2]])))
    l12 <- merge(x=l[[1]], y=l[[2]], by=col[1:2], all.x=F, all.y=F, suffixes = 1:2)
    jc <- paste0(col[3], 1:2)
    if(oprt=='+') l12[,col[3]] <- l12[,jc[1]]+l12[,jc[2]] else l12[,col[3]] <- l12[,jc[1]]+l12[,jc[2]]
    
    if(length(l)==2){
      return(l12[,col])
    } else {
      l <- l[-(1:2)]
      l[[oprt]] <- l12[,col]
      return(ADD_TRANS_TAB(l=l, oprt=oprt))
    }
  }
}




#rm(idg);rm(tm);rm(sr);rm(cdday)
#idg=IDG;tm=H$all;sr=H$sr;cdday=H$cdday
#idg=LAST_DAY;tm=H$all;sr=H$sr;cdday=H$cdday
#idg=S[S$dt==7,];tm=H$all;sr=H$sr;cdday=H$cdday
#LAST_DAY=EVOL_LAST_DAY(idg, H$all, H$sr, H$cdday)
#LAST_DAY=EVOL_LAST_DAY(IDG, H$household, H$sr, H$cdday)
#LAST_DAY=EVOL_LAST_DAY(LAST_DAY, H$household, H$sr, H$cdday)
#S=SIMU_PERIOD(n=3, idg=IDG, tm=H$household, sr=H$sr, cdday=H$cdday)
EVOL_LAST_DAY <- function(idg, tm, sr, cdday){
  # Test
  #if(idg$dt[1]%%10==0){
  #  print(
  #    paste0(
  #      'Start calculating day :', idg$dt[1], '. IDG has ', dim(idg)[1], ' rows and ',
  #      sum(duplicated(idg[, c('liv', 'age', 'dday', 'dt')])),' duplicated records. Sum of people :', sum(idg$qt)))
  #}
  col <- colnames(idg)[1:5]
  
  #------------------------------------
  # 1 - New death
  #------------------------------------
  
  # 1.1:  Report of the global survival risk by age
  idg <- merge(idg, sr, on='age', all.x=T, all.y=F)
  idg$sr[is.na(idg$sr)] <- 0
  # 1.2:  Report of the transmission and death coefficient by disease date
  idg <- merge(idg, cdday, on=c('liv', 'dday'), all.x=T, al.y=F)
  idg$ctrans[is.na(idg$ctrans)] <- 0
  idg$pdeath[is.na(idg$pdeath)] <- 0
  # 1.3:  Death quantity is a binomial random draw
  idg$new_death <- apply(X=idg[, c('qt', 'sr', 'pdeath')], MARGIN=1, function(r) rbinom(1,r[[1]],r[[2]]*r[[3]]))
  # 1.4: Update of demography
  idg$new_qt <- idg$qt-idg$new_death
  
  #------------------------------------
  # 2 - New infected
  #------------------------------------
  
  # 2.1:  Report of the transmissibility coefficient by age
  idg2 <-  merge(x=idg[,c('liv', 'age', 'dday', 'ctrans', 'qt')], y=tm, by.x='age', by.y='agent', all.x=T, all.y=F)
  idg2$ntrans[is.na(idg2$ntrans)] <- 0
  # 2.3:  Infected quantity is a Poisson random draw 
  idg2$new_inf <- idg2$qt*sapply(X=idg2$ntrans*idg2$ctrans, FUN=rpois, n=1)
  # 2.4:  Aggregation by age of target people
  idg2 <- data.frame(dday=0, aggregate(new_inf~liv+target, data=idg2, FUN=sum))
  # 2.5:  Report in the main table of the random number of infested people
  idg <- merge(
    x=idg, y=idg2,
    by.x=c('liv', 'age', 'dday'), by.y=c('liv', 'target', 'dday'),
    all.x=T, all.y=F
  )
  idg$new_inf[is.na(idg$new_inf)] <- 0
  # 2.6:  Final infected quantity are min between the draw and the recevable population
  idg$new_inf <- pmin(idg$new_inf, idg$new_qt)
  # 2.7:  Update of demography
  idg$new_qt <- idg$new_qt-idg$new_inf
  
  #------------------------------------
  # 3 - Infected takes a day
  #------------------------------------
  i <- idg$dday>0
  if(sum(i)>0) idg$dday[i] <- idg$dday[i]+1
    
  #------------------------------------
  # 3 - Final return
  #------------------------------------
  out <- list(SAME=idg[, c('liv', 'age', 'new_qt', 'dday', 'dt')])
  out[['DEAD']] <- data.frame(liv=F, aggregate(new_death~age, data=idg, FUN=sum), dday=1, idg$dt[1])
  out[['INFECTED']] <- data.frame(aggregate(new_inf~liv+age, data=idg, FUN=sum), dday=1, idg$dt[1])
  for(i in 1:length(out)) colnames(out[[i]]) <- col

  out <- do.call(rbind, out)
  rownames(out) <- c()
  # Going to the newt day
  out$dt <- out$dt+1
  return(out[out$qt>0,])
}

SIMU_PERIOD <- function(n, idg, tm, sr, cdday){
  LAST_DAY <- EVOL_LAST_DAY(idg=idg, tm=tm, sr=sr, cdday=cdday)
  if(n==1) return(LAST_DAY) else return(rbind(SIMU_PERIOD(n=n-1, idg=LAST_DAY, tm=tm, sr=sr, cdday=cdday), LAST_DAY))
}