# ------------------------------------------------------------------------------------------------------------
# Get a 'sce' data.frame with two date columns inside called 'dtstart', 'dtend'
# Get two date limites 'dtstart' and 'dtend'
# Return the scenario timeline between this limite

GET_SCENAR_TIME <- function(sce, dtstart, dtend){
  i <- is.na(sce$dtstart)
  if(sum(i)>0) sce$dtstart[i] <- pmin(dtstart, sce$dtend[i], na.rm = T)
  i <- is.na(sce$dtend)
  if(sum(i)>0) sce$dtend[i] <- pmax(dtend, sce$dtstart[i], na.rm=T)
  i <- sce$dtstart>dtend|sce$dtend<dtstart
  if(sum(i)>0) sce <- sce[!i,]
  sce$dtstart <- pmax(sce$dtstart, dtstart)
  sce$dtend <- pmin(sce$dtend, dtend)
  out <- seq(from=sce$dtstart[1], to=sce$dtend[1], by='day')
  out <- data.frame(dt=out, sce[1, -(1:2)], row.names=1:length(out))
  i <- format(out$dt, '%u') %in% 1:5
  # Coefficient 0 to working days if specified
  if(!(sce$flgwd[1]) & sum(i)>0) out$coeff[i] <- 0
  # Coefficient 0 to week-end days if specified
  if(!sce$flgwe[1] & sum(!i)>0) out$coeff[!i] <- 0
  if(dim(sce)[1]>1) return(rbind(out, GET_SCENAR_TIME(sce=sce[-1,], dtstart=dtstart, dtend=dtend))) else return(out)
}


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
# Take a small data.fame 'hp' containing parameters of the distribution
# Return a data.fram with daily risk of death and transmission coefficient

GET_DDAY<- function(pdday){
  out <- dweibull(0:99+0.5, shape=pdday$shape[1], scale=pdday$scale[1])
  i <- which(out>10^(-9))
  out <- out[i]
  out <- data.frame(liv=T, dday=i, pdeath=out/sum(out))
  out$ctrans <- pmin(
    c(tail(cumsum(out$pdeath), -pdday$lag[1]), rep(1, pdday$lag[1])),
    c(rep(1,pdday$lag[1]), 1-head(cumsum(out$pdeath), -pdday$lag[1]))
  )
  return(rbind(data.frame(liv=T, dday=0, pdeath=0, ctrans=0), out)) 
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

# Take a data.frame 'data' with cumulated quantity over a date column 'dtcol'
# Return the same data.frame but with new quantities each date

UNCUMUL <- function(data, dtcol, qtcol, bycol=NULL){
  n <- dim(data)[1]
  if(is.null(bycol)){
    data <- data[order(data[,c(dtcol,bycol)]),]
    data[2:n, qtcol] <- data[2:n, qtcol]-data[1:(n-1), qtcol]
  } else {
    data <- do.call(rbind, by(data = data, INDICES = data[,bycol], FUN=UNCUMUL, dtcol=dtcol, qtcol=qtcol))
    rownames(data) <- 1:n
  }
  return(data)
}

# Take a demography and an initialisation data.Frames
# Return a 'IDG' (INfecto-DemoGraphy) data.frame
INIT_IDG <- function(demo, init){
  idg <- cbind(liv=T, demo, dday=0, dt=init$dt[1])
  idg <- merge(x=idg, y=init, by='age', all.x=T, all.y=F, suffixe=c('', '_init'))
  idg$qt_init[is.na(idg$qt_init)] <- 0
  idg$qt <- idg$qt-idg$qt_init
  idg <- rbind(idg[, 1:5],H$init[,colnames(idg)[1:5]])
  return(idg)
}

#rm(idg);rm(tm);rm(sr);rm(pdday)
#idg=IDG;tm=H$all;sr=H$sr;pdday=H$pdday
#idg=LAST_DAY;tm=H$all;sr=H$sr;pdday=H$pdday
#idg=S[S$dt==7,];tm=H$all;sr=H$sr;pdday=H$pdday
#LAST_DAY=EVOL_LAST_DAY(IDG, H$all, H$sr, H$pdday)
#LAST_DAY=EVOL_LAST_DAY(LAST_DAY, H$all, H$sr, H$pdday)
#S=SIMU_PERIOD(n=30, idg=IDG, tm=H$all, sr=H$sr, pdday=H$pdday)
EVOL_LAST_DAY <- function(idg, tm, sr, pdday){
  # Test
  #if(idg$dt[1]%%10==0){
  #  print(
  #    paste0(
  #      'Start calculating day :', idg$dt[1], '. IDG has ', dim(idg)[1], ' rows and ',
  #      sum(duplicated(idg[, c('liv', 'age', 'dday', 'dt')])),' duplicated records. Sum of people :', sum(idg$qt)))
  #}
  col <- colnames(idg)[1:5]
  cdday <- GET_DDAY(pdday=pdday)
  
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
  # 3 - Infected people take a day
  #------------------------------------
  i <- idg$dday>0
  if(sum(i)>0) idg$dday[i] <- idg$dday[i]+1
    
  #------------------------------------
  # 3 - Final return
  #------------------------------------
  out <- list(
    qt=idg[,colnames(idg)!='qt'],
    death=data.frame(liv=F, dday=1, dt=idg$dt[1], aggregate(new_death~age, data=idg, FUN=sum)),
    inf=data.frame(dday=1, dt=idg$dt[1], aggregate(new_inf~liv+age, data=idg, FUN=sum))
  )
  for(i in c('qt', 'death', 'inf')){
    colnames(out[[i]])[colnames(out[[i]])==paste0('new_', i)] <- 'qt'
    out[[i]] <- out[[i]][,col]
  }
  
  out <- do.call(rbind, out)
  rownames(out) <- c()
  # Going to the newt day
  out$dt <- out$dt+1
  return(out[out$qt>0,])
}

# Recursive function to calculate EVOL_LAST_DAY from 'idg' to 'dt_end' (range of date)

SIMU_PERIOD <- function(dt_end, idg, tm, sr, pdday){
  if(dt_end<=idg$dt[1]){
    stop('You need to provide a further date of end.')
  } else {
    LAST_DAY <- EVOL_LAST_DAY(idg=idg, tm=tm, sr=sr, pdday=pdday)
    if(LAST_DAY$dt[1]==dt_end){
      return(LAST_DAY)
    } else {
      return(rbind(SIMU_PERIOD(dt_end=dt_end, idg=LAST_DAY, tm=tm, sr=sr, pdday=pdday), LAST_DAY))
    }
  }
}