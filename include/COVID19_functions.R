# Try: NEXT_DAY_SIR(S=0.999, I=0.001, beta=0.2, lambda=21, dt=as.Date('2020-02-15'), n.max=3)
NEXT_DAY_SIR <- function(S, I, beta, lambda, dt0, n.max=1, R.max=0){
  if((S+I)>1|(S<0)|(I<0)) stop('S and I are to be parameters in [0,1] with sum <1')
  if((beta<0)|(lambda<0)) stop('Beta and lambda are to be positive parameters')
  # New healthy people
  nS <- S-beta*S*I
  # New infected people
  nI <- I+beta*S*I-I/lambda
  # New cured people
  nR <- (1-S-I)+I/lambda
  out <- rbind(data.frame(dt=dt0, type=c('S', 'I', 'R'), p=c(nS, nI, nR)))
  if(n.max<=1|nR>R.max) return(out) else return(rbind(NEXT_DAY_SIR(S=nS, I=nI, beta=beta, lambda=lambda, dt0=dt0+1, n.max=n.max-1, R.max=R.max), out))
}


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

UNPIV_TRIANGL_CSV <- function(file, col, header=T, sep='\t', dec=',', encoding='utf-8'){
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
# Try: GET_DDAY(H$pdday, H$beta)
# pdday=H$pdday;beta=H$beta;l=60
GET_DDAY<- function(pdday, beta, l=60, lockdown=F){
  out <- dweibull(1:l-0.5, shape=pdday$shape[1], scale=pdday$scale[1])
  out <- data.frame(liv=T, dday=1:l, pdeath=out/sum(out))
  out$ctrans <- pmin(outer(out$pdeath, beta$power, '^')%*%beta$value,1)
  if(lockdown) out <- out[out$pdeath<0.05,]
  return(rbind(data.frame(liv=T, dday=0, pdeath=0, ctrans=0), out)) 
}
# out <- GET_DDAY(H$pdday, H$beta);plot(data=out, ctrans~dday, type='l', col='green');lines(data=out, pdeath~dday, col='blue')

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
# idg=INIT_IDG(demo=D$dFR, init=H$init); pdday=H$pdday; beta=H$beta; sr=H$sr; randommode=F
# idg=LAST_DAY;pdday=H$pdday; beta=H$beta; sr=H$sr; randommode=F
# LAST_DAY=EVOL_LAST_DAY(INIT_IDG(demo=D$dFR, init=H$init), H$pdday, H$beta, H$sr)
# LAST_DAY=EVOL_LAST_DAY(LAST_DAY, H$pdday, H$beta, H$sr)
#S=SIMU_PERIOD(n=30, idg=IDG, tm=H$all, sr=H$sr, pdday=H$pdday)
EVOL_LAST_DAY <- function(idg, pdday, beta, sr, randommode=F){
  # Test
  #if(as.integer(format(idg$dt[1], '%d'))%%10==0){
  #  print(
  #    paste0(
  #      'Start calculating day :', idg$dt[1], '. IDG has ', dim(idg)[1], ' rows and ',
  #      sum(duplicated(idg[, c('liv', 'age', 'dday', 'dt')])),' duplicated records. Sum of people :', sum(idg$qt)))
  #}
  col <- colnames(idg)[1:5]
  cdday <- GET_DDAY(pdday=pdday, beta=beta)
  
  #------------------------------------
  # 1 - New death
  #------------------------------------
  
  # Is there living people ?
  i <- idg$liv
  if(sum(i)>0){
    # 1.1:  Report of the global survival risk by age
    idg <- merge(idg, sr, on='age', all.x=T, all.y=F)
    idg$sr[is.na(idg$sr)] <- 0
    # 1.2:  Report of the transmission and death coefficient by disease date
    idg <- merge(idg, cdday, on=c('liv', 'dday'), all.x=T, al.y=F)
    idg$ctrans[is.na(idg$ctrans)] <- 0
    idg$pdeath[is.na(idg$pdeath)] <- 0
    # 1.3:  Death quantity is a binomial random draw
    #idg$new_death <- apply(X=idg[, c('qt', 'sr', 'pdeath')], MARGIN=1, function(r) rbinom(1,r[[1]],r[[2]]*r[[3]]))
    if(randommode){
      idg$new_death <- apply(X=idg[, c('qt', 'sr', 'pdeath')], MARGIN=1, function(r) rbinom(1,r[[1]],r[[2]]*r[[3]]))
    } else {
      idg$new_death <- idg$qt*idg$sr*idg$pdeath
    }
    # 1.4: Update of demography
    idg$new_qt <- idg$qt-idg$new_death
    
    #------------------------------------
    # 2 - New infected
    #------------------------------------
    
    # Is there infectable people ?
    i <- idg$liv&idg$dday<=0
    if(sum(i)>0){
      # 2.1:  Cartesian product. For each infected age, how many people can be infected ?
      idg2 <-  merge(
        x=cbind(j=T, idg[idg$liv&idg$dday>0, c('qt', 'ctrans')]),
        y=cbind(j=T, idg[idg$liv&idg$dday<=0, c('age', 'qt')]),
        by='j', all=T, suffixes=c('_from', '_to')
      )
      # 2.3:  Infected quantity is a Poisson random draw
      if(randommode){
        idg2$new_inf <- apply(X=idg2[, c('qt_from', 'qt_to', 'ctrans')], MARGIN=1, function(r) rbinom(1,r[[1]]*r[[2]],r[[3]]))
      } else {
        idg2$new_inf <- idg2$qt_from*idg2$qt_to*idg2$ctrans
      }
      # 2.4:  Aggregation by age of target people
      idg2 <- data.frame(liv=T, dday=0, aggregate(new_inf~age, data=idg2, FUN=sum))
      # 2.5:  Report in the main table of the random number of infested people
      idg <- merge(
        x=idg, y=idg2,
        by=c('liv', 'age', 'dday'),
        all.x=T, all.y=F
      )
      idg$new_inf[is.na(idg$new_inf)] <- 0
      # 2.6:  Final infected quantity are min between the draw and the recevable population
      idg$new_inf <- pmin(idg$new_inf, idg$new_qt)
      # 2.7:  Update of demography
      idg$new_qt <- idg$new_qt-idg$new_inf
    } else {
      idg$new_inf <- 0
    }
  } else {
    idg$new_death <- 0
    idg$new_inf <- 0
  }
  
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
# Try: SIMU_PERIOD(as.Date('2020-01-01'), INIT_IDG(demo=D$dFR, init=H$init), H$pdday, H$beta, H$sr)

# dt_end=S$dend; sce=S$sce; idg=INIT_IDG(demo=D$dFR, init=H$init); pdday=H$pdday; beta=H$beta; sr=H$sr
# idg=EVOL_LAST_DAY(idg=idg, pdday=pdday, beta=beta, sr=sr)
SIMU_PERIOD <- function(dt_end, sce, idg, pdday, beta, sr){
  if(dt_end<=idg$dt[1]){
    stop('You need to provide a further date of end.')
  } else {
    i <- sce$dt==idg$dt[1]
    if(sum(i)<1){
      stop(paste0('Date ', idg$dt[1], ' is missing in the scenario parameter'))
    } else {
      # Effect of the scenario
      beta$value <- beta$value*sce$coeff[i][1]
      LAST_DAY <- EVOL_LAST_DAY(idg=idg, pdday=pdday, beta=beta, sr=sr)
      if(LAST_DAY$dt[1]==dt_end) return(LAST_DAY) else return(rbind(LAST_DAY, SIMU_PERIOD(dt_end, sce, LAST_DAY, pdday, beta, sr)))
      }
    }
}
