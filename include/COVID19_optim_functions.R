#y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')]; idg=IDG; sce=D$sce; sr=H$sr; pdday=H$pdday
LOSS_F <- function(y, idg, sce, sr, pdday){
  y_ <- SIMU_PERIOD(dt_end=max(y$dt), idg=idg, sce=sce, sr=sr, pdday=pdday)
  y_$kpi <- NA
  y_$kpi[y_$liv&y_$dday==1] <- 'infected'
  y_$kpi[y_$liv&y_$dday==(H$ddaymax+1)] <- 'immunized'
  y_$kpi[(!y_$liv)&(y_$dday==1)] <- 'death'
  y_ <- aggregate(data=y_, qt~dt+kpi, FUN=sum)
  
  y <- merge(x=y, y=y_, by=c('dt', 'kpi'), all.x=T, all.y=F)
  y$qt[is.na(y$qt)] <- 0
  y$gap <- y$qt-y$qt_min
  y$cost <- (y$gap<0)*(y$gap^4)+(y$gap>=0)*(y$gap^2)
  return(sum(y$cost))
}

# theta='init';nmax=5;eps=100;nb=5;target=OPT$L0;init=H$init;y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')];demo=D$dFR;sce=S$sce;sr=H$sr;pdday=H$pdday
OPTIMIZE_PARAM <- function(theta, nmax, eps, nb, pace, target, init, y, demo, sce, sr, pdday){
  if(nmax>0){
    if(theta=='init'){
      theta_ <- MOVE_INIT(n=nb, init=init, demo=demo, pace=pace,  as.list=T)
      LOSS <- lapply(
        X= theta_,
        FUN=function(a) data.frame(a, loss=LOSS_F(y=y, idg=INIT_IDG(demo, a), sce=sce, sr=sr, pdday=pdday))
      )
    } else if(theta=='tm') {
      theta_ <- MOVE_TM(n=nb, sce=sce, pace=pace)
      LOSS <- lapply(
        X=theta_,
        FUN=function(a) data.frame(a$coef, loss=LOSS_F(y=y, idg=INIT_IDG(demo = demo, init=init), sce=a$tab, sr=sr, pdday=pdday))
      )
    } else if(theta=='pdday'){
      theta_ <- MOVE_PDDAY(n=nb, pdday=pdday, pace=pace, as.list=T)
      LOSS <- lapply(
        X=theta_,
        FUN=function(a) data.frame(a[1,], loss=LOSS_F(y=y, idg=INIT_IDG(demo = demo, init=init), sce=sce, sr=sr, pdday=a))
      )
    } else {
      return(NULL)
    }
    LOSS <- do.call(rbind, LOSS)
    newloss <- min(LOSS$loss)
    # If the new optimum is better that we already know, we continue
    # If not, we return the result
    if((newloss-eps)<target){
      if(theta=='init') init <- theta_[[which.min(LOSS$loss)]]
      if(theta=='sce') sce <- theta_[[which.min(LOSS$loss)]]$tab
      NEXT_LOSS=OPTIMIZE_PARAM(theta=theta, nmax=nmax-1, eps=eps, nb=nb, pace=pace, target=newloss, init=init, y=y, demo=demo, sce=sce, sr=sr, pdday=pdday)
      return(rbind(LOSS, NEXT_LOSS))
    } else {
      return(LOSS)
    }
  } else {
    return(NULL)
  }
}


#------------------------------------------------------------------------
# 1.0 - MOVING FUNCTION
#------------------------------------------------------------------------

MOVE_INIT <- function(n, init, demo, pace=0.5, as.list=F){
  nth <- runif(n)
  nth <- rgeom(n, 1/pace)*(-1*(nth<0.5)+1*(nth>=0.5))
  out <- data.frame(
    liv=T,
    dt=min(init$dt)+nth,
    age=sample(unique(demo$age), n, replace=T),
    dday=1,
    qt=1
  )
  out <- out[!duplicated(out),]
  return(if(as.list) split(out, 1:dim(out)[1]) else out)
}

# Move transmisson number matrix by age

MOVE_TM <- function(n, sce, pace=0.5){
  out <- data.frame(b0=runif(n,min = 1-pace/2,max=1+pace/2), b1=runif(n,min = 1-pace/2,max=1+pace/2))
  out <- split(out, 1:dim(out)[1])
  out <- lapply(out, function(a) list(coef=a, tab=data.frame(sce[,c('dt', 'agent', 'target')], ntrans=exp(a$b0+a$b1*sce$ntrans))))
  return(out)
}

# Move transmission coefficient and death probability by day during

MOVE_PDDAY <- function(n, pdday, pace=1, as.list=F){
  out <- sapply(c(pdday$shape[1], pdday$scale[1]), function(a) pmax(0,a+runif(n=n, -pace, pace)))
  out <- data.frame(pdday='weibull', shape=out[,1], scale=out[,2])
  return(if(as.list) split(out, 1:dim(out)[1]) else out)
}