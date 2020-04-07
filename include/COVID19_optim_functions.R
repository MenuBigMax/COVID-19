#y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')]; idg=IDG; tm=H$all; sr=H$sr; cdday=H$cdday
LOSS_F <- function(y, idg, tm, sr, cdday){
  y_ <- SIMU_PERIOD(dt_end=max(y$dt), idg=idg, tm=tm, sr=sr, cdday=cdday)
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

OPTIMIZE_PARAM <- function(theta, nmax, eps, p, target, init, y, demo, tm, sr, cdday){
  if(nmax>0){
    if(theta=='init'){
      theta_ <- MOVE_INIT(p, init=init, demo=demo, as.list=T)
      LOSS <- lapply(
        X= theta_,
        FUN=function(a) data.frame(a, loss=LOSS_F(y=y, idg=INIT_IDG(demo, a), tm=tm, sr=sr, cdday=cdday))
      )
    } else if(theta=='tm') {
      theta_ <- MOVE_TM(p, tm=tm)
      LOSS <- lapply(
        X=theta_,
        FUN=function(a) data.frame(a$coef, loss=LOSS_F(y=y, idg=INIT_IDG(demo = demo, init=init), tm=a$tab, sr=sr, cdday=cdday))
      )
    } else{
      return(NULL)
    }
    LOSS <- do.call(rbind, LOSS)
    newloss <- min(LOSS$loss)
    # If the new optimum is better that we already know, we continue
    # If not, we return the result
    if((newloss-eps)<target){
      if(theta=='init') init <- theta_[[which.min(LOSS$loss)]]
      if(theta=='tm') tm <- theta_[[which.min(LOSS$loss)]]$tab
      NEXT_LOSS=OPTIMIZE_PARAM(theta=theta, nmax=nmax-1, eps=eps, p=p, target=newloss, init=init, y=y, demo=demo, tm=tm, sr=sr, cdday=cdday)
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

MOVE_INIT <- function(n, init, demo, as.list=F){
  nth <- runif(n)
  nth <- rgeom(n, 0.5)*(-1*(nth<0.5)+1*(nth>=0.5))
  out <- data.frame(
    liv=T,
    dt=min(init$dt)+nth,
    age=sample(unique(demo$age), n, replace=T),
    dday=1+rgeom(n, 0.5),
    qt=1+rgeom(n,0.5)
  )
  out <- out[!duplicated(out),]
  return(if(as.list) split(out, 1:dim(out)[1]) else out)
}

# Move transmisson number matrix by age

MOVE_TM <- function(n, tm, max=2){
  out <- data.frame(b0=pmax(min(tm$ntrans), runif(n,0,max)), b1=runif(n,0,max), b2=runif(n,0,max))
  out <- split(out, 1:dim(out)[1])
  out <- lapply(out, function(a) list(coef=a, tab=data.frame(tm[,c('agent', 'target')], ntrans=a$b0+a$b1*tm$ntrans+a$b2*tm$ntrans^2)))
  return(out)
}