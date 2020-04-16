source('include/COVID19_optim_functions.R')

#------------------------------------------------------------------------
# 0.0 - STATE OF THE ART
#------------------------------------------------------------------------
OPT <- list()

OPT$L0 <- LOSS_F(
  y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')],
  idg=INIT_IDG(D$dFR, H$init),
  sce=S$sce, sr=H$sr, pdday=H$pdday
)


#------------------------------------------------------------------------
# 1.0 - ESTIMATION INITIALISATION
#------------------------------------------------------------------------

OPT$init <- OPTIMIZE_PARAM(
  theta='init', nmax=10, eps=100, nb=10, pace=10, target=OPT$L0,
  init=H$init,
  y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')],
  demo=D$dFR, sce=S$sce, sr=H$sr, pdday=H$pdday)
OPT$init <- OPT$init[order(OPT$init$loss, na.last=NA), ]
head(OPT$init)
tail(OPT$init)
#OPT$init <- read.csv(file='data/OPTI_init.csv', sep=',', dec='.')
OPT$graph_init <- aggregate(data=OPT$init, loss~dt+age, FUN=median)
ggplot(data=OPT$graph_init, aes(x=age, y=dt))+geom_tile(aes(fill=loss), colour='white')+scale_fill_gradient(low='green', high='red')+ggtitle('Loss function by date and age')
#write.csv(x=data.frame(dt_cal=Sys.Date(), OPT$init), file='data/OPTI_init.csv', quote=F, append=T, row.names=F)


#------------------------------------------------------------------------
# 2.0 - ESTIMATION TRANSMISSION MATRIX (WORK iN PROGRESS)
#------------------------------------------------------------------------
OPT$tm <- OPTIMIZE_PARAM(
  theta='tm', nmax=10, eps=100, nb=10, pace=10, target=OPT$L0,
  init=H$init,
  y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')],
  demo=D$dFR, sce=S$sce, sr=H$sr, pdday=H$pdday
  )
 
OPT$tm<- OPT$tm[order(OPT$tm$loss), ]
head(OPT$tm)
tail(OPT$tm)
ggplot(data=OPT$tm, aes(x=b0, y=b1, color=loss)+geom_point(shape=19)+scale_fill_gradient(low='green', high='red')+ggtitle('Loss function by coefficient')
#write.csv(x=data.frame(dt_cal=Sys.Date(), OPT$tm), file='data/OPTI_TM.csv', quote=F, append=T, row.names=F)

#------------------------------------------------------------------------
# 3.0 - ESTIMATION of DISEASE DAY PARAMETERS
#------------------------------------------------------------------------
OPT$pdday <- OPTIMIZE_PARAM(
  theta='pdday', nmax=10, eps=100, nb=10, pace=10, target=OPT$L0,
  init=H$init,
  y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')],
  demo=D$dFR, sce=S$sce, sr=H$sr, pdday=H$pdday)

OPT$pdday <- OPT$pdday[order(OPT$pdday$loss), ]
head(OPT$pdday)
tail(OPT$pdday)
test <- GET_DDAY(OPT$pdday[1, ])
plot(test$dday, test$pdeath, type='l', col='blue', ylim=c(0,1))
lines(test$dday, test$ctrans, col='green')
#write.csv(x=data.frame(dt_cal=Sys.Date(), OPT$pdday), file='data/OPTI_PDDAY.csv', quote=F, append=T, row.names=F)

