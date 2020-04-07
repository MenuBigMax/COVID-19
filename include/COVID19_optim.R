source('include/COVID19_optim_functions.R')

#------------------------------------------------------------------------
# 0.0 - STATE OF THE ART
#------------------------------------------------------------------------

L0 <- LOSS_F(
  y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')],
  idg=INIT_IDG(D$dFR, H$init),
  tm=H$all, sr=H$sr, cdday=H$cdday
)

OPT <- list()
#------------------------------------------------------------------------
# 1.0 - ESTIMATION INITIALISATION
#------------------------------------------------------------------------

OPT$init <- OPTIMIZE_PARAM(
  theta='init', nmax=100, eps=100, p=100, target=L0,
  init=H$init,
  y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')],
  demo=D$dFR, tm=H$all, sr=H$sr, cdday=H$cdday)
OPT$init <- OPT$init[order(OPT$init), ]
head(OPT$init)
tail(OPT$init)
#write.csv(x=data.frame(dt_cal=Sys.Date(), OPT$init), file='data/OPTI_init.csv', quote=F, append=T, row.names=F)


#------------------------------------------------------------------------
# 2.0 - ESTIMATION TRANSMISSION MATRIX (WORK iN PROGRESS)
#------------------------------------------------------------------------
OPT$tm <- OPTIMIZE_PARAM(
  theta='tm', nmax=5, eps=100, p=20, target=L0,
  init=H$init,
  y=D$KPI[!D$KPI$cumul,c('dt', 'kpi', 'qt_min')],
  demo=D$dFR, tm=H$all, sr=H$sr, cdday=H$cdday)

OPT$tm<- OPT$tm[order(OPT$tm$loss), ]
head(OPT$tm)
tail(OPT$tm)
#write.csv(x=data.frame(dt_cal=Sys.Date(), OPT$tm), file='data/OPTI_TM.csv', quote=F, append=T, row.names=F)

