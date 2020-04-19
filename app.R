rm(list=ls())

library(shiny)
library(reshape)
library(ggplot2)

#Windows: setwd("C:/Users/Maxime/OneDrive/400_IT_PROJECT/RSHINY_APP/COVID-19")
#Mac: setwd('/Users/maxime/OneDrive/400_IT_PROJECT/RSHINY_APP/COVID-19')

source('include/COVID19_functions.R')

#------------------------------------------------------------------------
# 1.0 - HYPOTHESIS
#------------------------------------------------------------------------

H <- list(
  sr=read.csv('data/20200406_DEATHRISK.csv', header=T, sep=';', dec=',', encoding='utf-8'),
  pdday=read.csv('data/20200408_DDAY_PARAM.csv', header=T, sep=';', dec=',', encoding='utf-8'),
  sce=read.csv(file='data/PARAM_TRANSMISS_SCENARIO.csv',header=T, sep=';', dec='.', encoding='utf-8'),
  home=UNPIV_TRIANGL_CSV(file='data/20200403_TRANSMISS_HOME.csv', col=c('agent', 'target', 'ntrans'), sep=';', dec=',', encoding='utf-8'),
  activity=UNPIV_TRIANGL_CSV(file='data/20200403_TRANSMISS_ACTIVITY.csv', col=c('agent', 'target', 'ntrans'), sep=';', dec=',', encoding='utf-8'),
  outdoor=UNPIV_TRIANGL_CSV(file='data/20200416_TRANSMISS_OUTDOOR.csv', col=c('agent', 'target', 'ntrans'), sep=';', dec=',', encoding='utf-8'),
  init=read.csv('data/PARAM_INIT_EPID.csv', header=T, sep=';', dec=',', encoding='utf-8')
)

H$ddaymax <- max(GET_DDAY(H$pdday)$dday)+1



# Conversaion of date field saved as characters in .csv
H$init$dt <- as.Date(H$init$dt, format='%Y-%m-%d')
H$init$liv <- H$init$liv=='VRAI'
if(length(unique(H$init$dt))>1) warning('Be careful, this application doses not manage multi-date initialisation.')
H$sce$dtstart <- as.Date(H$sce$dtstart, format='%Y-%m-%d')
H$sce$dtend <- as.Date(H$sce$dtend, format='%Y-%m-%d')
#------------------------------------------------------------------------
# 1.0 - DATA
#------------------------------------------------------------------------

# One country for the moment: FRANCE
D <- list(dFR=read.csv('data/20190101_INSEE_demo_FR.csv', header=T, sep=';', dec=','))
D$dFR <- aggregate(data=D$dFR, qt~age, FUN=sum)


try(download.file(url = 'https://www.data.gouv.fr/fr/datasets/r/0b66ca39-1623-4d9c-83ad-5434b7f9e2a4', destfile='data/chiffres-cles.csv'))
# Main KPI (death, infected, immunized) reported
D$KPI <- read.csv('data/chiffres-cles.csv', header=T, sep=',')
D$KPI <- aggregate(
  data=D$KPI[D$KPI$maille_code=='FRA',c('date', 'deces', 'deces_ehpad', 'hospitalises', 'gueris')],
  .~date,
  FUN=sum,
  na.rm=T, na.action=na.pass
  )
D$KPI$tout_deces <- with(D$KPI, deces-deces_ehpad)
D$KPI <- D$KPI[, c('date', 'tout_deces', 'hospitalises', 'gueris')]
colnames(D$KPI) <- c('dt', 'death', 'infected', 'immunized')
D$KPI$dt <- as.Date(D$KPI$dt, format='%Y-%m-%d')

# Unpivot
D$KPI <- melt(D$KPI, id.vars='dt')
colnames(D$KPI) <- c('dt', 'kpi', 'qt')
D$KPI <- rbind(
  cbind(D$KPI, cumul=T),
  cbind(UNCUMUL(data=D$KPI, dtcol = 'dt', qtcol = 'qt', bycol='kpi'), cumul=F))

#------------------------------------------------------------------------
# 2.0 - INFECTIO-DEMOGRAPHY MATRIX : INITIALISATION
#------------------------------------------------------------------------

# The spread of virus is summed up by a matrix containing for each quadruplet
# (flg_alive, age, dday, dt) the number of individual concerned (qt)

# At the initialization (dt=1), nobody is:
#   - alive (liv=T)
#   - in his age cluster (given by demography)
#   - not infected (dday=0)


#------------------------------------------------------------------------
# 3.0 - OPTIMIZATION
#------------------------------------------------------------------------

O <- list(
  init=read.csv(file='data/OPTI_init.csv', sep=',', dec='.'),
  tm=read.csv(file='data/OPTI_TM.csv', sep=',', dec='.')
)
O$graph_init <- aggregate(data=O$init, loss~dt+age, FUN=median)

#------------------------------------------------------------------------
# 4.0 - PREDICTION
#------------------------------------------------------------------------

# New transmission matrix
#H$all <- ADD_TRANS_TAB(H[c('household', 'activity')])

# Main calulation
S <- list(dend=Sys.Date()+45)
S$sce.time <- GET_SCENAR_TIME(sce=H$sce, dtstart=min(H$init$dt)-90, dtend=S$dend)
S$all.trans <- rbind(
  data.frame(name='3-home', H$home),
  data.frame(name='2-activity', H$activity),
  data.frame(name='1-outdoor', H$outdoor)
  )
S$sce <- merge(S$sce.time, y=S$all.trans, all.x=T, all.y=F, by='name')
S$sce$ntrans <- S$sce$coeff*S$sce$ntrans
S$sce <- aggregate(data=S$sce, ntrans~dt+agent+target, FUN=sum)
S$res <- SIMU_PERIOD(S$dend, idg=INIT_IDG(demo=D$dFR, init=H$init), sce=S$sce, sr=H$sr, pdday=H$pdday)

# PLOT PREPARATION

# Rows in result
PL <- list(
  new.died=list(N='new died', S=(!S$res$liv)&(S$res$dday==1), A=(!D$KPI$cumul)&(D$KPI$kpi=='death')),
  tot.died=list(N='total died', S=(!S$res$liv), A=(D$KPI$cumul)&(D$KPI$kpi=='death')),
  tot.livi=list(N='total living', S=S$res$liv, A=NULL),
  new.ifctd=list(N='new infected', S=S$res$liv&S$res$dday==1, A=(!D$KPI$cumul)&(D$KPI$kpi=='infected')),
  tot.ifctd=list(N='total infected', S=S$res$liv&S$res$dday>0&S$res$dday<=H$ddaymax, A=(D$KPI$cumul)&(D$KPI$kpi=='infected')),
  tot.ifctb=list(N='total infectable', S=S$res$liv&S$res$dday==0, A=NULL),
  new.imnzd=list(N='new immunized', S=S$res$liv&S$res$dday==(H$ddaymax+1), A=(!D$KPI$cumul)&(D$KPI$kpi=='immunized')),
  tot.imnzd=list(N='total immunized', S=S$res$liv&S$res$dday>H$ddaymax, A=(D$KPI$cumul)&(D$KPI$kpi=='immunized')),
  tot.imnzb=list(N='total immunizable', S=S$res$liv&S$res$dday<=H$ddaymax, A=NULL)
)
PL <- lapply(PL, FUN=function(a) {
  out <- ggplot()
  if(sum(a$S)>0){
    tb <- aggregate(qt~age+dt, data=S$res[a$S,], FUN = sum)
    out <- out+geom_col(data=tb, aes(x=dt, y=qt, fill=age, color=age), na.rm=TRUE)
  }
  if(sum(a$A)>0){
    tb <- aggregate(qt~dt, data=D$KPI[a$A,], FUN=sum)
    out <- out+geom_line(data=tb, aes(x=dt, y=qt), size = 2)
  }
  # Title
  out <- out+ggtitle(paste0('Evolution of ', a$N, ' over time'))
  # Start and end of analysys
  out <- out+xlim(H$init$dt[1], S$dend)
  # Current day
  out <- out+geom_vline(aes(xintercept=Sys.Date()), linetype='dashed', size=2)
  return(out)
  })


#------------------------------------------------------------------------
# 4.0 - SHINY APP: DESIGN
#------------------------------------------------------------------------

# Graphical settings
G <- list(
  tithyp=lapply(levels(S$sce.time$name), function(a) ggtitle(paste0('Mean daily new transmission by age in ', a))),
  color=list(
    tm=scale_fill_gradient(low='blue', high='red', limit=c(0,1))
  )
)


ui <- navbarPage(
  title='COVID-19 PANDEMIA EVOLUTION',
  tabPanel('About',
           h1( 'Purpose'),
           'Predict the evolution of COVID-19 pandemia as accurate as possible.',
           'Better understand what are the most efficient strategy to tackle Covid-19.',
           h1('Step line'),
           h2(style='color:green;','1 - Simulation'),
           'Simulate the evolution of COVID-19 pandemia by using several hypothesis (sheet 1), unanimous data (sheet 2).',
           h2(style='color:orange;','2 - Data collection'),
           'Collect the most real (unanimous) data about COVID-19 as possible. Sources: ',
           a(href='https://www.data.gouv.fr/fr/datasets/chiffres-cles-concernant-lepidemie-de-covid19-en-france/', 'data.gouv'),
           ' ,',
           a(href='https://www.insee.fr/fr/statistiques/1892088?sommaire=1912926', 'INSEE'),
           h2(style='color:orange;', '3 - Estimation'),
           'Fit the hypothesis paramaters to real data by machine learning/optimization.',
           h2(style='color:green;','4 - Planning'),
           'With this optimized parameters, predictions of different metrics about pandemia',
           h2(style='color:red;', '5 - Analysis'),
           'Analyse effects of hypothesis parameters variation',
           textOutput('TEST')
           
  ),
  tabPanel('Estimation',
    tabsetPanel(
      type='tabs',
      tabPanel('Death', plotOutput('s1'), plotOutput('s2'), plotOutput('s3')),
      tabPanel('Infection', plotOutput('s4'), plotOutput('s5'), plotOutput('s6')),
      tabPanel('Immunization',  plotOutput('s7'), plotOutput('s8'), plotOutput('s9'))
    )
  ),
  tabPanel('Hypothesis',
      tabsetPanel(
        type='tabs',
        tabPanel('Initialization',
                 h3('First infection ?'),
                 plotOutput('o5')
                 ),
        tabPanel('Social',
                 h3('Scenario of diffusion in the area'), plotOutput('h11'),
                 h3('How many people do each agent infect...'),
                 plotOutput('h12'),plotOutput('h13'),plotOutput('h14')
        ),
        tabPanel('Virus',
                 h3('Total risk of death when infected'), plotOutput('h21'),
                 h3('Transmission coefficient and risk of death when infected'), plotOutput('h22'),
                 'pdeath: probability of death during the disease',
                 'ctrans: transmission coefficient. Factor to evaluate the contagiousnes'
        )
      )
  )
)


#------------------------------------------------------------------------
# 5.0 - SHINY APP: SERVER
#------------------------------------------------------------------------
server <- function(input, output){
  # simulation output
  output$s1 <- renderPlot(PL[[1]])
  output$s2 <- renderPlot(PL[[2]])
  output$s3 <- renderPlot(PL[[3]])
  output$s4 <- renderPlot(PL[[4]])
  output$s5 <- renderPlot(PL[[5]])
  output$s6 <- renderPlot(PL[[6]])
  output$s7 <- renderPlot(PL[[7]])
  output$s8 <- renderPlot(PL[[8]])
  output$s9 <- renderPlot(PL[[9]])
  
  # Hypothesis output
  output$h11 <- renderPlot({
    ggplot(data=S$sce.time, aes(x=dt, y=coeff, fill=name, color=name))+geom_col(na.rm=TRUE)+G$vert$init+ggtitle('Scenario of the transmission type over time')
  })
  output$h12 <- renderPlot({
    ggplot(data=H$home, aes(x=agent, y=target))+geom_tile(aes(fill=ntrans), colour='white')+G$color$tm+G$tithyp[[3]]
  })
  output$h13 <- renderPlot({
    ggplot(data=H$activity, aes(x=agent, y=target))+geom_tile(aes(fill=ntrans), colour='white')+G$color$tm+G$tithyp[[2]]
  })
  output$h14 <- renderPlot({
    ggplot(data=H$activity, aes(x=agent, y=target))+geom_tile(aes(fill=ntrans), colour='white')+G$color$tm+G$tithyp[[1]]
  })
  output$h21 <- renderPlot({
    ggplot(data=H$sr, aes(x=age, y=sr))+geom_col(na.rm=TRUE)+ggtitle('Death risk when infected by age')
  })
  output$h22 <- renderPlot({
    ggplot(data=melt(GET_DDAY(H$pdday)[,c('dday', 'pdeath', 'ctrans')], id.vars='dday'), aes(x=dday, y=value, fill=variable))+geom_bar(position = "dodge", stat="identity")+ggtitle('Transmission coefficient and probabilty to death by disease day')
  })
  output$h6 <- renderTable(H$init)
  # Optimization output
  output$o5 <- renderPlot({
    ggplot(data=O$graph_init, aes(x=age, y=dt))+geom_tile(aes(fill=loss), colour='white')+scale_fill_gradient(low='green', high='red')+ggtitle('Loss function by date and age')
    })
}

shinyApp(ui, server)
#library(rsconnect);deployApp()