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
  pdday=read.csv('data/PARAM_DDAY.csv', header=T, sep=';', dec=',', encoding='utf-8'),
  beta=read.csv('data/PARAM_BETA.csv', header=T, sep=';', dec=',', encoding='utf-8', stringsAsFactors=FALSE),
  sce=read.csv(file='data/PARAM_TRANSMISS_SCENARIO.csv',header=T, sep=';', dec='.', encoding='utf-8'),
  home=UNPIV_TRIANGL_CSV(file='data/20200403_TRANSMISS_HOME.csv', col=c('agent', 'target', 'ntrans'), sep=';', dec=',', encoding='utf-8'),
  activity=UNPIV_TRIANGL_CSV(file='data/20200403_TRANSMISS_ACTIVITY.csv', col=c('agent', 'target', 'ntrans'), sep=';', dec=',', encoding='utf-8'),
  outdoor=UNPIV_TRIANGL_CSV(file='data/20200416_TRANSMISS_OUTDOOR.csv', col=c('agent', 'target', 'ntrans'), sep=';', dec=',', encoding='utf-8'),
  init=read.csv('data/PARAM_INIT_EPID.csv', header=T, sep=';', dec=',', encoding='utf-8')
)

H$ddaymax <- max(GET_DDAY(H$pdday, H$beta)$dday)+1



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
# 3.0 - OPTIMIZATION
#------------------------------------------------------------------------

O <- list(
  pdday=list(N='risk of death per day', D=read.csv(file='data/OPTI_PDDAY.csv', sep=',', dec='.')),
  tm=list(N='transmission matrix ', D=read.csv(file='data/OPTI_TM.csv', sep=',', dec='.')),
  init=list(N='patient 0', D=read.csv(file='data/OPTI_init.csv', sep=',', dec='.'))
)
O$pdday$G <- ggplot(data=O$pdday$D, aes(x=shape, y=scale, z=loss))
O$tm$G <- ggplot(data=O$tm$D, aes(x=b0, y=b1, z=loss))
O$init$G <- ggplot(data=aggregate(data=O$init$D, loss~dt+age, FUN=median), aes(x=age, y=dt, z=loss))
O <- lapply(O, FUN=function(a) a$G+stat_summary_2d(fun=median)+ggtitle(paste0('Accuracy depending of parameters of ', a$N))+scale_fill_gradient(low='green', high='red'))
#------------------------------------------------------------------------
# 4.0 - PREDICTION
#------------------------------------------------------------------------

# New transmission matrix
#H$all <- ADD_TRANS_TAB(H[c('household', 'activity')])

# Main calulation
S <- list(dstart=min(H$init$dt)-90, dend=Sys.Date()+45)
S$sce <- aggregate(
  data=GET_SCENAR_TIME(sce=H$sce, dtstart=S$dstart, dtend=S$dend),
  coeff~dt, FUN=sum)
S$res <- SIMU_PERIOD(
  dt_end=S$dend,
  sce=S$sce,
  idg=INIT_IDG(demo=D$dFR, init=H$init),
  pdday=H$pdday,
  beta=H$beta,
  sr=H$sr)

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

G <- lapply(
  list(home='home', activity='activity', outdoor='outdoor'),
  function(a){
    out <- ggplot(data=H[[a]], aes(x=agent, y=target, z=ntrans))
    out <- out+stat_summary_2d(fun=median)+scale_fill_gradient(low='blue', high='red', limit=c(0,1))
    out <- out+ggtitle(paste0('In the context of ', a, ', how many people does each sick person infect every day?'))
    return(out)
  })                    

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
           ', ',
           a(href='https://www.insee.fr/fr/statistiques/1892088?sommaire=1912926', 'INSEE'),
           h2(style='color:orange;', '3 - Estimation'),
           'Fit the hypothesis paramaters to real data by machine learning/optimization.',
           h2(style='color:green;','4 - Planning'),
           'With this optimized parameters, predictions of different metrics about pandemia',
           h2(style='color:red;', '5 - Analysis'),
           'Analyse effects of hypothesis parameters variation',
           textOutput('TEST')
           
  ),
  tabPanel(
    'SIR model',
    sliderInput("S0",'Healthy population at time 0', min = 0, max = 1, value = 0.99),
    sliderInput("beta", 'Infection per interaction', min=0, max=1, value=0.2),
    sliderInput("lambda", 'Days to be cured', min=0, max=120, value=21),
    textOutput('R0'),
    plotOutput('evol')
  ),
  tabPanel('Estimation',
    tabsetPanel(
      type='tabs',
      tabPanel('Death', plotOutput('s1'), plotOutput('s2'), plotOutput('s3')),
      tabPanel('Infection', plotOutput('s4'), plotOutput('s5'), plotOutput('s6')),
      tabPanel('Immunization',  plotOutput('s7'), plotOutput('s8'), plotOutput('s9'))
    )
  ),
  tabPanel('Parameters',
      tabsetPanel(
        type='tabs',
        tabPanel('Virus',
                 h2('Hypothesis for simulation'),
                 h3('Total risk of death when infected'), plotOutput('h21'),
                 h3('Transmission coefficient and risk of death when infected'), plotOutput('h22'),
                 'pdeath: probability of death during the disease (sum=1)',
                 'ctrans: transmission coefficient (between 0 and 1). Factor to evaluate the contagiousnes',
                 h2('Estimation'), plotOutput('e1')
        ),
        tabPanel('People',
                 h2('Hypothesis for simulation'),
                 h3('Scenario of diffusion in the area'), plotOutput('h11'),
                 h3('How many people do each agent infect...'),
                 plotOutput('h1'),plotOutput('h2'),plotOutput('h3'),
                 h2('Estimation'), plotOutput('e2')
        ),
        tabPanel('Initialization',
                 h2('Hypothesis for simulation'), textOutput('h_init'),
                 h2('Estimation'), plotOutput('e3')
                 )
      )
  )
)


#------------------------------------------------------------------------
# 5.0 - SHINY APP: SERVER
#------------------------------------------------------------------------
server <- function(input, output){
  # SIR model plot
  output$evol <- renderPlot({
    ggplot(
      data=NEXT_DAY_SIR(S=input$S0, I=(1-input$S0), beta=input$beta, lambda=input$lambda, dt0=as.Date('2020-02-15'), n.max=365, R.max=0.95),
      aes(x=dt, y=p, fill=type, col=type)
    )+geom_line()
  })
  output$R0 <- renderText(paste0('With this settings, R0=', input$beta*input$lambda))
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
  output$h_init <- renderText(paste0('The first infected people was ', H$init$age[1], ' years old and arrived the ', as.character(H$init$dt[1]), ', sick for ', H$init$dday[1], ' day(s)'))
  output$h11 <- renderPlot({
    ggplot(data=S$sce, aes(x=dt, y=coeff, color=coeff))+geom_col(na.rm=TRUE)+G$vert$init+ggtitle('Scenario of the transmission type over time')
  })
  output$h1 <- renderPlot({G$home})
  output$h2 <- renderPlot({G$activity})
  output$h3 <- renderPlot({G$outdoor})
  output$h21 <- renderPlot({
    ggplot(data=H$sr, aes(x=age, y=sr))+geom_col(na.rm=TRUE)+ggtitle('Death risk when infected by age')
  })
  output$h22 <- renderPlot({
    ggplot(data=melt(GET_DDAY(H$pdday, H$beta)[,c('dday', 'pdeath', 'ctrans')], id.vars='dday'), aes(x=dday, y=value, fill=variable))+geom_bar(position = "dodge", stat="identity")+ggtitle('Transmission coefficient and probabilty to death by disease day')
  })
  output$h6 <- renderTable(H$init)
  # Optimization output
 output$e1 <- renderPlot(O[[1]])
 output$e2 <- renderPlot(O[[2]])
 output$e3 <- renderPlot(O[[3]])
}

shinyApp(ui, server)
#library(rsconnect);deployApp()
