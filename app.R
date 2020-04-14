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
  sce=read.csv(file='data/20200414_TRANSMISS_SCENARIO.csv',header=T, sep=';', dec='.', encoding='utf-8'),
  household=unpiv_triangl_csv(file='data/20200403_TRANSMISS_HOUSEHOLD.csv', col=c('agent', 'target', 'ntrans'), sep=';', encoding='utf-8'),
  activity=unpiv_triangl_csv(file='data/20200403_TRANSMISS_ACTIVITY.csv', col=c('agent', 'target', 'ntrans'), sep=';', encoding='utf-8'),
  init=read.csv('data/20200406_INIT_EPID.csv', header=T, sep=';', dec=',', encoding='utf-8')
)

H$ddaymax <- max(GET_DDAY(H$pdday)$dday)+1



# Conversaion of date field saved as characters in .csv
H$init$dt <- as.Date(H$init$dt, format='%d/%m/%y')
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
colnames(D$KPI) <- c('dt', 'kpi', 'qt_min')
D$KPI <- rbind(
  cbind(D$KPI, cumul=T),
  cbind(UNCUMUL(data=D$KPI, dtcol = 'dt', qtcol = 'qt_min', bycol='kpi'), cumul=F))

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
H$all <- ADD_TRANS_TAB(H[c('household', 'activity')])

# Main calulation
S <- list(dend=Sys.Date()+21)
S$sce.time <- GET_SCENAR_TIME(sce=H$sce, dtstart=min(H$init$dt), dtend=S$dend)
S$all.trans <- rbind(data.frame(name='3-house', H$household), data.frame(name='2-aactivity', H$activity))
S$sce <- merge(S$sce, y=S$all.trans, all.x=T, all.y=F, by='name')
S$res <- SIMU_PERIOD(S$dend, idg=INIT_IDG(demo=D$dFR, init=H$init), tm=H$all, sr=H$sr, pdday=H$pdday)

# PLOTS

PL <- list(
  new_infected=S$res$liv&S$res$dday==1,
  total_infectable=S$res$liv&S$res$dday==0,
  total_infected=S$res$liv&S$res$dday>0&S$res$dday<=H$ddaymax,
  new_immunized=S$res$liv&S$res$dday==(H$ddaymax+1),
  total_immunizable=S$res$liv&S$res$dday<=H$ddaymax,
  total_immunized=S$res$liv&S$res$dday>H$ddaymax,
  new_death=(!S$res$liv)&(S$res$dday==1),
  total_death=(!S$res$liv),
  total_living=S$res$liv
)
PL <- lapply(PL, FUN=function(a) if(sum(a)>0) aggregate(qt~dt+age, data=S$res[a,], FUN = sum) else S$res[0, c('qt', 'dt', 'age')])



#------------------------------------------------------------------------
# 4.0 - SHINY APP: DESIGN
#------------------------------------------------------------------------

# Graphical settings
G <- list(
  xlim=xlim(H$init$dt[1], S$dend),
  titsim=lapply(names(PL), function(a) ggtitle(paste0('Evolution of ',a, ' by age'))),
  vert=geom_vline(xintercept = Sys.Date(), linetype='dotted', size=2) 
  )


ui <- navbarPage(
  title='COVID-19 PANDEMIA EVOLUTION',
  tabPanel('Home',
           h1( 'Purpose'),
           'Predict the evolution of COVID-19 pandemia as accurate as possible.',
           'Better understand what are the most efficient strategy to tackle Covid-19.',
           h1('Step line'),
           h2(style='color:green;','1 - Simulation'),
           'Simulate the evolution of COVID-19 pandemia by using several hypothesis (sheet 1), unanimous data (sheet 2).',
           h2(style='color:orange;','2 - Data collection'),
           'Collect the most real (unanimous) data about COVID-19 as possible',
           h2(style='color:orange;', '3 - Estimation'),
           'Fit the hypothesis paramaters to real data by machine learning/optimization.',
           h2(style='color:green;','4 - Planning'),
           'With this optimized parameters, predictions of different metrics about pandemia',
           h2(style='color:red;', '5 - Analysis'),
           'Analyse effects of hypothesis parameters variation',
           textOutput('TEST')
           
  ),
  tabPanel('Hypothesis',
      tabsetPanel(
        type='tabs',
        tabPanel('Social',
                 h3('Level of different type of diffusion '), plotOutput('h1'),
                 h3('How many people do each agent infect...'),
                 plotOutput('h2'),plotOutput('h3')
        ),
        tabPanel('Virus',
                 h3('Total risk of death when infected'), plotOutput('h4'),
                 h3('Transmission coefficient and risk of death when infected'), plotOutput('h5'),
                 'pdeath: probability of death during the disease',
                 'ctrans: transmission coefficient. Factor to evaluate the contagiousnes'
        ),
        tabPanel('Day 1',
                 h3('Who was the first people infected in the country ?'), tableOutput('h6')
      )
      )
  ),
  tabPanel('Actual data',
           tabsetPanel(
             type='tabs',
             tabPanel('Virus',
                      h2('Different actual measures of the pandemia'),
                      h3('New by day'), plotOutput('d2'),
                      h3('Cumulated'), plotOutput('d3'),
                      'Source', 
                      a(href='https://www.data.gouv.fr/fr/datasets/chiffres-cles-concernant-lepidemie-de-covid19-en-france/', 'data.gouv')
             ),
             tabPanel('Demography', 
                      h3('Demography by age'),
                      plotOutput('d1'), 'Source',
                      a(href='https://www.insee.fr/fr/statistiques/1892088?sommaire=1912926', 'INSEE')
               )
             
           )
           
  ),
  tabPanel('Estimation',
          tabsetPanel(
             type='tabs',
             tabPanel('Social',
                      h2('How many people do each agent infect...'),
                      h3('...in household ?'),
                      h3('... during their outdoor activity?')
             ),
             tabPanel('Virus',
                      h3('Total risk of death when infected'),
                      h3('Transmission coefficient and risk of death when infected')
             ),
             tabPanel('Day 1',
                      h3('Who was the first people infected in the country ?'), plotOutput('o5')
             )
           )
          ),
  tabPanel(
    'Planning',
    tabsetPanel(
      type='tabs',
      tabPanel('Infection', plotOutput('s1'), plotOutput('s2'), plotOutput('s3')),
      tabPanel('Immunization', plotOutput('s4'), plotOutput('s5'), plotOutput('s6')),
      tabPanel('Death',  plotOutput('s7'), plotOutput('s8'), plotOutput('s9'))
      )
    )
  
)


#------------------------------------------------------------------------
# 5.0 - SHINY APP: SERVER
#------------------------------------------------------------------------
server <- function(input, output){
  
  # Hypothesis output
  output$h1 <- renderPlot({
    ggplot(data=S$sce.time, aes(x=dt, y=coeff, fill=name, color=name, order=))+geom_col(na.rm=TRUE)+ggtitle('Scenario of the transmission type over time')
  })
  output$h2 <- renderPlot({
    ggplot(data=H$household, aes(x=agent, y=target))+geom_tile(aes(fill=ntrans), colour='white')+scale_fill_gradient(low='blue', high='red', limit=c(0,2))+ggtitle('Average daily new transmission by age into house')
  })
  output$h3 <- renderPlot({
    ggplot(data=H$activity, aes(x=agent, y=target))+geom_tile(aes(fill=ntrans), colour='white')+scale_fill_gradient(low='blue', high='red', limit=c(0,2))+ggtitle('Average daily new transmission by age during activity')
  })
  output$h4 <- renderPlot({
    ggplot(data=H$sr, aes(x=age, y=sr))+geom_col(na.rm=TRUE)+ggtitle('Death risk when infected by age')
  })
  output$h5 <- renderPlot({
    ggplot(data=melt(GET_DDAY(H$pdday)[,c('dday', 'pdeath', 'ctrans')], id.vars='dday'), aes(x=dday, y=value, fill=variable))+geom_bar(position = "dodge", stat="identity")+ggtitle('Transmission coefficient and probabilty to death by disease day')
  })
  output$h6 <- renderTable(H$init)
  # Optimization output
  output$o5 <- renderPlot({
    ggplot(data=O$graph_init, aes(x=age, y=dt))+geom_tile(aes(fill=loss), colour='white')+scale_fill_gradient(low='green', high='red')+ggtitle('Loss function by date and age')
    })
  # Data output
  output$d1 <- renderPlot({
    ggplot(data=D$dFR, aes(x=age, y=qt))+geom_col(na.rm=TRUE)+ggtitle('Number of french people by age')
  })
  output$d2 <- renderPlot({
    ggplot(data=D$KPI[!D$KPI$cumul,], aes(x=dt, y=qt_min, color=kpi))+geom_line(size = 2)
  })
  output$d3 <- renderPlot({
    ggplot(data=D$KPI[D$KPI$cumul,], aes(x=dt, y=qt_min, color=kpi))+geom_line(size = 2)
  })
 

  # simulation output
  output$s1 <- renderPlot(ggplot(data=PL[[1]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[1]]+G$vert)
  output$s2 <- renderPlot(ggplot(data=PL[[2]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[2]]+G$vert)
  output$s3 <- renderPlot(ggplot(data=PL[[3]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[3]]+G$vert)
  output$s4 <- renderPlot(ggplot(data=PL[[4]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[4]]+G$vert)
  output$s5 <- renderPlot(ggplot(data=PL[[5]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[5]]+G$vert)
  output$s6 <- renderPlot(ggplot(data=PL[[6]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[6]]+G$vert)
  output$s7 <- renderPlot(ggplot(data=PL[[7]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[7]]+G$vert)
  output$s8 <- renderPlot(ggplot(data=PL[[8]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[8]]+G$vert)
  output$s9 <- renderPlot(ggplot(data=PL[[9]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+G$xlim+G$titsim[[9]]+G$vert)
}

shinyApp(ui, server)
#library(rsconnect);deployApp()
#library(rsconnect);deployApp()