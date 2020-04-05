rm(list=ls())

library(shiny)
library(reshape)
library(ggplot2)

#setwd("C:/Users/Maxime/OneDrive/400_IT_PROJECT/RSHINY_APP/COVID-19")

source('include/COVID19_functions.R')

#------------------------------------------------------------------------
# 1.0 - DATA
#------------------------------------------------------------------------

#----------------------------
# 1.1 - HYPOTHESIS
#----------------------------

H <- list(
  sr=read.csv('data/20200328_SURVIVALRISK.csv', header=T, sep=';', dec=',', encoding='utf-8'),
  cdday=read.csv('data/20200403_DDAY_COEFF.csv', header=T, sep=';', dec=',', encoding='utf-8'),
  scenario=read.csv(file='data/20200403_TRANSMISS_POLITIC_COEFF.csv',header=T, sep=';', dec=',', encoding='utf-8'),
  household=unpiv_triangl_csv(file='data/20200403_TRANSMISS_HOUSEHOLD.csv', col=c('agent', 'target', 'ntrans'), sep=';', encoding='utf-8'),
  activity=unpiv_triangl_csv(file='data/20200403_TRANSMISS_ACTIVITY.csv', col=c('agent', 'target', 'ntrans'), sep=';', encoding='utf-8')
)
# Sum of 2 transmisson matrix
H[['all']] <- ADD_TRANS_TAB(H[c('household', 'activity')])
# Conversion of boolean field saved as characters in .csv
H$cdday$liv <- H$cdday$liv=='VRAI'

#----------------------------
# 1.2 - DEMOGRAPHY
#----------------------------

# One country for the moment: FRANCE
D <- list(dFR=read.csv('data/20190101_INSEE_demo_FR.csv', header=T, sep=';', dec=','))
D$dFR <- aggregate(data=D$dFR, qt~age, FUN=sum)

#------------------------------------------------------------------------
# 2.0 - INFECTIO-DEMOGRAPHY MATRIX
#------------------------------------------------------------------------

# The spread of virus is summed up by a matrix containing for each quadruplet
# (flg_alive, age, dday, dt) the number of individual concerned (qt)

# At the initialization (dt=1), nobody is:
#   - alive (liv=T)
#   - in his age cluster (given by demography)
#   - not infected (dday=0)

IDG <- cbind(liv=T, D$dFR, dday=0, dt=1)

#------------------------------------------------------------------------
# 3.0 - SIMULATION
#------------------------------------------------------------------------

# To simulate the spread of the virus, we simulate a person infected
# So the day 1, we pull somebody in one age cluster and put in as dday=1
IDG <- rbind(IDG, data.frame(liv=T, age='50-59', qt=1, dday=1, dt=1))
i <- IDG$flg_alive&IDG$age=='50-59'&IDG$dday==0&IDG$dt==1
if(sum(i)>0) IDG$qt[i] <- IDG$qt[i]-1

# Main calulation
N <- 40
S <- SIMU_PERIOD(n=N, idg=IDG, tm=H$all, sr=H$sr, cdday=H$cdday)

# PLOTS
H$ddaymax <- max(H$cdday$dday[H$cdday$pdeath==0|H$cdday$ctrans>0])
PL <- list(
  new_infected=S$liv&S$dday==1,
  total_infectable=S$liv&S$dday==0,
  total_infected=S$liv&S$dday>0&S$dday<=H$ddaymax,
  new_immunized=S$liv&S$dday==(H$ddaymax+1),
  total_immunizable=S$liv&S$dday<=H$ddaymax,
  total_immunized=S$liv&S$dday>H$ddaymax,
  new_death=(!S$liv)&(S$dday==1),
  total_death=(!S$liv),
  total_living=S$liv
)
PL <- lapply(PL, FUN=function(a) if(sum(a)>0) aggregate(qt~dt+age, data=S[a,], FUN = sum) else S[0, c('qt', 'dt', 'age')])



#------------------------------------------------------------------------
# 4.0 - SHINY APP: DESIGN
#------------------------------------------------------------------------

ui <- navbarPage(
  title='COVID-19 PANDEMIA EVOLUTION',
  tabPanel('Home',
           h1( 'Purpose'),
           'Predict the evolution of COVID-19 pandemia as accurate as possible.',
           'Better understand what are the most efficient strategy to tackle Covid-19.',
           h1('Step line'),
           h2(style='color:green;','1 - Simulation'),
           'Simulate the evolution of COVID-19 pandemia by using several hypothesis (sheet 1), unanimous data (sheet 2).',
           h2(style='color:red;','2 - Data collection'),
           'Collect the most real (unanimous) data about COVID-19 as possible',
           h2(style='color:red;', '3 - Estimation'),
           'Fit the hypothesis paramaters to real data by machine learning/optimization.',
           h2(style='color:red;','4 - Prediction'),
           'With this optimized parameters, predictions of different metrics about pandemia',
           h2(style='color:red;', '5 - Analysis'),
           'Analyse effects of hypothesis parameters variation',
           textOutput('TEST')
           
  ),
  tabPanel('Hypothesis',
      tabsetPanel(
        type='tabs',
        tabPanel('Social',
                 h2('How many people do each infected agent infect...'),
                 h3('...in household ?'), plotOutput('h1'),
                 h3('... during their outdoor activity?'), plotOutput('h2')
        ),
        tabPanel('Virus',
                 h3('Total risk of death when infected'), plotOutput('h3'),
                 h3('Transmission coefficient and risk of death when infected'), plotOutput('h4')
        ),
        tabPanel('Initialisation: infected day 1', tableOutput('h5'))
      )
  ),
  tabPanel('Data',
           plotOutput('d1'), 'Source', a(href='https://www.insee.fr/fr/statistiques/1892088?sommaire=1912926', 'INSEE')
  ),
  tabPanel(
    'Simulation',
    tabsetPanel(
      type='tabs',
      tabPanel('Infection', plotOutput('s1'), plotOutput('s2'), plotOutput('s3')),
      tabPanel('Immunization', plotOutput('s4'), plotOutput('s5'), plotOutput('s6')),
      tabPanel('Death',  plotOutput('s7'), plotOutput('s8'), plotOutput('s9')),
      tabPanel('New simulation',
               sliderInput(
                 inputId = 'test',
                 label = 'How many days for the new simulation ?',
                 min = 1, max = 200, value = 15),
               'Work in progress...'
      )
    )
  ),
  tabPanel('Estimation', 'Work in progress ...')
)


#------------------------------------------------------------------------
# 5.0 - SHINY APP: SERVER
#------------------------------------------------------------------------
server <- function(input, output){
  
  # Hypothesis output
  output$h1 <- renderPlot({
    ggplot(data=H$household, aes(x=agent, y=target))+geom_tile(aes(fill=ntrans), colour='white')+scale_fill_gradient(low='blue', high='red', limit=c(0,2))+ggtitle('Average daily new transmission by age')
  })
  output$h2 <- renderPlot({
    ggplot(data=H$activity, aes(x=agent, y=target))+geom_tile(aes(fill=ntrans), colour='white')+scale_fill_gradient(low='blue', high='red', limit=c(0,2))+ggtitle('Average daily new transmission by age')
  })
  output$h3 <- renderPlot({
    ggplot(data=H$sr, aes(x=age, y=sr))+geom_col(na.rm=TRUE)+ggtitle('Death risk when infected by age')
  })
  output$h4 <- renderPlot({
    ggplot(data=melt(H$cdday[,2:4], id.vars='dday'), aes(x=dday, y=value, fill=variable))+geom_bar(position = "dodge", stat="identity")+ggtitle('Transmission coefficient and probabilty to death by disease day')
  })
  output$h5 <- renderTable(aggregate(qt~dday+age, data=IDG, FUN=sum))
  # Data output
  output$d1 <- renderPlot({
    ggplot(data=D$dFR, aes(x=age, y=qt))+geom_col(na.rm=TRUE)+ggtitle('Number of french people by age')
  })
  # simulation output
  output$s1 <- renderPlot(ggplot(data=PL[[1]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[1], ' by age')))
  output$s2 <- renderPlot(ggplot(data=PL[[2]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[2], ' by age')))
  output$s3 <- renderPlot(ggplot(data=PL[[3]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[3], ' by age')))
  output$s4 <- renderPlot(ggplot(data=PL[[4]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[4], ' by age')))
  output$s5 <- renderPlot(ggplot(data=PL[[5]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[5], ' by age')))
  output$s6 <- renderPlot(ggplot(data=PL[[6]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[6], ' by age')))
  output$s7 <- renderPlot(ggplot(data=PL[[7]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[7], ' by age')))
  output$s8 <- renderPlot(ggplot(data=PL[[8]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[8], ' by age')))
  output$s9 <- renderPlot(ggplot(data=PL[[9]], aes(x=dt,y=qt, fill=age, color=age))+geom_col(na.rm=TRUE)+xlim(0,N)+ggtitle(paste0('Evolution of ',names(PL)[9], ' by age')))
  
  # ---------------------------------
  # Test
  # ---------------------------------
  output$TEST <- renderText("")
  
}

shinyApp(ui, server)

#library(rsconnect);deployApp()