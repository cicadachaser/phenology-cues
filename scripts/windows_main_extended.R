#Script for carrying out simulations of evolution of emergence decisions
#Uses windows_subs.R, windows_save.R, and windows_plot.R

#clear all variables
rm(list=ls())
set_wrkdir<-function(){
  #function for setting working directory to the right place given the current computer/user
  if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin" || Sys.getenv("USERNAME")=="Collin.work"){ #If it's collin
    if(Sys.info()["nodename"]=="DESKTOP-D6QSU8F"){
      setwd("G:\\Repos\\phenology-cues") #desktop
    }else{
      setwd("C:\\Repos\\phenology-cues") #desktop
    }
  }else{
    if(Sys.getenv("COMPUTERNAME")=="ENT-YANG01"){
      setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")#desktop
    }else{
      setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")} #laptop
  }
}

ptm <-proc.time()
#########################
# Simulation parameters #
#########################
runType="standard" ##THIS DETERMINES WHAT KIND OF YEARS WE'RE USING!
traits=c("day","temp","precip")
numsims=5# number of simulations of each type to do
runsnames=c("-earlylate50-","-punctual50-") #string without spaces (for simplicity)
#traits=c("day","temp","precip","cutemp","cuprecip","daysq","tempsq","precipsq","cutempsq","cuprecipsq")
#unitTestConst is for running the population through a unit test with the same gaussian fitness every year
#and constant environmental conditions
#unitTestRand will be for running the populations through a
#unit test with the same gaussian fitness every year and random envi conditions
#standard is for running the populations through a set of replications of the first 10 good years of the davis data
plotExtra=TRUE # do we plot snapshots of emergence through time?
plotPheno=FALSE # do we plot snapshots of phenotype through time?
viewLength=500 #for comparisons of simulation types,
#  how many generations (starting from the final and working backwards) to plot/compare
duration=10 #number of days organizm is emerged.
N=100 #number of individuals
numYears=1000 #number of years to simulate
burnIn=200 #number of years to not plot (to avoid scale issues from broad initial population traits)
best.temp=20; sd.temp=5; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=10; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
mutdist=.01 #What fraction of the total "cue space" should mutations (on average) traverse (kinda).
# if we're looking at b.day, the cue space is 365 days. So a mutdist of .01 means that each mutation will
# be drawn from a normal distribution with mean 0 and sd 3.65
set_wrkdir()
years.indlist=read.csv("enviromental histories/earlylate50.csv")
years.indlist=years.indlist$x-1913 #kludge to turn 1900s values into 1-100
years.indmat=matrix(sample(years.indlist,size=numYears*numsims,replace=TRUE),ncol=numYears,nrow=numsims) # This is the list of which year.list data to use for each generation of the model
#######################################
# Handling libraries and source files #
#######################################

#libraries
library(timeDate)
# library(Cairo) #I'm not sure if we need this with the plotting removed
library(zoo)
#Set appropriate working directory
set_wrkdir()
#Load sources file(s)
source("scripts/windows_subs.R")
###############################
set_wrkdir()
source("scripts/rate_setup.R")  #this sets up mutation rates, distances, etc.
###############################
# Generate environmental data #
###############################
if(runType=="standard"){
  years.list=yeargen.davis(best.temp = best.temp,sd.temp = sd.temp,
                           best.precip = best.precip,sd.precip = sd.precip)
} else if(runType=="unitTestConst"){
  out=yeargen.const(numYears)
  years.list=out[["years.list"]]
  years.index=rep(1,numYears)
} else if (runType=="unitTestRand"){
  out=yeargen.rand(numYears)
  years.list=out[["years.list"]]
}

store.mean=store.max=matrix(0,nrow=numsims*length(runsnames),ncol=numYears)
#matrices for storing mean and maximum possible fitness
store.names=rep(0,numsims*length(runsnames)) #vector for storing run names, corresponds to rows of store.mean
finalpops=NULL #for storing the final populations of each run.

count=1 #for tracking year in years.indmat
for(i.sim in 1:numsims){
  years.index=years.indmat[count,]
  runName=sprintf("%s%d",runsnames[1],i.sim)

  set_wrkdir()
  source("scripts/sim_runner.R")  #this sets up mutation rates, distances, etc.
  #####################
  #Saving our results #
  #####################
  #Set appropriate working directory
  set_wrkdir()
  #We have a "save data" script called windows_save.R
  source("scripts/windows_save.R")

  ############
  # Plotting #
  ############
  set_wrkdir()
  source("scripts/windows_plot.R")
  store.mean[i.sim,]=meanfit
  store.max[i.sim,]=maxfit
  store.names[i.sim]=runName
  temppop=cbind(run=rep(runName,N),pophistory[[numYears]])
  finalpops=rbind(finalpops,temppop)
  count=count+1
}

#################################################################################################
#################################################################################################
#################################################################################################
######################################
#HERE WE MAKE CHANGES FOR RUN TYPE 2!#
######################################
years.indlist=read.csv("enviromental histories/punctual50.csv")
years.indlist=years.indlist$x-1913 #kludge to turn 1900s values into 1-100
years.indmat=matrix(sample(years.indlist,size=numYears*numsims,replace=TRUE),ncol=numYears,nrow=numsims) # This is the list of which year.list data to use for each generation of the model
plotPheno=FALSE

#######################
# Reworking prep work #
#######################
set_wrkdir()
source("scripts/rate_setup.R") #this sets up mutation rates, distances, etc.
########################################################################################

#Okay, run the second run type
count=1 #for tracking row in years.indmat
for(i.sim in (numsims+1):(2*numsims)){
  years.index=years.indmat[count,]
  runName=sprintf("%s%d",runsnames[2],i.sim-numsims)
  set_wrkdir()
  source("scripts/sim_runner.R")  #this sets up mutation rates, distances, etc.
  #####################
  #Saving our results #
  #####################
  #Set appropriate working directory
  set_wrkdir()
  #We have a "save data" script called windows_save.R
  source("scripts/windows_save.R")

  ############
  # Plotting #
  ############
  set_wrkdir()
  source("scripts/windows_plot.R")
  store.mean[i.sim,]=meanfit
  store.max[i.sim,]=maxfit
  store.names[i.sim]=runName
  temppop=cbind(run=rep(runName,N),pophistory[[numYears]])
  finalpops=rbind(finalpops,temppop)
  count=count+1
}

set_wrkdir()
source("scripts/compare_runs.R") #plots comparisons between the two runs