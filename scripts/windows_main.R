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
#unitTestConst is for running the population through a unit test with the same gaussian fitness every year
#and constant environmental conditions
#unitTestRand will be for running the populations through a
#unit test with the same gaussian fitness every year and random envi conditions
#standard is for running the populations through a set of replications of the first 10 good years of the davis data
plotExtra=TRUE
runName="-distest" #string without spaces (for simplicity)
duration=10
N=100 #number of individuals
numYears=1000 #number of years to simulate
burnIn=100 #number of years to not plot (to avoid scale issues from broad initial population traits
best.temp=30; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=10; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
start<-data.frame(  #this represents the min and max values used when randomly assigning initial values to the population
  daymin=1/365*100,daymax=100,
  tempmin=0,tempmax=0,
  precipmin=0,precipmax=0)
sds<-data.frame( #standard deviations for trait mutations. Currently set so that variance = max initial trait value
  day=.1,
  temp=.1,
  precip=.1)
mutrate<-data.frame( #probability of each trait mutating in an individual. Mutations are independent of one another
  const=.01,
  day=.01,
  temp=.01,
  precip=.01)
years.index=rep(50:100,length.out=numYears) # This is the list of which year.list data to use for each generation of the model
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
# Generate environmental data #
###############################
#Based on the value of "runType", generate the appropriate type of data.
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
#######################
# initializing population
#######################
##intialize a population of N individuals
# Their min and max values are determined by the start$ parameters
# b.day<-runif(n=N,min=start$daymin,max=start$daymax)
b.day=100/runif(n=N,min=1,max=365)
b.temp<-runif(n=N,min=start$tempmin,max=start$tempmax)
b.precip<-runif(n=N,min=start$precipmin,max=start$precipmax)
newpop<-data.frame(b.day,b.temp,b.precip)
pop<-selection(newpop,duration,year=years.list[[1]],N)
###########################
## Running the Simulation #
###########################
pophistory=runSim(startpop=pop,years.list=years.list,
                  years.ind=years.index,N=N,duration=duration,
                  sds=sds,mutrate=mutrate,generations=length(years.index[-1]))
#Note: we've already used year 1 in initiating the pop

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

proc.time()-ptm
