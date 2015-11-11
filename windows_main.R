#CURRENTLY IN THE "SLOW AND CORRECT" STAGE
#optimize AFTER we confirm it works

#clear all variables
rm(list=ls())
#Load libraries
library(Cairo)

#Set appropriate working directory
if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin"){ #If it's collin
  if(Sys.info()[1]=="Linux"){
    setwd("/home/collin/Dropbox/Grad school/research projects/yang_cue")
  }else{
    setwd("C:\\Repos\\phenology-cues") #desktop
  }
}else{
  if(Sys.getenv("COMPUTERNAME")=="ENT-YANG01"){
    setwd("C:\\Users\\lhyang.ent-yang01\\SkyDrive\\Phenology simulation\\phenology-cues")#desktop
  }else{  
    setwd("C:\\Users\\lhyang\\Skydrive\\Phenology simulation\\phenology-cues")} #laptop
}
#Load sources file(s)
source("windows_subs.R")

#########################
# Simulation parameters #
#########################
generations=24
duration=10
best.temp=15; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=55; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
N=40 #number of individuals 
start<-data.frame(  #this represents the min and max values used when randomly assigning initial values to the population 
  constmin=0,constmax=50,
  daymin=0,daymax=5,
  tempmin=0,tempmax=20,
  precipmin=0,precipmax=5) 
sds<-data.frame( #standard deviations for trait mutations. Currently set so that variance = max initial trait value
  const=sqrt(start$constmax),
  day=sqrt(start$daymax),
  temp=sqrt(start$tempmax),
  precip=sqrt(start$precipmax))
mutrate<-data.frame( #probability of each trait mutating in an individual. Mutations are independent of one another
  const=.1,
  day=.1,
  temp=.1,
  precip=.1)



######################################################
# Import sequence of years - LOUIE'S STUFF GOES HERE #
######################################################
years.list=NULL #Replace this with code to grab a list of data frames. Each data frame is a year.
# Each year data frame has $day, $precip, $tmean, $tmax, $tmin
# This will be the same list for all configurations of years - this is essentially just our year database
years.index=NULL # This is the list of which year.list data to use for each generation of the model



######################
# Fitness generation #
######################
#For now, daily incremental fitness will be found by multiplying two gaussian functions together:
#  one for temp, that's maximized at best.temp with sd tempsd
#  the other for precip that's maximized at best.precip with sd precipsd
# We will then normalize the results to vary from 0 to 1
for(i.year in 1:length(years.list)){
  daily.fit=dnorm(years.list[[i.year]]$tmax,mean=best.temp,sd=sd.temp)*dnorm(years.list[[i.year]]$precip,mean=best.precip,sd=sd.precip)
  daily.fit=(daily.fit-min(daily.fit))/(max(daily.fit)-min(daily.fit))
  years.list[[i.year]]=cbind(years.list[[i.year]], fit.daily=daily.fit)
}

#######################
# initializing population
#######################
##intialize a population of N individuals
# Their min and max values are determined by the start$ parameters
b.const<-runif(n=N,min=start$constmin,max=start$constmax)
b.day<-runif(n=N,min=start$daymin,max=start$daymax)
b.temp<-runif(n=N,min=start$tempmin,max=start$tempmax)
b.precip<-runif(n=N,min=start$precipmin,max=start$precipmax)
pop<-data.frame(b.const,b.day,b.temp,b.precip)
pop<-selection(newpop,duration,cur.year,N)

## Run Simulation
pophistory=runSim(startpop=pop,years.list=years.list,
                  years.ind=years.ind,N=N,duration=duration,
                  sds=sds,mutrate=mutrate,generations=generations)