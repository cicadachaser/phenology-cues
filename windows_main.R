#CURRENTLY IN THE "SLOW AND CORRECT" STAGE
#optimize AFTER we confirm it works

#clear all variables
rm(list=ls())
library(Cairo)
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
source("windows_subs.R")

##################
# Run parameters #
##################
generations=24
duration=10
best.temp=15; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=55; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
N=40 #number of individuals 
start<-as.data.frame(
  constmin=0,constmax=50,
  daymin=0,daymax=5,
  tempmin=0,tempmax=20,
  precipmin=0,precipmax=5
) #this represents the min and max values used when randomly assigning initial values to the population 
sds<-( #standard deviations for trait mutations. Currently set so that variance = max initial trait value
  const=sqrt(start$constmax),
  day=sqrt(start$daymax),
  temp=sqrt(start$tempmax),
  precip=sqrt(start$precipmax)
)
mutrate<-( #probability of each trait mutating in an individual. Mutations are independent of one another
  const=.1,
  day=.1,
  temp=.1,
  precip=.1
)



############################
# Import sequence of years #
############################
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
# $b.const, $b.day, $b.temp, $b.precip  
b.const<-sample(seq(start$constmin,start$constmax,0.1),N,replace=T)
b.day<-sample(seq(start$daymin,start$daymax,0.1),N,replace=T)
b.temp<-sample(seq(start$tempmin,start$tempmax,0.1),N,replace=T)
b.precip<-sample(seq(start$precipmin,start$precipmax,0.1),N,replace=T)
pop<-as.data.frame(b.const,b.day,b.temp,b.precip)
pop<-selection(pop,duration,W)
pophistory<-list(pop) #initialize the population history

## Run Simulation
pophistory=runSim(pop=pop,y1=y1,y2=y2,month=month,y1.opt=y2.opt,y2.opt=y2.opt,duration=duration,generations=generations, N=N, years.list=years.list,years.index=years.index)

#plot t.start vs. t. duration
dev.new(height=8, width=8)
par(mar=c(1, 1, 1, 1) + 0.1, mfrow=c(6,4))
for(h in 1:generations){
  with(pophistory[[h]],hist(t.start,breaks=20,main=paste("hist of start time, generation",h)))
}
