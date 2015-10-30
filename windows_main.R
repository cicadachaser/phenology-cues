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






############################
# Import sequence of years #
############################
years.list=NULL #Replace this with code to grab a list of data frames. Each data frame is a year.
# Each year data frame has $day, $precip, $tmean, $tmax, $tmin



######################
# Fitness generation #
######################
#For now, daily incremental fitness will be found by multiplying two gaussian functions together:
#  one for temp, that's maximized at best.temp with sd tempsd
#  the other for precip that's maximized at best.precip with sd precipsd
# We will then normalize the results to vary from 0 to 1
for(i.year in 1:length(years.list)){
  temp.fit=dnorm(years.list[[i.year]]$tmean,mean=best.temp,sd=sd.temp)*dnorm(years.list[[i.year]]$precip,mean=best.precip,sd=sd.precip)
  temp.fit=(temp.fit-min(temp.fit))/(max(temp.fit)-min(temp.fit))
  years.list[[i.year]]=cbind(years.list[[i.year]], fit.daily=temp.fit)
}

#######################
# initializing population
#######################


##intialize a population of N individuals
N<-40
t.start<-sample(seq(0,12,0.1),N,replace=T)
pop<-as.data.frame(t.start)
pop<-selection(pop,duration,W)
pophistory<-list(pop) #initialize the population history

## Run Simulation
pophistory=runSim(pop=pop,y1=y1,y2=y2,month=month,y1.opt=y2.opt,y2.opt=y2.opt,duration=duration,generations=generations)

#plot t.start vs. t. duration
dev.new(height=8, width=8)
par(mar=c(1, 1, 1, 1) + 0.1, mfrow=c(6,4))
for(h in 1:generations){
  with(pophistory[[h]],hist(t.start,breaks=20,main=paste("hist of start time, generation",h)))
}
