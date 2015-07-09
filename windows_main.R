#clear all variables
rm(list=ls())
library(Cairo)
if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin"){ #If it's collin
  if(Sys.info()[1]=="Linux"){
    setwd("/home/collin/Dropbox/Grad school/research projects/yang_cue")
  }else{
    setwd("C:\\Repos\\phenology-cues") #desktop
  }
}else{setwd("C:\\Users\\lhyang\\Dropbox\\Phenology simulation")} #laptop

source("windows_subs.R")
######################################
# Choosing threshold values
######################################

y1.lo<-12 #lowest `healthy temperature' value
y1.hi<-21 #Highest `healthy temperature' value
y2.lo<-30 # minimum `healthy rainfall'
y2.hi<-85 #max `healthy reainfall'
graphics=FALSE
generations=24

#input data, this is just a placeholder of monthly climate data for now
dat<-read.csv("davis.csv", header=T)

##first niche dimension (e.g. temperature)
y1<-dat$tmean
month<-dat$month
model1<-loess(y1~month, span=.35);
xv1<-seq(0,12,0.0001)
yv1<-predict(model1,data.frame(month=xv1))

#fitness function
y1.opt<-mean(c(y1.lo,y1.hi)) #optimal temp value is at mid-point
W1r<-dnorm(yv1,mean=y1.opt,sd=y1.opt)
W1<-(W1r-min(W1r))/(max(W1r)-min(W1r)) #rescaled between 0 to 1

##second niche dimension (e.g. precipitation)
y2<-dat$precip
model2<-loess(y2~month, span=.35);
xv2<-seq(0,12,0.0001)
yv2<-predict(model2,data.frame(month=xv2))

#fitness function
y2.opt<-mean(c(y2.lo,y2.hi)) #optimal value is at mid-point
W2r<-dnorm(yv2,mean=y2.opt,sd=y2.opt)
W2<-(W2r-min(W2r))/(max(W2r)-min(W2r)) #rescaled between 0 to 1

##combining two fitness curves
Wr<-W1*W2 #raw fitness
W<-2*(Wr-mean(range(Wr)))/(max(Wr)-min(Wr)) #defining the combined fitness landscale, rescaled between -1 and 1 to prevent long-lived strategies

##intialize a population of N individuals
N<-40
t.start<-sample(seq(0,12,0.1),N,replace=T)
t.duration<-sample(seq(0.1,6,0.1),N,replace=T) #random durations
#! I'm thinking maybe swap to treating t.duration as constant?
pop<-as.data.frame(cbind(t.start,t.duration))
pop<-pop[order(t.start),]
pop<-selection2(pop,W)
pophistory<-list(pop) #initialize the population history

## Run Simulation
pophistory=runSim(pop=pop,y1=y1,y2=y2,month=month,y1.opt=y2.opt,y2.opt=y2.opt,generations=generations)

#plot t.start vs. t. duration
dev.new(height=8, width=8)
par(mar=c(1, 1, 1, 1) + 0.1, mfrow=c(6,4))
for(h in 1:generations){
  with(pophistory[[h]],plot(t.start,t.duration,xlim=c(0,12),ylim=c(0,6)))
}
