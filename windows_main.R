#clear all variables
rm(list=ls())
library(Cairo)
if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin"){ #If it's collin
  if(Sys.info()[1]=="Linux"){
    setwd("/home/collin/Dropbox/Grad school/research projects/yang_cue")
  }else{
    setwd("C:\\Users\\Collin\\Dropbox\\Grad school\\research projects\\yang_cue") #desktop
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
#y1<-(-cos+1)/2 #for other options: names(dataset)
y1<-dat$tmean
month<-dat$month
model1<-loess(y1~month, span=.35);
xv1<-seq(0,12,0.0001)
yv1<-predict(model1,data.frame(month=xv1))
if(graphics==TRUE){
  dev.new(1, height=12, width=8)
  par(mar=c(4, 5, 1.5, 1.5) + 0.1, mfrow=c(5,1))
  plot(month,y1,type="p",xaxp=c(0, 12, 12))
  lines(xv1,yv1)
}

#lower threshold vertical
if(graphics==TRUE){
  xco<-c(min(xv1[yv1>=y1.lo]),xv1,max(xv1[yv1>=y1.lo]))
  yco<-c(min(yv1),(as.numeric(yv1>=y1.lo)*yv1)+(as.numeric(yv1<=y1.lo)*min(yv1)),min(yv1))
  polygon(xco,yco,col="#FFA50075", border=NA)
  xlines<-xv1[abs(yv1-y1.lo)<0.001]
  ylines<-yv1[abs(yv1-y1.lo)<0.001]
  arrows(xlines,ylines,y1=min(yv1),lty=5,length=.1)
  #abline(v=xlines,lty=3)
  abline(h=y1.lo,lty=3)
}

#upper threshold vertical
if(graphics==TRUE){
  xco<-xv1
  xco<-c(min(xv1[yv1<=y1.hi]),xv1,max(xv1[yv1<=y1.hi]))
  yco<-c(min(yv1[yv1<=y1.hi]),(as.numeric(yv1<=y1.hi)*yv1)+(as.numeric(yv1>=y1.hi)*min(yv1[yv1<=y1.hi])),min(yv1[yv1<=y1.hi]))
  polygon(xco,yco,col="#87CEFA75", border=NA)
  xlines<-xv1[abs(yv1-y1.hi)<0.001]
  ylines<-yv1[abs(yv1-y1.hi)<0.001]
  arrows(xlines,ylines,y1=min(yv1),lty=5,length=.1)
  #abline(v=xlines,lty=3)
  abline(h=y1.hi,lty=3)
}
#fitness function
y1.opt<-mean(c(y1.lo,y1.hi)) #optimal temp value is at mid-point
W1r<-dnorm(yv1,mean=y1.opt,sd=y1.opt)
W1<-(W1r-min(W1r))/(max(W1r)-min(W1r)) #rescaled between 0 to 1
if(graphics==TRUE){
  plot(xv1,W1,type="l", xaxp=c(0, 12, 12))
}


##second niche dimension (e.g. precipitation)
y2<-dat$precip
model2<-loess(y2~month, span=.35);
xv2<-seq(0,12,0.0001)
yv2<-predict(model2,data.frame(month=xv2))
if(graphics==TRUE){
  plot(month,y2,type="p",xaxp=c(0, 12, 12))
  lines(xv2,yv2)
  #lower threshold vertical
  xco<-c(min(xv2[yv2>=y2.lo]),xv2,max(xv2[yv2>=y2.lo]))
  yco<-c(min(yv2),(as.numeric(yv2>=y2.lo)*yv2)+(as.numeric(yv2<=y2.lo)*min(yv2)),min(yv2))
  polygon(xco,yco,col="#FFA50075", border=NA)
  xlines<-xv2[abs(yv2-y2.lo)<0.001]
  ylines<-yv2[abs(yv2-y2.lo)<0.001]
  arrows(xlines,ylines,y1=min(yv2),lty=5,length=.1)
  abline(h=y2.lo,lty=3)
  #upper threshold vertical
  xco<-xv2
  xco<-c(min(xv2[yv2<=y2.hi]),xv2,max(xv2[yv2<=y2.hi]))
  yco<-c(min(yv2[yv2<=y2.hi]),(as.numeric(yv2<=y2.hi)*yv2)+(as.numeric(yv2>=y2.hi)*min(yv2[yv2<=y2.hi])),min(yv2[yv2<=y2.hi]))
  polygon(xco,yco,col="#87CEFA75", border=NA)
  xlines<-xv2[abs(yv2-y2.hi)<0.001]
  ylines<-yv2[abs(yv2-y2.hi)<0.001]
  arrows(xlines,ylines,y1=min(yv2),lty=5,length=.1)
  abline(h=y2.hi,lty=3)
}

#fitness function
y2.opt<-mean(c(y2.lo,y2.hi)) #optimal value is at mid-point
W2r<-dnorm(yv2,mean=y2.opt,sd=y2.opt)
W2<-(W2r-min(W2r))/(max(W2r)-min(W2r)) #rescaled between 0 to 1
if(graphics==TRUE){
  plot(xv2,W2,type="l", xaxp=c(0, 12, 12))
}

##combining two fitness curves
Wr<-W1*W2 #raw fitness
W<-2*(Wr-mean(range(Wr)))/(max(Wr)-min(Wr)) #defining the combined fitness landscale, rescaled between -1 and 1 to prevent long-lived strategies
if(graphics==TRUE){
  plot(xv2,W,type="l", xaxp=c(0, 12, 12))
  abline(h=0,lty=3)
}

##intialize a population of N individuals
N<-40
t.start<-sample(seq(0,12,0.1),N,replace=T)
t.duration<-sample(seq(0.1,6,0.1),N,replace=T) #random durations
#! I'm thinking maybe swap to treating t.duration as constant?
pop<-as.data.frame(cbind(t.start,t.duration))
pop<-pop[order(t.start),]
pop<-selection2(pop,W)
pophistory<-list(pop) #initialize the population history
if(graphics==TRUE){
  ##re-plot W  in year 1
  dev.new(2, height=4, width=8)
  par(mar=c(4, 5, 1.5, 1.5) + 0.1)
  plot(xv2,W,type="l", xaxp=c(0, 12, 12))
  abline(h=0,lty=3)
  
  #add arrows for each individual
  xlines<-pop$t.start
  ylines<-2*(pop$Ws-mean(c(min(pop$Ws),max(pop$Ws))))/(max(pop$Ws)-min(pop$Ws)) #rescale height of arrow between -1 and 1
  arrows(xlines,ylines,x1=(pop$t.start+pop$t.duration),lty=1,length=.1)
}

## Run Simulation
pophistory=runSim(pop=pop,y1=y1,y2=y2,month=month,y1.opt=y2.opt,y2.opt=y2.opt,generations=generations)

#plot t.start vs. t. duration
dev.new(height=8, width=8)
par(mar=c(1, 1, 1, 1) + 0.1, mfrow=c(6,4))
for(h in 1:generations){
  with(pophistory[[h]],plot(t.start,t.duration,xlim=c(0,12),ylim=c(0,6)))
}
