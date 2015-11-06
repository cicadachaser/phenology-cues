#clear all variables
rm(list=ls())

#libraries
library(timeDate)
library(Cairo) #I'm not sure if we need this with the plotting removed

#System identification
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

######################
# Setting parameters #
######################

y1.lo<-12 #lowest `healthy temperature' value
y1.hi<-21 #Highest `healthy temperature' value
y2.lo<-30 # minimum `healthy rainfall'
y2.hi<-85 #max `healthy reainfall'
generations=24
duration=1

#input data
davis.daily<-read.csv("davis-data/626713.csv", header=T, na.strings="-9999")
davis.daily$PRCP<-davis.daily$PCRP/10 #precips are reported in tenths of mm
davis.daily$TMAX<-davis.daily$TMAX/10 #temps are reported in tenths of degree C
davis.daily$TMIN<-davis.daily$TMIN/10 #temps are reported in tenths of degree C
davis.daily$DATE2<-as.Date(as.character(davis.daily$DATE),format="%Y %m %d") #DATE2 is date formatted
davis.daily$JULIAN<-julian(davis.daily$DATE2,origin=as.Date("1892-12-31")) #1893-01-01 is day 1...
davis.daily$YEAR<-as.numeric(substr(davis.daily$DATE,1,4)) #simple field for year
davis.daily$MONTH<-as.numeric(substr(davis.daily$DATE,5,6)) #simple field for month
davis.daily$DAY<-as.numeric(substr(davis.daily$DATE,7,8)) #simple field for day
davis.daily<-davis.daily[,c("DATE2","JULIAN", "YEAR","MONTH","DAY","PRCP","TMAX","TMIN")] #simplified dataframe

davis.yearlist<-split(davis.daily,davis.daily$YEAR) #list of each year separated

davis.yearnames<-unique(davis.daily$YEAR) #gives a list of all the years in the data

#calculates the "day of year", i.e. Jan 1 is 1, and 12/31 is 365 
#adds a DAY.OF.YEAR column to each dataframe in the year list
for (i in 1:length(davis.yearnames)){
  davis.yearlist[[i]]$DAY.OF.YEAR<-julian(davis.yearlist[[i]]$DATE2, origin=as.Date(paste(davis.yearnames[i],"01","01",sep="-")))+1 #add +1 so that the first day of the year is 1, not zero. 
}

davis.daily<-unsplit(davis.yearlist,davis.daily$YEAR)
davis.daily.means<-aggregate(cbind(TMAX,TMIN,PRCP)~DAY.OF.YEAR, data=davis.daily, mean)

#### stopped here, still need to intialize TMAX.SS and make other sum of sq vectors....

for (i in 1:length(davis.yearnames)){
  comparison<-merge(davis.yearlist[[i]],davis.daily.means,by="DAY.OF.YEAR")
  TMAX.SS[i]<-sum(comparison$TMAX.x-comparison$TMAX.y)^2
}



#dat<-read.csv("davis.csv", header=T) #this is just a placeholder of monthly climate data for now

##first niche dimension (e.g. temperature)
y1<-dat$tmean
month<-dat$month
model1<-loess(y1~month, span=.35);
xv1<-seq(0,12,0.0001)
yv1<-predict(model1,data.frame(month=xv1))

#fitness function
y1.opt<-mean(c(y1.lo,y1.hi)) #optimal temp value is at mid-point
W1r<-dnorm(yv1,mean=y1.opt,sd=y1.opt) #higher W values for conditions near the opt 
W1<-(W1r-min(W1r))/(max(W1r)-min(W1r)) #rescaled between 0 to 1

##second niche dimension (e.g. precipitation)
y2<-dat$precip
model2<-loess(y2~month, span=.35);
xv2<-seq(0,12,0.0001)
yv2<-predict(model2,data.frame(month=xv2))

#fitness function
y2.opt<-mean(c(y2.lo,y2.hi)) #optimal value is at mid-point
W2r<-dnorm(yv2,mean=y2.opt,sd=y2.opt) #higher W values for conditions near the opt
W2<-(W2r-min(W2r))/(max(W2r)-min(W2r)) #rescaled between 0 to 1

##combining two fitness curves
Wr<-W1*W2 #raw fitness
W<-2*(Wr-mean(range(Wr)))/(max(Wr)-min(Wr)) #defining the combined fitness landscale, rescaled between -1 and 1 to prevent long-lived strategies

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
