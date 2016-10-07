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
set_wrkdir()
runsnames=c("allyears lagtest d-t-p","strange.TMAX.25 d-t-p-ct-cp")
start<-proc.time()
runsname=runsnames[1]
#simulation<-function(runsname){
rm(list=setdiff(ls(), c("set_wrkdir","runsname","runsname.old","runsnames.first","runsname.second","runsnames","start")))
follow.run=FALSE
require(timeDate)
require(zoo)
#Set appropriate working directory
set_wrkdir()
#Load sources file(s)
source("scripts/windows_subs.R")
set_wrkdir()
source(paste("parameters/",runsname,".R",sep="")) #read in all the parameters
#reading in fitness shape file
source(paste("fitcurve/",fitshape,".R",sep=""))
years.list=NULL

if(runType=="standard"){
  years.stuff=yeargen.davis(best.temp = best.temp,sd.temp = sd.temp,
                            best.precip = best.precip,sd.precip = sd.precip)
  years.list=years.stuff[[1]]
  years.indlist=years.stuff[[2]]
} else if(runType=="ithaca"){
  years.stuff=yeargen.ithaca(best.temp = best.temp,sd.temp = sd.temp,
                             best.precip = best.precip,sd.precip = sd.precip)
  years.list=years.stuff[[1]]
  years.indlist=years.stuff[[2]]
} else if(runType=="unitTestConst"){
  out=yeargen.const(numYears)
  years.list=out[["years.list"]]
  years.index=rep(1,numYears)
} else if (runType=="unitTestRand"){
  out=yeargen.rand(numYears)
  years.list=out[["years.list"]]
}
source("scripts/rate_setup.R")  #this sets up mutation rates, distances, etc.

#Setting up year generation
#check to see if there is a set of years already generated with the appropriate name
# This will let us compare different runs across the same set of years.
set_wrkdir()
years.name=paste(runType,"-yearsmat-",yearSet,"-nyrs",numYears,"-nsims",numsims,"-label_",yearLabel,".csv",sep="")
years.index=years.indlist-1913 #kludge to turn 1900s values into 1-100


yearPrecip=yearTemp=matrix(0,nrow=length(years.index),ncol=365)
count=1
for(i in as.character(years.indlist)){
  curTemp=years.list[[i]]$temp
  if(length(curTemp)==366){curTemp=curTemp[-366]} #to handle leap years, remove last day
  yearTemp[count,]=curTemp;
  curPrecip=years.list[[i]]$precip
  if(length(curPrecip)==366){curPrecip=curPrecip[-366]} #to handle leap years, remove last day
  yearPrecip[count,]=curPrecip;
  count=count+1
}
#Calculate the mean fitness accrued each day
par(mar=c(5,5,4,3))
meanTemp=apply(yearTemp,2,mean)
meanPrecip=apply(yearPrecip,2,mean)
precipProb=apply(yearPrecip==0,2,sum)
smoothTemp=loess.smooth(1:365,meanTemp,evaluation = 365,span=1/6)
smoothPrecip=scatter.smooth(1:365,meanPrecip,evaluation = 365,span=1/3)
smoothPrecip=loess.smooth(1:365,meanPrecip,evaluation = 365,span=1/6)
#Decided to simply use a constant precip, since it has an incredibly high
#  day-to-day variance and minimal yearly pattern
meanYr=cbind(smoothTemp$x,smoothTemp$y,rep(0,365))
colnames(meanYr)<-c("day","temp","precip")
locName=runType;if(runType=="standard"){locName="davis"}

save(list=c("meanYr"),file=paste("data-years/",locName,"Meanyr.Rdata"))
