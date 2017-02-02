##Relevant place for playing with stuff is in the last few lines
#Plus the line after "#Choose year data"

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
#Choose year data
runsname="allyears lagtest d-t-p"
set_wrkdir()
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
#############################################
#############################################
#############################################
#Actual exploration stuff below here: #
#############################################
#############################################
i.year=10
cur.year=years.list[[i.year]]
out=acf(cur.year$temp,lag.max=200)
plot(out,type='l')
k=100;m=365-k;plot(cur.year$temp[1:m],cur.year$temp[(1+k):(m+k)])
