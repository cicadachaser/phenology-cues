#Script to use exp_fitfun() to look at fitness functions

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
source("scripts/exploratory/explore_funs.R")
cat("fitness files you might use:")
list.files("fitcurve/")
cat("Climate data files we can use:")
list.files("data-years/")
x11()
#Parameters for exp_fitfun:
# fit.file: name of fitness file to use (fitness file is within "fitcurve" folder
# fit.parms: list of fitness parameters to use. It looks like R is just using the ones defined in the fitness file, though
# dat.file: file of the climate data to use
# other.name: name of other parameter to use. Default is "moist", probably keep it that way.
# n.plotyears: number of years to use when plotting points.

#Example call of the function:
exp_fitfun(fit.file="skewgauss.R",dat.file="ithacaDat.Rdata",max.temp=60,fit.parms = fit.parms)
