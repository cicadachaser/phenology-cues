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
#List of parameters:
# dat.file: name of climate file to use, found in data-years folder
# fit.file: name of fitness file to use, found in fitcurve folder
# decay=.2: decay rate for calculating moisture
# interests: vector of strings for the name of the cues we're interested in
# baseTemp: base temperature to use for calucating cumulative temperature, defaults to zero
# numyears: number of years to plot for the individual plots
# other.name: second climate variable to use when calculating fitness, default is "moist"
# duration: duration used when calculating fitness, default is 10
# lag: lag used when calculating fitness, default is 1

#Example call of the function:
x11()
exp_covars(dat.file="davisDat.Rdata",fit.file = "standgauss.R",decay=.1,interests=c("temp","photo","fit.tot"))
