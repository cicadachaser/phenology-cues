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
#Parameters for the function
# dat.file: climate file to use, found in years-data folder
# traits: vector of strings for the cues to plot
# text.pos: Ignore this for now - optiona, used for choosing custom legend positions, but have disabled legends
# fit.file: fitness file/function to use, found in fitcurve folder
# decay: decay rate parameter for calculating moisture. Default is .2
# baseTemp: baseline temperature used when calculating cumulative temp, default is 0
# numyears: number of years to plot for the individual plots
# other.name: second climate variable to use when calculating fitness, default is "moist"
# duration: duration to use when calculating fitness, default is 10
# lag: lag to use when calculating fitness, default is 1
#Example call of the function:
x11()
exp_corr(dat.file="Davisdat.Rdata", traits=c("temp", "cutemp"), fit.file = "standgauss.R",decay=.1)
