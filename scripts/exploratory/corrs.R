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
#Example call of the function:
exp_corr(dat.file="ithacaDat.Rdata")

#Full list of traits:
x11()
exp_corr(dat.file="RapidCityDat.Rdata",
         traits=c("day","year", "precip", "temp",   "moist",  "cutemp", "cuprecip",
                  "daysq",  "tempsq", "precipsq", "cutempsq", "cuprecipsq"),
         fit.file = "skewgauss.R"
)
