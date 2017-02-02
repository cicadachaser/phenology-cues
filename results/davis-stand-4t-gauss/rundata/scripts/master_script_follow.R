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
runsnames.first=c("allyears d-t-p")
runsnames.second=c("strange.TMAX.25 d-t-p")
start<-proc.time()
runsname=runsnames.first[1]
source("scripts/simulation.R")
source("scripts/analytic.R")
source("scripts/analytic_image.R")
set_wrkdir()
runsname.old=runsname
runsname=runsnames.second[1]
source("scripts/simulation_follow.R")
source("scripts/analytic.R")
source("scripts/analytic_image.R")
proc.time()-start

