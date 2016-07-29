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
runsnames=c("allyears d-ct-cp","allyears d-t-p-ct-cp")
runsname=runsnames[1]
source("scripts/simulation_par.R")
source("scripts/analytic.R")
source("scripts/analytic_image.R")
runsname=runsnames[2]
source("scripts/simulation.R")
source("scripts/analytic.R")
source("scripts/analytic_image.R")
set_wrkdir()
source("scripts/compare_runs.R")


