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
<<<<<<< HEAD
runsnames=c("cu-stand_par")
=======
runsnames=c("Davis-stand-gauss")
runsnames=c("greatfalls-stand-fat",
            "greatfalls-stand-gauss", "Honolulu-stand-fat", "Honolulu-stand-gauss",
            "ithaca-stand-fat", "ithaca-stand-gauss", "rapidcity-stand-fat",
            "rapidcity-stand-gauss", "davis-stand-fat","davis-stand-gauss")
>>>>>>> b808c14... fixed code, removed bad results
start.time<-proc.time()
runsname=runsnames[1]
# source("scripts/simulation.R")
source("scripts/simulation_par.R")
# source("scripts/analytic.R")
# source("scripts/analytic_image.R")
# runsname=runsnames[2]i
# source("scripts/simulation.R")
# source("scripts/analytic.R")
# source("scripts/analytic_image.R")
# set_wrkdir()
# source("scripts/compare_runs.R")
proc.time()-start.time

