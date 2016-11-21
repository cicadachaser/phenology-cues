runFile="debug_100yr.R"

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
# source("scripts/var_exper/fig1_calculator_opt_comp.R")
source("scripts/var_exper/fig1_calculator_opt.R")
tempPath=getwd()
set_wrkdir()
file.copy(from=paste("scripts/var_exper/runfiles/",runFile,sep = ""),
          to=paste(tempPath,"/pars-version-",runFile,sep="")
          )
