simulation<-function(runsname){
  rm(list=setdiff(ls(), c("set_wrkdir","runsname")))
  require(timeDate)
  require(zoo)
  #Set appropriate working directory
  set_wrkdir()
  #Load sources file(s)
  source("scripts/windows_subs.R")
  ###############################
  set_wrkdir()
  source(paste("parameters/",runsname,".R",sep="")) #read in all the parameters
  if(runType=="standard"){
    years.list=yeargen.davis(best.temp = best.temp,sd.temp = sd.temp,
                             best.precip = best.precip,sd.precip = sd.precip)
  } else if(runType=="unitTestConst"){
    out=yeargen.const(numYears)
    years.list=out[["years.list"]]
    years.index=rep(1,numYears)
  } else if (runType=="unitTestRand"){
    out=yeargen.rand(numYears)
    years.list=out[["years.list"]]
  }
  #setting up all the mutation rates etc
  source("scripts/rate_setup.R")  #this sets up mutation rates, distances, etc.

  #Setting up year generation
  #check to see if there is a set of years already generated with the appropriate name
  # This will let us compare different runs across the same set of years.
  years.name=paste("yearsmat-nyrs",numYears,"-nsims",numsims,"-label_",yearLabel,".csv",sep="")
  if(file.exists(paste("yearinds/",years.name,sep=""))){
    years.indmat=read.csv(paste("yearinds/",years.name,sep=""))
  }else{
    years.indlist=read.csv("enviromental histories/earlylate50.csv")
    years.indlist=years.indlist$x-1913 #kludge to turn 1900s values into 1-100
    years.indmat=matrix(sample(years.indlist,size=numYears*numsims,replace=TRUE),ncol=numYears,nrow=numsims) # This is the list of which year.list data to use for each generation of the model
    write.csv(years.indmat,paste("yearinds/",years.name,sep=""))
  }


  #Stuff for storing summary stats of multiple sims
  store.mean=store.max=matrix(0,nrow=numsims,ncol=numYears)
  #matrices for storing mean and maximum possible fitness
  store.names=rep(0,numsims) #vector for storing run names, corresponds to rows of store.mean
  finalpops=NULL #for storing the final populations of each run.

  #Create folder for saving results, making a backup folder to save scripts into

  dir.create(paste("results/",runsname,sep=""),showWarnings = FALSE)
  dir.create(paste("results/",runsname,"/rundata",sep=""))
  file.copy(from = paste("parameters/",runsname,".R",sep=""),
            to = paste("results/",runsname,"/rundata/",runsname,".R",sep="")) #save parameter file
  file.copy("scripts",paste("results/",runsname,"/rundata/scripts",sep=""),recursive=TRUE) #save scripts used
  write.csv(years.indmat,paste("results/",runsname,"/rundata/",years.name,sep=""))
  file.copy(from = paste("fitcurve/",fitshape,".R"),
            to = paste("results/",runsname,"/rundata/",fitshape,".R",sep=""))


  count=1 #for tracking year in years.indmat
  for(i.sim in 1:numsims){
    years.index=years.indmat[count,]
    runName=sprintf("%s%d",runsname[1],i.sim)

    set_wrkdir()
    source("scripts/sim_runner.R")  #this sets up mutation rates, distances, etc.
    #####################
    #Saving our results #
    #####################
    #Set appropriate working directory
    set_wrkdir()
    #We have a "save data" script called windows_save.R
    source("scripts/windows_save.R")

    ############
    # Plotting #
    ############
    set_wrkdir()
    source("scripts/windows_plot.R")
    store.mean[i.sim,]=meanfit
    store.max[i.sim,]=maxfit
    store.names[i.sim]=runName
    temppop=cbind(run=rep(runName,N),pophistory[[numYears]])
    finalpops=rbind(finalpops,temppop)
    count=count+1
  }
  set_wrkdir()
  save(list=c("store.mean","store.names","finalpops"),file=paste("results/",runsname,"/",runsname,"_summary.RData",sep=""))
}