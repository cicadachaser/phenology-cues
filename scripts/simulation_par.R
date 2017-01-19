#simulation<-function(runsname){
rm(list=setdiff(ls(), c("set_wrkdir","runsname","runsname.old","runsnames.first","runsname.second","runsnames","start.time")))
follow.run=FALSE
require(timeDate)
require(zoo)
#Set appropriate working directory
set_wrkdir()
#Load sources file(s)
source("scripts/windows_subs.R")
###############################
set_wrkdir()
source(paste("parameters/",runsname,".R",sep="")) #read in all the parameters
#reading in fitness shape file
source(paste("fitcurve/",fitshape,".R",sep=""))
years.list=NULL
if(runType=="unitTestConst"){
  out=yeargen.const(numYears)
  years.list=out[["years.list"]]
  years.index=rep(1,numYears)
} else if (runType=="unitTestRand"){
  out=yeargen.rand(numYears)
  years.list=out[["years.list"]]
} else{
  years.stuff=yeargen(dat.file=runType, fit.parms=fit.parms,
                      baseTemp=baseTemp,
                      other.name=other.name,
                      decay=decay,
                      moist.norm=moist.norm)
  years.list=years.stuff[[1]]
  years.indlist=years.stuff[[2]]

}
#setting up all the mutation rates etc
source("scripts/rate_setup.R")  #this sets up mutation rates, distances, etc.

#Setting up year generation
#check to see if there is a set of years already generated with the appropriate name
# This will let us compare different runs across the same set of years.
set_wrkdir()
years.name=paste(runType,"-yearsmat-",yearSet,"-nyrs",numYears,"-nsims",numsims,"-label_",yearLabel,".csv",sep="")
if(file.exists(paste("yearinds/",years.name,sep=""))){
  years.indmat=read.csv(paste("yearinds/",years.name,sep=""),header=TRUE)
}else{
  years.indlist=years.indlist-1913 #kludge to turn 1900s values into 1-100
  years.indmat=matrix(sample(years.indlist,size=numYears*numsims,replace=TRUE),ncol=numYears,nrow=numsims) # This is the list of which year.list data to use for each generation of the model
  write.csv(years.indmat,paste("yearinds/",years.name,sep=""),row.names=FALSE)
}


#Stuff for storing summary stats of multiple sims
#matrices for storing mean and maximum possible fitness

#Create folder for saving results, making a backup folder to save scripts into

dir.create(paste("results/",runsname,sep=""),showWarnings = FALSE)
dir.create(paste("results/",runsname,"/rundata",sep=""),showWarnings = FALSE)
file.copy(from = paste("parameters/",runsname,".R",sep=""),
          to = paste("results/",runsname,"/rundata/",runsname,".R",sep="")) #save parameter file
#dir.create(paste("results/",runsname,"/rundata/scripts/",sep=""),showWarnings = FALSE)
file.copy(from = "scripts",
          to = paste("results/",runsname,"/rundata",sep=""),recursive=TRUE,overwrite=TRUE) #save scripts used
write.csv(years.indmat,paste("results/",runsname,"/rundata/",years.name,sep=""))
file.copy(from = paste("fitcurve/",fitshape,".R"),
          to = paste("results/",runsname,"/rundata/",fitshape,".R",sep=""))

require(doSNOW); require(parallel); require(doParallel)
nClust=detectCores(all.tests=FALSE,logical=TRUE)
c1<-makeCluster(min(nClust-1,5))
registerDoParallel(c1)
res=foreach(i.sim = 1:numsims, .packages=c("timeDate","zoo","vegan","scatterplot3d")) %dopar% {
  years.index=as.numeric(years.indmat[i.sim,])
  runName=sprintf("%s%d",runsname[1],i.sim)
  set_wrkdir()
  pophistory=simrunner(N=N,traits=traits,start=start,years.list=years.list,years.index=years.index,
                       duration=duration,sds=sds,mutrate=mutrate,numyears=numyears,fattail=fattail)
  #####################
  #Saving our results #
  #####################
  #Set appropriate working directory
  set_wrkdir()
  #We have a "save data" script called windows_save.R
  list.all=list(mget(ls(all.names=TRUE)))
  windows_save(list.all)
  set_wrkdir()
  plotres=windows_plot(list.all) #this function returns act.eff for later use
  attach(plotres)
  dev.off()

  # Note that act.eff has been calculated already for windows_plot.R
  wrk.acteff=(act.eff[act.eff[,"gen"]>(numYears-50),]) #choose only the final values of act.eff
  store.coEff=apply(wrk.acteff[,sprintf("b.%s",traits)],2,mean)
  temppop=cbind(run=rep(runName,N),pophistory[[numYears]])
  return(list(meanfit=meanfit,maxfit=maxfit,runName=runName,wrk.acteff=wrk.acteff,store.coEff=store.coEff,finalpops=temppop))
}
stopCluster(c1)
# set_wrkdir()
save(list=c("res"),file=paste("results/",runsname,"/",runsname,"_summary.RData",sep=""))
#}

