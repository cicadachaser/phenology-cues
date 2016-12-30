#Setting up for parallel processing

#Small function to simplify rewriting factor to number
fac2num=function(vec){return(as.numeric(as.character(vec)))}

#Getting basic information
source(paste("scripts/var_exper/runfiles/",runFile,sep = "")) #call the runFile which as all the variables defined
startTime=proc.time()
#Creating combinations of year and day variations using matrix trick
yrstdmat=matrix(seq(0,yearstdMax,length=numpts),ncol=numpts,nrow=numpts,byrow = TRUE)
meanYr
daystdmat=matrix(seq(0,daystdMax,length=numpts),ncol=numpts,nrow=numpts,byrow = FALSE)
yearstds=yrstdmat[1:(numpts^2)]
daystds=daystdmat[1:(numpts^2)]

require(zoo) #for use in fitness calculations

#Reading in the appropriate fitness function
set_wrkdir()
source(file = paste("fitcurve/",fitfile,".R",sep=""))

#Reading in mean year, assigning to appropriate variables
load(paste("data-years/",locName,"Dat.Rdata",sep = ""),envir= e <- new.env())
precip.means=e$daily.means$PRCP.means
temp.means=e$daily.means$TMAX.means
remove(list="e")

#smoothing across three repetitions of the year, using middle one (so that boundaries touch)
smooth.tempyrs=loess.smooth(rep(1:(365*3)),rep(temp.means[1:365],3),evaluation = 365*3,span=1/18)
smooth.temp=smooth.tempyrs$y[366:(366+364)]

#Smoothing again, this time across precip
smooth.precipyrs=loess.smooth(rep(1:(365*3)),rep(precip.means[1:365],3),evaluation = 365*3,span=1/18)
smooth.precip=smooth.precipyrs$y[366:(366+364)]

#create mean year dataframe, using either the smoothed precip or all zeros for precip
if(actprecip==FALSE){
  meanYr=as.data.frame(cbind(1:365,smooth.temp,rep(0,365)))
}else{
  meanYr=as.data.frame(cbind(1:365,smooth.temp,smooth.precip))
}
colnames(meanYr)=c("day","temp","precip")

#Calculate the "base temp" for cumulative temp calculations, taking baseTempQth quantile
baseTemp=sort(meanYr$temp)[round(365/baseTempQ)]

#Pull in the appropriate functions
source("scripts/var_exper/prepstuff.R")
mutdist=0 #to avoid error in the part of rate_setup.R that we're not using
source("scripts/rate_setup.R")
source("scripts/windows_subs.R")

seedvect=sample(1:10^5,length(yearstds))
require(doSNOW); require(parallel); require(doParallel)
nClust=detectCores(all.tests=FALSE,logical=TRUE)
c1<-makeCluster(nClust-2)
if(numYears>200) c1<-makeCluster(2)
registerDoParallel(c1)

#Make a dataframe to store the results of each run
totnum=length(yearstds)*length(traitslist)*slownum #total number of runs
#Iterate through each each combination of stdev'
res=foreach(i.stdev = 1:length(yearstds)) %dopar% {
  #Re-load each library needed
  require(zoo) #for use in choosing initial points to check
  ###################################################
  # Produce a set of years of appropriate variances #
  ###################################################
  yearstd=yearstds[i.stdev]
  daystd=daystds[i.stdev]
  set.seed(seedvect[i.stdev])

  #generate years of climate
  years.list=create_yrs(meanYr,numYears,yearstd,daystd)
  #####################
  #create initial starting points, check them
  resmat=NULL
  #FOR B.DAY
  start.opt=proc.time()
  # dayres=opt_day(years.list)
  dayres=opt_day(years.list)
  resmat=rbind(resmat,c(daystd,yearstd,"day",dayres))
  # tempres=opt_temp(years.list)
  tempres=opt_temp(years.list)
  resmat=rbind(resmat,c(daystd,yearstd,"temp",tempres))
  # cures=opt_cutemp(years.list)
  cures=opt_cutemp(years.list)
  resmat=rbind(resmat,c(daystd,yearstd,"cutemp",cures))
  time.opt=proc.time()-start.opt
}


# df <- data.frame(matrix(unlist(res), ncol=5, byrow=F))
overall.res = do.call(rbind.data.frame, res)
names(overall.res)=c("daystd","yearstd","trait","geofit","traitval")
overall.res$daystd=fac2num(overall.res$daystd)
overall.res$yearstd=fac2num(overall.res$yearstd)
overall.res$trait=as.character(overall.res$trait)
overall.res$geofit=fac2num(overall.res$geofit)
overall.res$traitval=fac2num(overall.res$traitval)



# Save results, make figures
set_wrkdir()
dir.create(paste("results/fig1/",runnum,sep=""))
save.image(file=paste("results/fig1/",runnum,"/fig1dat-version",runnum,".Rdata",sep=""))
source("scripts/var_exper/fig1_plot_opt.R")
runTime=proc.time()-startTime;print(runTime)
stopCluster(c1)
