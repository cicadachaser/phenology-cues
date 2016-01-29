#CURRENTLY IN THE "SLOW AND CORRECT" STAGE
#optimize AFTER we confirm it works

#NEXT STEP: WORK WITH THE "ACTUAL EFFECT SIZE" object and make graphs!

#clear all variables
rm(list=ls())
#########################
# Simulation parameters #
#########################
runType="unitTestRand" ##THIS DETERMINES WHAT KIND OF YEARS WE'RE USING!
#unitTestConst is for running the population through a unit test with the same gaussian fitness every year
#and constant environmental conditions
#unitTestRand will be for running the populations through a
#unit test with the same gaussian fitness every year and random envi conditions
#standard is for running the populations through a set of replications of the first 10 good years of the davis data
runNumber=12
duration=10
numYears=5000
best.temp=15; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=55; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
N=40 #number of individuals
start<-data.frame(  #this represents the min and max values used when randomly assigning initial values to the population
  daymin=0,daymax=100,
  tempmin=0,tempmax=10,
  precipmin=0,precipmax=10)
sds<-data.frame( #standard deviations for trait mutations. Currently set so that variance = max initial trait value
  day=sqrt(start$daymax),
  temp=sqrt(start$tempmax),
  precip=sqrt(start$precipmax))
mutrate<-data.frame( #probability of each trait mutating in an individual. Mutations are independent of one another
  const=.5,
  day=.5,
  temp=.5,
  precip=.5)

#######################################
# Handling libraries and source files #
#######################################

#libraries
library(timeDate)
library(Cairo) #I'm not sure if we need this with the plotting removed

#Set appropriate working directory
if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin"){ #If it's collin
  if(Sys.info()[1]=="Linux"){
    setwd("/home/collin/Dropbox/Grad school/research projects/yang_cue")
  }else{
    if(Sys.info()["nodename"]=="DESKTOP-D6QSU8F"){
      setwd("G:\\Repos\\phenology-cues") #desktop
    }else{
    setwd("C:\\Repos\\phenology-cues") #desktop
    }
  }
}else{
  if(Sys.getenv("COMPUTERNAME")=="ENT-YANG01"){
    setwd("C:\\Users\\lhyang.ent-yang01\\SkyDrive\\Phenology simulation\\phenology-cues")#desktop
  }else{
    setwd("C:\\Users\\lhyang\\Skydrive\\Phenology simulation\\phenology-cues")} #laptop
}
#Load sources file(s)
source("windows_subs.R")



###############################
# Generate environmental data #
###############################


if(runType=="standard"){
  out=yeargen.davistest(numYears,best.temp = best.temp,sd.temp = sd.temp,best.precip = best.precip,sd.precip = sd.precip)
  years.list=out[["years.list"]]
  years.index=out[["years.index"]]
} else if(runType=="unitTestConst"){
  out=yeargen.const(numYears)
  years.list=out[["years.list"]]
  years.index=out[["years.index"]]
} else if (runType=="unitTestRand"){
  out=yeargen.rand(numYears)
  years.list=out[["years.list"]]
  years.index=out[["years.index"]]
}
#######################
# initializing population
#######################
##intialize a population of N individuals
# Their min and max values are determined by the start$ parameters
b.day<-runif(n=N,min=start$daymin,max=start$daymax)
b.temp<-runif(n=N,min=start$tempmin,max=start$tempmax)
b.precip<-runif(n=N,min=start$precipmin,max=start$precipmax)
newpop<-data.frame(b.day,b.temp,b.precip)
pop<-selection(newpop,duration,year=years.list[[1]],N)

###########################
## Running the Simulation #
###########################
pophistory=runSim(startpop=pop,years.list=years.list,
                  years.ind=years.index,N=N,duration=duration,
                  sds=sds,mutrate=mutrate,generations=length(years.index[-1]))
#Note: we've already used year 1 in initiating the pop

#####################
#Saving our results #
#####################
#Set appropriate working directory
if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin"){ #If it's collin
  if(Sys.info()[1]=="Linux"){
    setwd("/home/collin/Dropbox/Grad school/research projects/yang_cue")
  }else{
    if(Sys.info()["nodename"]=="DESKTOP-D6QSU8F"){
      setwd("G:\\Repos\\phenology-cues") #desktop
    }else{
      setwd("C:\\Repos\\phenology-cues") #desktop
    }
  }
}else{
  if(Sys.getenv("COMPUTERNAME")=="ENT-YANG01"){
    setwd("C:\\Users\\lhyang.ent-yang01\\SkyDrive\\Phenology simulation\\phenology-cues")#desktop
  }else{
    setwd("C:\\Users\\lhyang\\Skydrive\\Phenology simulation\\phenology-cues")} #laptop
}
#We have a "save data" script called windows_save.R
source("windows_save.R")

############
# Plotting #
############

yearFit=NULL
for(i in years.index){
  curfits=years.list[[i]]$fit.daily
  if(length(curfits)==366){curfits=curfits[-366]} #to handle leap years, remove least useful day
  yearFit=rbind(yearFit,curfits)
}
par(mar=c(5,5,4,3))
meanFit=apply(yearFit,2,mean)
meanFitSum=NULL
for(i.day in 1:365){
  meanFitSum=c(meanFitSum,sum(rep(meanFit,2)[i.day:(i.day+duration-1)]))
}

x11(width=9,height=6)
for(curgen in seq(2,length(years.index),length=5)){
  curgen=round(curgen)
  arheight=rep(max(meanFit)*1.1,N)
  emergeDay=pophistory[[curgen]]$emerge
  plot(meanFit,type='l',ylim=c(0,max(meanFit)*1.2))
  arrows(y0=jitter(arheight,factor=1.5),x0=emergeDay,x1=emergeDay+duration-1,length=.1)
  dev.print(pdf,paste("dailyfit-run",runNumber,"-gen",curgen,"-meanfit.pdf",sep=""))
  plot(meanFitSum,type='l',ylim=c(0,max(meanFitSum)*1.2),
       main=paste("Mean fitness gained, gen",curgen),
       ylab="Fitness gained",
       xlab="Julian date",
       cex.lab=1.3,
       cex.main=1.3)
  arheight=jitter(rep(max(meanFitSum)*1.05,N),factor=.8)
  arrows(y0=arheight+.05*max(meanFitSum),x0=emergeDay,y1=arheight,length=.1)
  dev.print(pdf,paste("dailyfitSum-run",runNumber,"-gen",curgen,"-meanfit.pdf",sep=""))

  #now calculate the fitSum for THIS YEAR ONLY
  FitSum=NULL
  for(i.day in 1:365){
    FitSum=c(FitSum,sum(rep(years.list[[curgen]]$fit.daily,2)[i.day:(i.day+duration-1)]))
  }
  plot(FitSum,type='l',ylim=c(0,max(FitSum)*1.2),
       main=paste("Fitness gained this year, gen",curgen),
       ylab="Fitness gained",
       xlab="Julian date",
       cex.lab=1.3,
       cex.main=1.3)
  arheight=jitter(rep(max(FitSum)*1.05,N),factor=.8)
  arrows(y0=arheight+.05*max(FitSum),x0=emergeDay,y1=arheight,length=.1)
  dev.print(pdf,paste("dailyfitSum-run",runNumber,"-gen",curgen,"-actualfit.pdf",sep=""))
}

#Calculating changes in mean fitness through time
meanfit=rep(0,length(years.index))
maxfit=meanfit
for(curgen in 1:length(years.index)){
  meanfit[curgen]=mean(pophistory[[curgen]]$Wi)
  cur.fitness=years.list[[years.index[curgen]]]$fit.daily
  cur.fitness.durated=rep(0,365)
  for(i.day in 1:length(cur.fitness.durated)){cur.fitness.durated[i.day]=sum(cur.fitness[i.day:min(i.day+duration-1,365)])}
  maxfit[curgen]=max(cur.fitness.durated)
}
plot(maxfit,type='l',col='red',
     main=paste("Mean fitness through time for run",runNumber),
     ylab="generation",
     xlab="Raw mean fitness",
     sub="red is maximum possible",
     ylim=c(0,max(maxfit))
)
points(1:length(meanfit),meanfit,type='l',ylim=c(0,.04))
dev.print(pdf,paste("meanfitThroughTime_wMax-run",runNumber,"-gen",curgen,".pdf",sep=""))

plot(meanfit,type="l",
     main=paste("Mean fitness through time for run",runNumber),
     ylab="generation",
     xlab="Raw mean fitness"
)
dev.print(pdf,paste("meanfitThroughTime-run",runNumber,"-gen",curgen,".pdf",sep=""))

#Looking at coef changes through time
#  doing two things. exp.eff is the "expected" effect size by muliplying the coefficient of each individual by the mean environmental conditions for that generation
#  The act.eff is the actual effect size, found by multiplying the coefficient of each indiv by the environmental conditions of their day of emergence.
#    Those plots use crosses to represent individuals who didn't emerge until the final day, and circles for those that emerged on a normal day (ie their cue
exp.eff=meanTraitEff(years.index,years.list,pophistory)
act.eff=actTraitEff(years.index,years.list,pophistory)

x11(width=9,height=6)
par(mar=c(5,5,4,4))

traitplot(indivs=act.eff,trait="b.day")
dev.print(pdf,paste("coefs-day-expected-run",runNumber,".pdf",sep=""))
traitplot(indivs=act.eff,trait="b.temp")
dev.print(pdf,paste("coefs-temp-expected-run",runNumber,".pdf",sep=""))
traitplot(indivs=act.eff,trait="b.precip")
dev.print(pdf,paste("coefs-precip-expected-run",runNumber,".pdf",sep=""))

emergePlot(indivs=act.eff,trait="b.day")
dev.print(pdf,paste("coefs-day-actual-run",runNumber,".pdf",sep=""))
emergePlot(indivs=act.eff,trait="b.temp")
dev.print(pdf,paste("coefs-temp-actual-run",runNumber,".pdf",sep=""))
emergePlot(indivs=act.eff,trait="b.precip")
dev.print(pdf,paste("coefs-precip-actual-run",runNumber,".pdf",sep=""))

#Plot it all in one
x11()
par(mfrow=c(3,1))
traitplot(indivs=act.eff,trait="b.day")
traitplot(indivs=act.eff,trait="b.temp")
traitplot(indivs=act.eff,trait="b.precip")
dev.print(pdf,paste("coefs-all-expected-run",runNumber,".pdf",sep=""))

emergePlot(indivs=act.eff,trait="b.day")
emergePlot(indivs=act.eff,trait="b.temp")
emergePlot(indivs=act.eff,trait="b.precip")
dev.print(pdf,paste("coefs-all-actual-run",runNumber,".pdf",sep=""))