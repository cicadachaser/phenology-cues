#Script for carrying out simulations of evolution of emergence decisions
#Uses windows_subs.R, windows_save.R, and windows_plot.R

#clear all variables
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

ptm <-proc.time()
#########################
# Simulation parameters #
#########################
runType="standard" ##THIS DETERMINES WHAT KIND OF YEARS WE'RE USING!
traits=c("day","temp","precip")
numsims=10 # number of simulations of each type to do
runsnames=c("-daytempprecip-","-day-") #string without spaces (for simplicity)
#traits=c("day","temp","precip","cutemp","cuprecip","daysq","tempsq","precipsq","cutempsq","cuprecipsq")
#unitTestConst is for running the population through a unit test with the same gaussian fitness every year
#and constant environmental conditions
#unitTestRand will be for running the populations through a
#unit test with the same gaussian fitness every year and random envi conditions
#standard is for running the populations through a set of replications of the first 10 good years of the davis data
plotExtra=FALSE # do we plot snapshots of emergence through time?
plotPheno=TRUE # do we plot snapshots of phenotype through time?
viewLength=200 #for comparisons of simulation types,
#  how many generations (starting from the final and working backwards) to plot/compare
duration=10 #number of days organizm is emerged.
N=100 #number of individuals
numYears=5000 #number of years to simulate
burnIn=200 #number of years to not plot (to avoid scale issues from broad initial population traits)
best.temp=30; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=10; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
mutdist=.01 #What fraction of the total "cue space" should mutations (on average) traverse (kinda).
# if we're looking at b.day, the cue space is 365 days. So a mutdist of .01 means that each mutation will
# be drawn from a normal distribution with mean 0 and sd 3.65
maxcues<-list(#these are approximately the max values of the first five years of the Davis data.
  # use them for scaling the starting values, and for
  day=365,
  temp=45,
  precip=100,
  cutemp=9000,
  cuprecip=600,
  daysq=133225,
  tempsq=1971.36,
  precipsq=10000,
  cutempsq=250000,
  cuprecipsq=20000
)
start<-list(#these are used to generate the starting values of individuals. Starting values will produce
  # individuals who, if they relied solely on a single cue, would emerge uniformly across a range of
  # cues from 1 to x times approximately the max value, where x is the number of traits used.
  # The max value is selected to be approximately the maximum found in the first five years of the davis
  # data set, and the min is set to approximately the minimum found. Note that we can't use zero.
  day=c(1,maxcues$day*length(traits)),
  temp=c(5,maxcues$temp*length(traits)),
  precip=c(.01,maxcues$precip*length(traits)),
  cutemp=c(5,maxcues$cutemp*length(traits)),
  cuprecip=c(.01,maxcues$cuprecip*length(traits)),
  daysq=c(1,maxcues$daysq*length(traits)),
  tempsq=c(5,maxcues$tempsq*length(traits)),
  precipsq=c(.01,maxcues$precipsq*length(traits)),
  cutempsq=c(5,maxcues$cutempsq*length(traits)),
  cuprecipsq=c(.01,maxcues$cuprecipsq*length(traits))
)
temporary<-list( #temporary object for masking unused traits
  day=0,
  temp=0,
  precip=0,
  cutemp=0,
  cuprecip=0,
  daysq=0,
  tempsq=0,
  precipsq=0,
  cutempsq=0,
  cuprecipsq=0
)
for(i.trait in traits){
  temporary[i.trait]=start[i.trait]
}
start=temporary
####THE SDS BELOW NEED TO BE RESCALED FOR FAIRNESS.
sds<-list( #standard deviations for trait mutations.
  day=mutdist*maxcues$day,
  temp=mutdist*maxcues$temp,
  precip=mutdist*maxcues$precip,
  cutemp=mutdist*maxcues$cutemp,
  cuprecip=mutdist*maxcues$cuprecip, #max
  daysq=mutdist*maxcues$daysq,
  tempsq=mutdist*maxcues$tempsq,
  precipsq=mutdist*maxcues$precipsq,
  cutempsq=mutdist*maxcues$cutempsq,
  cuprecipsq=mutdist*maxcues$cuprecipsq)
mutrate<-list( #probability of each trait mutating in an individual. Mutations are independent of one another
  day=.01,
  temp=.01,
  precip=.01,
  cutemp=.01,
  cuprecip=.01,
  daysq=.01,
  tempsq=.01,
  precipsq=.01,
  cutempsq=.01,
  cuprecipsq=.01)
#Now ensure that mutrate is zero for any trait that ISN'T in the trait list
temporary<-list( #probability of each trait mutating in an individual. Mutations are independent of one another
  day=0,
  temp=0,
  precip=0,
  cutemp=0,
  cuprecip=0,
  daysq=0,
  tempsq=0,
  precipsq=0,
  cutempsq=0,
  cuprecipsq=0
)
for(i.trait in traits){
  temporary[i.trait]=mutrate[i.trait]
}
mutrate=temporary

years.indmat=matrix(sample(1:100,size=numYears*numsims,replace=TRUE),ncol=numYears,nrow=numsims) # This is the list of which year.list data to use for each generation of the model
#######################################
# Handling libraries and source files #
#######################################

#libraries
library(timeDate)
# library(Cairo) #I'm not sure if we need this with the plotting removed
library(zoo)
#Set appropriate working directory
set_wrkdir()
#Load sources file(s)
source("scripts/windows_subs.R")

###############################
# Generate environmental data #
###############################
#Based on the value of "runType", generate the appropriate type of data.
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

store.mean=store.max=matrix(0,nrow=numsims*length(runsnames),ncol=numYears)
  #matrices for storing mean and maximum possible fitness
store.names=rep(0,numsims*length(runsnames)) #vector for storing run names, corresponds to rows of store.mean
finalpops=NULL #for storing the final populations of each run.

count=1 #for tracking year in years.indmat
for(i.sim in 1:numsims){
  years.index=years.indmat[count,]
  runName=sprintf("%s%d",runsnames[1],i.sim)
  #######################
  # initializing population
  #######################
  ##intialize a population of N individuals
  # Their min and max values are determined by the start$ parameters
  #Start by setting them all equal to zero. Fill in with the working traits from the traits variable
  b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,N)
  for(i.trait in traits){
    if(start[[i.trait]][1]==0 & start[[i.trait]][2]==0){
      curvals=rep(0,N)
    }else{
      randnums=runif(n=N,min=start[[i.trait]][1],max=start[[i.trait]][2])
      randnums[randnums==0]=1/(10^10)
      curvals=randnums
    }
    curname=paste("b.",i.trait,sep="")
    assign(curname,curvals)
  }
  newpop<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
  pop<-selection(newpop,duration,year=years.list[[years.index[1]]],N,traits=traits)
  ###########################
  ## Running the Simulation #
  ###########################
  pophistory=runSim(startpop=pop,years.list=years.list,
                    years.ind=years.index,N=N,duration=duration,
                    sds=sds,mutrate=mutrate,generations=numYears,traits=traits)
  #Note: we've already used year 1 in initiating the pop

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

#################################################################################################
#################################################################################################
#################################################################################################
######################################
#HERE WE MAKE CHANGES FOR RUN TYPE 2!#
######################################
traits=c("day")
plotPheno=FALSE

#######################
# Reworking prep work #
#######################
maxcues<-list(#these are approximately the max values of the first five years of the Davis data.
    # use them for scaling the starting values, and for
  day=365,
  temp=45,
  precip=100,
  cutemp=9000,
  cuprecip=600,
  daysq=133225,
  tempsq=1971.36,
  precipsq=10000,
  cutempsq=250000,
  cuprecipsq=20000
)
start<-list(#these are used to generate the starting values of individuals. Starting values will produce
  # individuals who, if they relied solely on a single cue, would emerge uniformly across a range of
  # cues from 1 to x times approximately the max value, where x is the number of traits used.
  # The max value is selected to be approximately the maximum found in the first five years of the davis
  # data set, and the min is set to approximately the minimum found. Note that we can't use zero.
  day=c(1,maxcues$day*length(traits)),
  temp=c(5,maxcues$temp*length(traits)),
  precip=c(.01,maxcues$precip*length(traits)),
  cutemp=c(5,maxcues$cutemp*length(traits)),
  cuprecip=c(.01,maxcues$cuprecip*length(traits)),
  daysq=c(1,maxcues$daysq*length(traits)),
  tempsq=c(5,maxcues$tempsq*length(traits)),
  precipsq=c(.01,maxcues$precipsq*length(traits)),
  cutempsq=c(5,maxcues$cutempsq*length(traits)),
  cuprecipsq=c(.01,maxcues$cuprecipsq*length(traits))
)
temporary<-list( #temporary object for masking unused traits
  day=0,
  temp=0,
  precip=0,
  cutemp=0,
  cuprecip=0,
  daysq=0,
  tempsq=0,
  precipsq=0,
  cutempsq=0,
  cuprecipsq=0
)
for(i.trait in traits){
  temporary[i.trait]=start[i.trait]
}
start=temporary
####THE SDS BELOW NEED TO BE RESCALED FOR FAIRNESS.
sds<-list( #standard deviations for trait mutations.
  day=mutdist*maxcues$day,
  temp=mutdist*maxcues$temp,
  precip=mutdist*maxcues$precip,
  cutemp=mutdist*maxcues$cutemp,
  cuprecip=mutdist*maxcues$cuprecip, #max
  daysq=mutdist*maxcues$daysq,
  tempsq=mutdist*maxcues$tempsq,
  precipsq=mutdist*maxcues$precipsq,
  cutempsq=mutdist*maxcues$cutempsq,
  cuprecipsq=mutdist*maxcues$cuprecipsq)
mutrate<-list( #probability of each trait mutating in an individual. Mutations are independent of one another
  day=.01,
  temp=.01,
  precip=.01,
  cutemp=.01,
  cuprecip=.01,
  daysq=.01,
  tempsq=.01,
  precipsq=.01,
  cutempsq=.01,
  cuprecipsq=.01)
#Now ensure that mutrate is zero for any trait that ISN'T in the trait list
temporary<-list( #probability of each trait mutating in an individual. Mutations are independent of one another
  day=0,
  temp=0,
  precip=0,
  cutemp=0,
  cuprecip=0,
  daysq=0,
  tempsq=0,
  precipsq=0,
  cutempsq=0,
  cuprecipsq=0
)
for(i.trait in traits){
  temporary[i.trait]=mutrate[i.trait]
}
mutrate=temporary

########################################################################################

#Okay, run the second run type
count=1 #for tracking row in years.indmat
for(i.sim in (numsims+1):(2*numsims)){
  years.index=years.indmat[count,]
  runName=sprintf("%s%d",runsnames[2],i.sim-numsims)
  #######################
  # initializing population
  #######################
  ##intialize a population of N individuals
  # Their min and max values are determined by the start$ parameters
  #Start by setting them all equal to zero. Fill in with the working traits from the traits variable
  b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,N)
  for(i.trait in traits){
    if(start[[i.trait]][1]==0 & start[[i.trait]][2]==0){
      curvals=rep(0,N)
    }else{
      randnums=runif(n=N,min=start[[i.trait]][1],max=start[[i.trait]][2])
      randnums[randnums==0]=1/(10^10)
      curvals=randnums
    }
    curname=paste("b.",i.trait,sep="")
    assign(curname,curvals)
  }
  newpop<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
  pop<-selection(newpop,duration,year=years.list[[years.index[1]]],N,traits=traits)
  ###########################
  ## Running the Simulation #
  ###########################
  pophistory=runSim(startpop=pop,years.list=years.list,
                    years.ind=years.index,N=N,duration=duration,
                    sds=sds,mutrate=mutrate,generations=numYears,traits=traits)
  #Note: we've already used year 1 in initiating the pop

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
setwd("results")
resultsdir=paste("compare",runsnames[1],"vs",runsnames[2])
dir.create(resultsdir,showWarnings = FALSE)
setwd(resultsdir)

cols=c('red','blue')
matplot(((numYears-viewLength+1):numYears),t(store.max[,-(1:(numYears-viewLength))]),type='l',col=c(rep('red',numsims),rep('blue',numsims)),
     main=paste("Mean fitness through time for all runs",runsnames[1],runsnames[2]),
     xlab="generation",
     ylab="Raw mean fitness",
     ylim=c(0,max(store.max))
)
matpoints(((numYears-viewLength+1):numYears),t(store.mean[,-(1:(numYears-viewLength))]),type='l',col=c(rep("chocolate",numsims),rep("cornflowerblue",numsims)))
legend(x="bottomright",legend=c(sprintf("max %s",runsnames),runsnames),
  fill=c(cols,"chocolate","cornflowerblue"),cex=2)
dev.print(pdf,paste("compare-allruns.pdf",sep=""))

matplot(((numYears-viewLength+1):numYears),t((store.mean/store.max)[,-(1:(numYears-viewLength))]),type='l',col=c(rep('red',numsims),rep('blue',numsims)),
        main=paste("Scaled mean fitness through time for all runs",runsnames[1],runsnames[2]),
        xlab="generation",
        ylab="Raw mean fitness",
        ylim=c(0,1)
)
abline(h=1)
legend(x="bottomright",legend=runsnames,fill=cols,cex=2)
dev.print(pdf,paste("compare-allruns-scaled.pdf",sep=""))

meanmeans.1=apply(store.mean[1:numsims,],2,mean)
meanmax.1=apply(store.max[1:numsims,],2,mean)
meanmeans.2=apply(store.mean[(1+numsims):(2*numsims),],2,mean)
meanmax.2=apply(store.max[(1+numsims):(2*numsims),],2,mean)

matplot(t(rbind(meanmeans.1,meanmeans.2)[,-(1:(numYears-viewLength))]),type='l',col=cols,
        main=paste("Mean fitness through time for mean of runs",runsnames[1],runsnames[2]),
        xlab="generation",
        ylab="Raw mean fitness"
)

plot((meanmeans.1-meanmeans.2)[-(1:(numYears-viewLength))],type='l',col=cols,
        main=paste("Fitness of ",runsnames[1],"minus",runsnames[2]),
        xlab="generation",
        ylab="Raw mean fitness"
)
abline(h=0)
legend(x="bottomright",legend=runsnames,fill=cols,cex=2)
dev.print(pdf,paste("difference-of-means.pdf",sep=""))

hist((meanmeans.1-meanmeans.2)[-(1:(numYears-viewLength))],breaks=30,
     main="histogram of differences",
     xlab="scaled difference",
     sub="green is mean")
abline(v=0,col='red',lwd=2)
abline(v=mean((meanmeans.1-meanmeans.2)[-(1:(numYears-viewLength))]),lwd=2,col='green')
dev.print(pdf,paste("hist-difference-of-means.pdf",sep=""))

plot((meanmeans.1/meanmax.1-meanmeans.2/meanmax.2)[-(1:(numYears-viewLength))],type='l',col=cols,
     main=paste("Fitness of scaled",runsnames[1],"minus scaled",runsnames[2]),
     xlab="generation",
     ylab="Raw mean fitness"
)
abline(h=0)
legend(x="bottomright",legend=runsnames,fill=cols,cex=2)
dev.print(pdf,paste("difference-of-means-scaled.pdf",sep=""))

hist((meanmeans.1/meanmax.1-meanmeans.2/meanmax.2)[-(1:(numYears-viewLength))],breaks=30,
     main="histogram of scaled differences",
     xlab="scaled difference",
     sub="green is mean");
abline(v=0,col='red',lwd=2)
abline(v=mean((meanmeans.1/meanmax.1-meanmeans.2/meanmax.2)[-(1:(numYears-viewLength))]),lwd=2,col='green')
dev.print(pdf,paste("hist-difference-of-means-scaled.pdf",sep=""))

matplot(t(rbind(meanmeans.1/meanmax.1,meanmeans.2/meanmax.2)[,-(1:(numYears-viewLength))]),type='l',col=cols,
        main=paste("Scaled fitness through time for mean of runs",runsnames[1],runsnames[2]),
        xlab="generation",
        ylab="Raw mean fitness",
        ylim=c(0,1)
)
abline(h=1)
legend(x="bottomright",legend=runsnames,fill=cols,cex=2)
dev.print(pdf,paste("compare-means-scaled.pdf",sep=""))

latefit=apply(store.mean[,(numYears-viewLength):numYears],1,sum)
latemax=apply(store.max[,(numYears-viewLength):numYears],1,sum)
plot(jitter(c(rep(1,numsims),rep(2,numsims)),factor=.1),latefit/latemax,xlim=c(.5,2.5),ylim=c(0,1),
     xaxt='n',
     main=paste("Comparing scaled sum fitness over final",viewLength, "years")
     )
axis(1,at=c(1,2),labels = runsnames)
abline(h=1,col='red')
dev.print(pdf,paste("latefitness-scaled.pdf",sep=""))

