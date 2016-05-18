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
traits=c("day","cutemp","cuprecip")
#traits=c("day","temp","precip","cutemp","cuprecip","daysq","tempsq","precipsq","cutempsq","cuprecipsq")
#unitTestConst is for running the population through a unit test with the same gaussian fitness every year
#and constant environmental conditions
#unitTestRand will be for running the populations through a
#unit test with the same gaussian fitness every year and random envi conditions
#standard is for running the populations through a set of replications of the first 10 good years of the davis data
plotExtra=TRUE # do we plot snapshots of emergence through time?
plotPheno=TRUE # do we plot snapshots of phenotype through time?
runName="-meetingtest" #string without spaces (for simplicity)
duration=10
N=100 #number of individuals
numYears=1000 #number of years to simulate
burnIn=100 #number of years to not plot (to avoid scale issues from broad initial population traits
best.temp=30; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=10; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
start<-list(  #this represents the min and max values used when randomly assigning initial values to the population
  day=c(1,1200),
  temp=c(0,0),
  precip=c(0,0),
  cutemp=c(1,27000),
  cuprecip=c(1,1800),
  daysq=c(0,0),
  tempsq=c(0,0),
  precipsq=c(0,0),
  cutempsq=c(0,0),
  cuprecipsq=c(0,0)
)
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
  temporary[i.trait]=start[i.trait]
}
start=temporary
####THE SDS BELOW NEED TO BE RESCALED FOR FAIRNESS.
sds<-list( #standard deviations for trait mutations.
  day=.1,
  temp=.1,
  precip=.1,
  cutemp=.1,
  cuprecip=.1, #max
  daysq=.1,
  tempsq=.1,
  precipsq=.1,
  cutempsq=.1,
  cuprecipsq=.1)
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



years.index=rep(50:100,length.out=numYears) # This is the list of which year.list data to use for each generation of the model
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
    curvals=100/randnums
  }
  curname=paste("b.",i.trait,sep="")
  assign(curname,curvals)
}
newpop<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
pop<-selection(newpop,duration,year=years.list[[years.index[1]]],N)
###########################
## Running the Simulation #
###########################
pophistory=runSim(startpop=pop,years.list=years.list,
                  years.ind=years.index,N=N,duration=duration,
                  sds=sds,mutrate=mutrate,generations=length(years.index))
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

proc.time()-ptm
