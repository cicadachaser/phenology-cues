#CURRENTLY IN THE "SLOW AND CORRECT" STAGE
#optimize AFTER we confirm it works

#clear all variables
rm(list=ls())

#libraries
library(timeDate)
library(Cairo) #I'm not sure if we need this with the plotting removed

#Set appropriate working directory
if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin"){ #If it's collin
  if(Sys.info()[1]=="Linux"){
    setwd("/home/collin/Dropbox/Grad school/research projects/yang_cue")
  }else{
    setwd("C:\\Repos\\phenology-cues") #desktop
  }
}else{
  if(Sys.getenv("COMPUTERNAME")=="ENT-YANG01"){
    setwd("C:\\Users\\lhyang.ent-yang01\\SkyDrive\\Phenology simulation\\phenology-cues")#desktop
  }else{  
    setwd("C:\\Users\\lhyang\\Skydrive\\Phenology simulation\\phenology-cues")} #laptop
}
#Load sources file(s)
source("windows_subs.R")

#########################
# Simulation parameters #
#########################
generations=24
duration=10
best.temp=15; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=55; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
N=40 #number of individuals 
start<-data.frame(  #this represents the min and max values used when randomly assigning initial values to the population 
  constmin=0,constmax=50,
  daymin=0,daymax=5,
  tempmin=0,tempmax=20,
  precipmin=0,precipmax=5) 
sds<-data.frame( #standard deviations for trait mutations. Currently set so that variance = max initial trait value
  const=sqrt(start$constmax),
  day=sqrt(start$daymax),
  temp=sqrt(start$tempmax),
  precip=sqrt(start$precipmax))
mutrate<-data.frame( #probability of each trait mutating in an individual. Mutations are independent of one another
  const=.1,
  day=.1,
  temp=.1,
  precip=.1)

#input data
davis.daily<-read.csv("davis-data/626713.csv", header=T, na.strings="-9999")
davis.daily$PRCP<-davis.daily$PRCP/10 #precips are reported in tenths of mm
davis.daily$TMAX<-davis.daily$TMAX/10 #temps are reported in tenths of degree C
davis.daily$TMIN<-davis.daily$TMIN/10 #temps are reported in tenths of degree C
davis.daily$DATE2<-as.Date(as.character(davis.daily$DATE),format="%Y %m %d") #DATE2 is date formatted
davis.daily$JULIAN<-julian(davis.daily$DATE2,origin=as.Date("1892-12-31")) #1893-01-01 is day 1...
davis.daily$YEAR<-as.numeric(substr(davis.daily$DATE,1,4)) #simple field for year
davis.daily$MONTH<-as.numeric(substr(davis.daily$DATE,5,6)) #simple field for month
davis.daily$DAY<-as.numeric(substr(davis.daily$DATE,7,8)) #simple field for day
davis.daily<-davis.daily[,c("DATE2","JULIAN", "YEAR","MONTH","DAY","PRCP","TMAX","TMIN")] #simplified dataframe

davis.yearlist<-split(davis.daily,davis.daily$YEAR) #list of each year separated
goodyears=NULL
for(iyear in davis.yearnames){
  nacount=sum(sum(is.na(davis.yearlist[[as.character(iyear)]])))
  daycount=dim(davis.yearlist[[as.character(iyear)]])[1]
  if(nacount==0 & daycount>364){goodyears=c(goodyears,iyear)}
}
davis.yearlist=davis.yearlist[[as.character(goodyears)]]

davis.yearnames<-goodyears #gives a list of all the years in the data

#calculates the "day of year", i.e. Jan 1 is 1, and 12/31 is 365 
#adds a DAY.OF.YEAR column to each dataframe in the year list
for (i in 1:length(davis.yearnames)){
  davis.yearlist[[i]]$DAY.OF.YEAR<-julian(davis.yearlist[[i]]$DATE2, origin=as.Date(paste(davis.yearnames[i],"01","01",sep="-")))+1 #add +1 so that the first day of the year is 1, not zero. 
}

davis.daily<-unsplit(davis.yearlist,davis.daily$YEAR)
davis.daily.means<-aggregate(cbind(TMAX,TMIN,PRCP)~DAY.OF.YEAR, data=davis.daily, mean)

davis.yearvar<-data.frame(row.names=davis.yearnames) #dataframe to hold environmental variability

for (i in 1:length(davis.yearnames)){
  #temporary dataframe to compare with mean conditions
  #this creates a VAR.x for each year and a VAR.y for the daily means
  comparison<-merge(davis.yearlist[[i]],davis.daily.means,by="DAY.OF.YEAR") 
  #number of complete cases (is.na=F) for each year
  davis.yearvar[i,"TMAX.N"]<-sum(complete.cases(davis.yearlist[[i]]$TMAX))
  davis.yearvar[i,"TMIN.N"]<-sum(complete.cases(davis.yearlist[[i]]$TMIN))
  davis.yearvar[i,"PRCP.N"]<-sum(complete.cases(davis.yearlist[[i]]$PRCP))
  #sum of squared differences with an average year - how weird is each year?
  #some years have incomplete data, so this is the mean SS per observed day
  davis.yearvar[i,"TMAX.SS"]<-(sum(comparison$TMAX.x-comparison$TMAX.y,na.rm=T)^2)/davis.yearvar[i,"TMAX.N"]
  davis.yearvar[i,"TMIN.SS"]<-(sum(comparison$TMIN.x-comparison$TMIN.y,na.rm=T)^2)/davis.yearvar[i,"TMIN.N"]
  davis.yearvar[i,"PRCP.SS"]<-(sum(comparison$PRCP.x-comparison$PRCP.y,na.rm=T)^2)/davis.yearvar[i,"PRCP.N"]
  #CV within years - how variable is each year?
  davis.yearvar[i,"TMAX.CV"]<-sd(comparison$TMAX.x,na.rm=T)/mean(comparison$TMAX.x,na.rm=T)
  davis.yearvar[i,"TMIN.CV"]<-sd(comparison$TMIN.x,na.rm=T)/mean(comparison$TMIN.x,na.rm=T)
  davis.yearvar[i,"PRCP.CV"]<-sd(comparison$PRCP.x,na.rm=T)/mean(comparison$PRCP.x,na.rm=T)
  #sum of differences (not squared) with an average year - how hot/wet is each year?
  #some years have incomplete data, so this is the mean difference per observed day
  davis.yearvar[i,"TMAX.DEL"]<-sum(comparison$TMAX.x-comparison$TMAX.y,na.rm=T)/davis.yearvar[i,"TMAX.N"]
  davis.yearvar[i,"TMIN.DEL"]<-sum(comparison$TMIN.x-comparison$TMIN.y,na.rm=T)/davis.yearvar[i,"TMIN.N"]
  davis.yearvar[i,"PRCP.DEL"]<-sum(comparison$PRCP.x-comparison$PRCP.y,na.rm=T)/davis.yearvar[i,"PRCP.N"]
}

#dat<-read.csv("davis.csv", header=T) #this is just a placeholder of monthly climate data for now

### End revisions 11/6/2015 ###

######################################################
# Import sequence of years - LOUIE'S STUFF GOES HERE #
######################################################
years.list=NULL #Replace this with code to grab a list of data frames. Each data frame is a year.
# Each year data frame has $day, $precip, $tmean, $tmax, $tmin
# This will be the same list for all configurations of years - this is essentially just our year database
years.index=NULL # This is the list of which year.list data to use for each generation of the model



######################
# Fitness generation #
######################
#For now, daily incremental fitness will be found by multiplying two gaussian functions together:
#  one for temp, that's maximized at best.temp with sd tempsd
#  the other for precip that's maximized at best.precip with sd precipsd
# We will then normalize the results to vary from 0 to 1
for(i.year in 1:length(years.list)){
  daily.fit=dnorm(years.list[[i.year]]$tmax,mean=best.temp,sd=sd.temp)*dnorm(years.list[[i.year]]$precip,mean=best.precip,sd=sd.precip)
  daily.fit=(daily.fit-min(daily.fit))/(max(daily.fit)-min(daily.fit))
  years.list[[i.year]]=cbind(years.list[[i.year]], fit.daily=daily.fit)
}

#######################
# initializing population
#######################
##intialize a population of N individuals
# Their min and max values are determined by the start$ parameters
b.const<-runif(n=N,min=start$constmin,max=start$constmax)
b.day<-runif(n=N,min=start$daymin,max=start$daymax)
b.temp<-runif(n=N,min=start$tempmin,max=start$tempmax)
b.precip<-runif(n=N,min=start$precipmin,max=start$precipmax)
pop<-data.frame(b.const,b.day,b.temp,b.precip)
pop<-selection(newpop,duration,cur.year,N)

## Run Simulation
pophistory=runSim(startpop=pop,years.list=years.list,
                  years.ind=years.ind,N=N,duration=duration,
                  sds=sds,mutrate=mutrate,generations=generations)