#data input script
#reads climate data, truncates to 101 years, and imputes missing values

#clear all variables
rm(list=ls())

#libraries
library(timeDate)
library(Amelia)
library(parallel)
library(ggplot2)
library(plyr)

# Set working directory ---------------------------------------------------

if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin"){ #If it's collin
  if(Sys.info()[1]=="Linux"){
    setwd("/home/collin/Dropbox/Grad school/research projects/yang_cue")
  }else{
    setwd("C:\\Repos\\phenology-cues") #desktop
  }
}else{
  if(Sys.getenv("COMPUTERNAME")=="ENT-YANG01"){
    setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")#desktop
  }else{  
    setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")} #SP4
}

# input data --------------------------------------------------------------

start.date<-as.Date("1914-01-01") #the dataset will start on 1914-01-01
end.date<-as.Date("2014-12-31") #the dataset will end on 2014-12-31
complete.dates<-data.frame(seq(start.date,end.date,1)) #all dates between the start and end dates
colnames(complete.dates)<-c("COMPLETE.DATES") #renaming the column

#daily<-read.csv("davis-data/626713.csv", header=T, na.strings="-9999")
daily<-read.csv("ithaca-data/627453.csv", header=T, na.strings="-9999")
daily$DATE2<-as.Date(as.character(daily$DATE),format="%Y %m %d") #DATE2 is date formatted
daily<-merge(daily,complete.dates,by.x="DATE2",by.y="COMPLETE.DATES",all.y=TRUE) #add missing rows

daily$YEAR<-as.numeric(format(daily$DATE2,"%Y")) #simple field for year
daily$MONTH<-as.numeric(format(daily$DATE2,"%m")) #simple field for month
daily$DAY<-as.numeric(format(daily$DATE2,"%d")) #simple field for day
daily$DAY.OF.YEAR<-as.POSIXlt(daily$DATE2)$yday+1 #day of the year

daily$JULIAN<-julian(daily$DATE2,origin=start.date-1) #1914-01-01 is day 1
daily<-daily[daily$YEAR>1913 & daily$YEAR<2015,] #truncates the data to 101 complete years between 1914 and 2014

daily$PRCP<-daily$PRCP/10 #precips are reported in tenths of mm, units changed to mm
daily$TMAX<-daily$TMAX/10 #temps are reported in tenths of degree C, units changed to degree C
daily$TMIN<-daily$TMIN/10 #temps are reported in tenths of degree C, units changed to degree C


daily<-daily[,c("DATE2","JULIAN","YEAR","MONTH","DAY","DAY.OF.YEAR","PRCP","TMAX","TMIN")] #simplified dataframe

daily.means<-aggregate(cbind(TMAX,TMIN,PRCP)~DAY.OF.YEAR, data=daily, mean)
colnames(daily.means)<-c("DAY.OF.YEAR","TMAX.means","TMIN.means","PRCP.means")
daily.sd<-aggregate(cbind(TMAX,TMIN,PRCP)~DAY.OF.YEAR, data=daily, sd)
colnames(daily.sd)<-c("DAY.OF.YEAR","TMAX.sd","TMIN.sd","PRCP.sd")

# imputation of missing data ----------------------------------------------

temp.bounds<-matrix(c(8, min(daily$TMAX[!is.na(daily$TMAX)]), max(daily$TMAX[!is.na(daily$TMAX)])), nrow = 1, ncol = 3) #matrix in the format [column lowerbound upperbound]

daily<-merge(daily,daily.means,by="DAY.OF.YEAR")
daily<-merge(daily,daily.sd,by="DAY.OF.YEAR")
daily<-daily[order(daily$DATE2),]

#use priors based on the daily mean and sd for TMAX across the dataset
temp.priors<-cbind(1:max(nrow(daily)),8,daily$TMAX.means,daily$TMAX.sd) #matrix [row col mean sd]

a.out<-amelia(daily,m=1,ts="DAY.OF.YEAR",cs="YEAR",idvars=c("DATE2","MONTH","DAY","JULIAN"),intercs=T,splinetime=3,parallel="snow", npcus=detectCores()-1,bounds=temp.bounds,priors=temp.priors,leads="TMAX",lags="TMAX",max.resample = 500)

daily.imp<-a.out$imputations[[1]]

daily.imp[daily.imp$PRCP<0,"PRCP"]<-0 #set all negative precip values to zero

plot(1:length(daily.imp[daily.imp$YEAR==1918,"TMAX"]),daily.imp[daily.imp$YEAR==1918,"TMAX"])
plot(1:length(daily.imp[daily.imp$YEAR==1924,"TMAX"]),daily.imp[daily.imp$YEAR==1924,"TMAX"])

# descriptive statistics --------------------------------------------------

qplot(x=DAY.OF.YEAR,y=TMAX,data=daily.means)
qplot(x=DAY.OF.YEAR,y=TMIN,data=daily.means)
qplot(x=DAY.OF.YEAR,y=PRCP,data=daily.means)

yearnames<-unique(daily.imp$YEAR) #vector of years
yearlist<-split(daily.imp,daily.imp$YEAR) #list of each year separated
yearvar<-data.frame(row.names=yearnames) #dataframe to hold environmental variability

for (i in 1:length(yearnames)){
  #dataframe to compare with mean conditions
  #this creates a VAR.x for each year and a VAR.y for the daily means
  comparison<-merge(yearlist[[i]],daily.means,by="DAY.OF.YEAR") 
  #sum of squared differences with an average year - how weird is each year?
  yearvar[i,"TMAX.SS"]<-sum((comparison$TMAX.x-comparison$TMAX.y)^2)
  yearvar[i,"TMIN.SS"]<-sum((comparison$TMIN.x-comparison$TMIN.y)^2)
  yearvar[i,"PRCP.SS"]<-sum((comparison$PRCP.x-comparison$PRCP.y)^2)
  #CV within years - how variable is each year?
  yearvar[i,"TMAX.CV"]<-sd(comparison$TMAX.x)/mean(comparison$TMAX.x)
  yearvar[i,"TMIN.CV"]<-sd(comparison$TMIN.x)/mean(comparison$TMIN.x)
  yearvar[i,"PRCP.CV"]<-sd(comparison$PRCP.x)/mean(comparison$PRCP.x)
  #sum of differences (not squared) with an average year - how hot/wet is each year?
  yearvar[i,"TMAX.DEL"]<-sum(comparison$TMAX.x-comparison$TMAX.y)
  yearvar[i,"TMIN.DEL"]<-sum(comparison$TMIN.x-comparison$TMIN.y)
  yearvar[i,"PRCP.DEL"]<-sum(comparison$PRCP.x-comparison$PRCP.y)
  #annual totals
  yearvar[i,"TMAX.TOT"]<-max(cumsum(comparison$TMAX.x))
  yearvar[i,"TMIN.TOT"]<-max(cumsum(comparison$TMIN.x))
  yearvar[i,"PRCP.TOT"]<-max(cumsum(comparison$PRCP.x))
  #spring totals
  yearvar[i,"TMAX.SPR"]<-max(cumsum(comparison$TMAX.x[1:120]))
  yearvar[i,"TMIN.SPR"]<-max(cumsum(comparison$TMIN.x[1:120]))
  yearvar[i,"PRCP.SPR"]<-max(cumsum(comparison$PRCP.x[1:120]))
}  

# generating some environmental histories ---------------------------------

strange.TMAX.25<-as.numeric(rownames(tail(yearvar[order(yearvar$TMAX.SS),],25))) #25 least normal TMAX years
normal.TMAX.25<-as.numeric(rownames(head(yearvar[order(yearvar$TMAX.SS),],25))) #25 most normal TMAX years

strange.PRCP.25<-as.numeric(rownames(tail(yearvar[order(yearvar$PRCP.SS),],25))) #25 least normal PRCP years
normal.PRCP.25<-as.numeric(rownames(head(yearvar[order(yearvar$PRCP.SS),],25))) #25 most normal PRCP years

hot.25<-as.numeric(rownames(tail(yearvar[order(yearvar$TMAX.DEL),],25))) #25 hottest TMAX years
cool.25<-as.numeric(rownames(head(yearvar[order(yearvar$TMAX.DEL),],25))) #25 coolest TMAX years

hot.cold.50<-c(hot.25,cool.25) #mix of hottest and coldest years
goldilocks.50<-as.numeric(rownames(head(yearvar[order(abs(yearvar$TMAX.DEL)),],50))) #most normal TMAX

wet.25<-as.numeric(rownames(tail(yearvar[order(yearvar$PRCP.TOT),]),25)) #25 wettest years
dry.25<-as.numeric(rownames(head(yearvar[order(yearvar$PRCP.TOT),]),25)) #25 dryest years

wet.dry.50<-c(wet.25,dry.25) #wet and dry years
meso.wet.50<-as.numeric(rownames(head(yearvar[order(abs(yearvar$PRCP.DEL)),],50))) #middle wetness years

early.25<-as.numeric(rownames(tail(yearvar[order(yearvar$TMAX.SPR),],25))) #25 years with warm springs
late.25<-as.numeric(rownames(head(yearvar[order(yearvar$TMAX.SPR),],25)))  #25 years with cool springs

early.late.50<-c(early.25,late.25) #early and late years
punctual.50<-as.numeric(rownames(yearvar[order(abs(yearvar$TMAX.SPR)),][26:75,])) #normal springs
