#data input script
#reads climate data, truncates to 101 years, and imputes missing values

#clear all variables
rm(list=ls())

#libraries
library(timeDate)
library(Amelia)
library(parallel)
library(ggplot2)

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

daily<-read.csv("davis-data/626713.csv", header=T, na.strings="-9999")
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


# imputation of missing data ----------------------------------------------

t.go<-proc.time()
a.out<-amelia(daily,m=5,ts="DAY.OF.YEAR",cs="YEAR",idvars=c("DATE2","MONTH","DAY","JULIAN"),intercs=T,splinetime=3,parallel = "snow", npcus=detectCores())
proc.time()-t.go

plot(a.out)

save(a.out, file = "imputations.RData") 



# descriptive statistics --------------------------------------------------

daily.means<-aggregate(cbind(TMAX,TMIN,PRCP)~DAY.OF.YEAR, data=daily, mean)

qplot(x=DAY.OF.YEAR,y=TMAX,data=daily.means)
qplot(x=DAY.OF.YEAR,y=TMIN,data=daily.means)
qplot(x=DAY.OF.YEAR,y=PRCP,data=daily.means)











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



#this is resource intensive - USE WITH CAUTION
overimpute(a.out, var = "TMAX")