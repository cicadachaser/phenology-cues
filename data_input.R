#data input script
#reads climate data, truncates to 101 years, and imputes missing values

#clear all variables
rm(list=ls())

#libraries
library(timeDate)
library(Amelia)

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
    setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")} #SP4
}

#input data
davis.daily<-read.csv("davis-data/626713.csv", header=T, na.strings="-9999")
davis.daily$YEAR<-as.numeric(substr(davis.daily$DATE,1,4)) #simple field for year
davis.daily<-davis.daily[davis.daily$YEAR>1913 & davis.daily$YEAR<2015,] #truncates the data to 101 complete years between 1914 and 2014

davis.daily$PRCP<-davis.daily$PRCP/10 #precips are reported in tenths of mm
davis.daily$TMAX<-davis.daily$TMAX/10 #temps are reported in tenths of degree C
davis.daily$TMIN<-davis.daily$TMIN/10 #temps are reported in tenths of degree C
davis.daily$DATE2<-as.Date(as.character(davis.daily$DATE),format="%Y %m %d") #DATE2 is date formatted
davis.daily$JULIAN<-julian(davis.daily$DATE2,origin=as.Date("1913-12-31")) #1914-01-01 is day 1
davis.daily$MONTH<-as.numeric(substr(davis.daily$DATE,5,6)) #simple field for month
davis.daily$DAY<-as.numeric(substr(davis.daily$DATE,7,8)) #simple field for day
davis.daily<-davis.daily[,c("DATE2","JULIAN", "YEAR","MONTH","DAY","PRCP","TMAX","TMIN")] #simplified dataframe

davis.yearlist<-split(davis.daily,davis.daily$YEAR) #list of each year separated

davis.yearnames<-unique(davis.daily$YEAR) #gives a list of all the years in the data

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

#temporary plots checking that the data are mostly complete
plot(davis.daily$JULIAN,davis.daily$TMAX) 
plot(davis.daily$JULIAN,is.na(davis.daily$TMAX))
#however, there are 143 missing rows
max(davis.daily$JULIAN)-length(davis.daily$JULIAN)


match(seq(1:max(davis.daily$JULIAN)),davis.daily$JULIAN)




#amelia scripting in progress
amelia(davis.daily,m=5,idvars=c("DATE2","JULIAN",)

