# rm(list=ls())

require(graphics)
require(fields)

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
set_wrkdir()
source("scripts/windows_subs.R")
fitshape="standgauss"
source(paste("fitcurve/",fitshape,".R",sep=""))
load(paste("data-years/","Davis","Dat.Rdata",sep = ""),envir= e <- new.env())
meanyr=e$daily.means
meanyr=meanyr[,-3]
names(meanyr)=c("day","temp","precip")
baseTempQ=10
baseTemp=sort(meanyr$temp)[round(365/baseTempQ)]

yr.list=yeargen.davis(best.precip,sd.precip,best.temp,sd.temp,baseTemp=baseTemp)
yrnum=5
yr=yr.list[[1]][[yrnum]]
# x11()
par(mar=c(5,5,5,3))
myplot<-function(x,y,xlab,ylab,main,col){
  plot(x,y,xlab=xlab,ylab=ylab,main=main,col=col,type='l',cex.lab=1.8,cex.main=1.8,lwd=2,cex.axis=1.5)
}

set_wrkdir()
setwd("documents/symposium 2016/figs")
dates=seq(as.Date(paste(yr.list[[2]][[yrnum]],"/01/01",sep="")),format="%d/%m",by="days",length=length(yr$day))
#precip
myplot(dates,yr$precip,
       xlab="Day",ylab="Precipitation (in milimeters)",
       main=paste("Davis precipitation cue,",yr.list[[2]][[yrnum]]),
       col="blue")
dev.print(pdf,"cue-precip.pdf")

#Temperature
myplot(dates,yr$temp,
       xlab="Day",ylab="Max temperature (C)",
       main=paste("Davis temperature cue,",yr.list[[2]][[yrnum]]),
       col="Red")
dev.print(pdf,"cue-temp.pdf")

#Cutemp
myplot(dates,yr$cutemp,
       xlab="Day",ylab="Cumulative temperature",
       main=paste("Davis cumulative temperature cue,",yr.list[[2]][[yrnum]]),
       col="orange")
dev.print(pdf,"cue-cutemp.pdf")

#Photoperiod
myplot(dates,yr$day,
       xlab="Day",ylab="Day",
       main=paste("Davis day cue,",yr.list[[2]][[yrnum]]),
       col="black")
dev.print(pdf,"cue-photoperiod.pdf")

#Show fitness through time
myplot(x=dates,y=yr$fit.daily,
       xlab="Day",
       ylab="Daily fitness",
       main=paste("Davis daily fitness,",yr.list[[2]][[yrnum]]),
       col="black"
)
dev.print(pdf,"dailyfit.pdf")
dayind=109
duration=10
lifeind=(dayind+1):(dayind+duration)
points(dates[dayind],yr$fit.daily[dayind],cex=1,pch=8,lwd=3,col='red')
dev.print(pdf,"dailyfit-2.pdf")
abline(v=dates[dayind+1],col='red',lwd=2,lty=2)
abline(v=dates[dayind+duration],col='red',lwd=2,lty=2)
points(dates[lifeind],yr$fit.daily[lifeind],col='red',type='l',lwd=3)
dev.print(pdf,"dailyfit-3.pdf")

require(zoo) #for use in fitness calculations
duration=10
fit.tot=c(rollapply(c(yr$fit.daily,rep(0,duration-1)),duration,by=1,sum))
myplot(x=dates,y=fit.tot,
       xlab="Day",
       ylab="Total fitness",
       main=paste("Davis total fitness,",yr.list[[2]][[yrnum]]),
       col="blue"
)
dev.print(pdf,"totfit.pdf")



#simple emergence example, b.temp=30
b.temp=30
dayind=min(which(yr$temp>=b.temp))
myplot<-function(x,y,xlab,ylab,main,col){
  plot(x,y,xlab=xlab,ylab=ylab,main=main,col=col,type='l',cex.lab=1.8,cex.main=1.8,lwd=1,cex.axis=1.5)
}
myplot(dates,yr$temp,
       xlab="Day",ylab="Max temperature (C)",
       main=paste("Davis temperature cue,",yr.list[[2]][[yrnum]]),
       col="Red")
dev.print(pdf,"E-example1-1.pdf")
abline(h=b.temp,col='blue')
dev.print(pdf,"E-example1-2.pdf")
points(dates[dayind],b.temp,cex=1,pch=8,lwd=3,col='blue')
segments(x0=dates[dayind],y0=0,y1=b.temp,col='blue')
dev.print(pdf,"E-example1-3.pdf")

#complex example, b.temp=30,b.day=100

#Show fitness function
numpts=500
tempmat=matrix(seq(0,60,length=numpts),ncol=numpts,nrow=numpts,byrow = TRUE)
precipmat=matrix(seq(0,60,length=numpts),ncol=numpts,nrow=numpts,byrow = FALSE)
tempv=tempmat[1:(numpts^2)]
precipv=precipmat[1:(numpts^2)]
fit=dnorm(tempmat,mean=best.temp,sd=sd.temp)*dnorm(precipmat,mean=best.precip,sd=sd.precip)
# image(x=seq(0,50,length=numpts),y=seq(0,50,length=numpts),z=fit)
image.plot(x=seq(0,60,length=numpts),
           y=seq(0,60,length=numpts),
           z=fit,
           xlab="Precipitation",
           ylab="Temperature",
           main="Fitness gained per day",
           #col=heat.colors,
           cex.lab=1.6,cex.main=1.6,cex.axis=1.5
)
dev.print(pdf,"Fitness-function.pdf")

years.df=NULL
for(i in 1:length(yr.list[[1]])){years.df=rbind(years.df,yr.list[[1]][[i]][1:365,])}
years.df=as.data.frame(years.df)
sampind=1:(365*10)
points(jitter((years.df$precip[sampind])),(years.df$temp[sampind]),lwd=2)
dev.print(pdf,"Fitness-function-wpoints.pdf")

require(ggtern)
set_wrkdir()
if(FALSE){
store.endeff=NULL
  for(i in 1:10){
    load(paste("results/EEB d-t-p/resRunEEB d-t-p",i,"/dat.RData",sep=""),envir=te<-new.env(),)
    numgen=length(te$years.index)
    act.eff=actTraitEff(te$years.index,te$years.list,te$pophistory,te$N,te$traits)
    res=act.eff[act.eff[,"gen"]==numgen,]
    res=cbind(res,sim=rep(i,te$N))
    store.endeff=rbind(store.endeff,res)
  }
store.endeff=as.data.frame(store.endeff)
save.image("rawdump.Rdat")
}


ggtern(data=store.endeff,aes(b.day,b.temp,b.precip)) +
  geom_point(fill=factor(store.endeff$sim),shape=21,size=3)+
  theme(text=element_text(size=18))
#  theme(plot.margin=unit(c(2,2,2,2),"cm"))

setwd("documents/symposium 2016/figs")
dev.print(pdf,"coeff-eff-tern.pdf")

#Plotting meanyear
myplot(dates,meanYr$temp,
       xlab='Dates',ylab='Temperature',main="Baseline year",col='black')
dev.print(pdf,"basyear.pdf")
#adding day to day noise
dayvar=meanYr;dayvar$temp=dayvar$temp+rnorm(365,sd=1.5)
myplot(dates,meanYr$temp,
       xlab='Dates',ylab='Temperature',main="+day-to-day noise",col='grey')
points(dates,dayvar$temp,
       xlab='Dates',ylab='Temperature',main="Baseline year",col='red',type='l')
dev.print(pdf,"basyear-plus-d2d.pdf")
#adding phenology shift
yrvar=dayvar
yrvar$temp= dayvar$temp[((1:365+20) %% 365)+1]
myplot(dates,meanYr$temp,
       xlab='Dates',ylab='Temperature',main="+seasonal shift",col='grey')
points(dates,dayvar$temp,
       xlab='Dates',ylab='Temperature',main="Baseline year",col='grey',type='l',lwd=1)
points(dates,yrvar$temp,
       xlab='Dates',ylab='Temperature',main="Baseline year",col='red',type='l',lwd=1)
dev.print(pdf,"basyear-plus-all.pdf")

#Working with fig1 stuff
set_wrkdir()
load("results/fig1/test-1000yr/fig1dat-versiontest-1000yr.Rdata")
varmeans=overall.res
varpts=unique(varmeans[,c("daystd","yearstd")])
winner=rep("init",nrow(varpts))
yearstds=daystds=rep(-10,nrow(varpts))
for(i.pt in 1:nrow(varpts)){
  cur.ind=which(varmeans$daystd==varpts$daystd[i.pt] & varmeans$yearstd==varpts$yearstd[i.pt])
  curvar=varmeans[cur.ind,]
  if(sum(curvar$geofit==min(curvar$geofit))>1){
    winner[i.pt]="tie"
  }else{
    winner[i.pt]=curvar[which.max(curvar$geofit),]$trait
  }
}
df=data.frame(x=varpts$daystd,
              y=varpts$yearstd,
              z=winner)
x11()
(ggplot(df, aes(x,y))+
  geom_tile(aes(fill=z))+
  labs(title="Winning trait",
       x="day to day variance",
       y="year to year variance")+
  theme(text=element_text(size=20)))
setwd("documents/symposium 2016/figs")
dev.print(pdf,"winners.pdf")

#get distances between each y points
ddiff=mean(diff(unique(varmeans$daystd)))
yrdiff=mean(diff(unique(varmeans$yearstd)))
nyrsd=varmeans$yearstd
varmeans$geofit=varmeans$geofit/numYears
nyrsd=nyrsd-yrdiff/5*as.numeric(varmeans$trait=='temp')+yrdiff/5*as.numeric(varmeans$trait=='cutemp')
ndysd=varmeans$daystd+ddiff/5*as.numeric(varmeans$trait=='day')
size=(varmeans$geofit-min(varmeans$geofit)+.025)*15

# x11()
# par(fig=c(0.1,1,0.1,0.8), new=TRUE)
temp.vals=unique(factor(varmeans$trait))
temp.names=unique(as.character(varmeans$trait))
par(mar=c(5.1,5.1,4.1,8.1),xpd=TRUE)
plot(x=ndysd,y=nyrsd,col=factor(varmeans$trait),cex=size,pch=19,
     xlab="Day to day variation",
     ylab="Year to year variation",
     main="Performance of traits",
     cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
legend("topright",inset=c(-.2,0),legend=temp.names,fill=factor(varmeans$trait),cex=1.4,bty="n")
legend(x=1,legend=temp.names,fill=factor(varmeans$trait))
dev.print(pdf,"trait-response-surface.pdf")
