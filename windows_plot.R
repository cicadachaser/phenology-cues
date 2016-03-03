setwd("results")
resultsdir=sprintf("resRun%s",runName)
dir.create(resultsdir,showWarnings = FALSE)
setwd(resultsdir)

yearFit=NULL
for(i in years.index){
  curfits=years.list[[i]]$fit.daily
  if(length(curfits)==366){curfits=curfits[-366]} #to handle leap years, remove last day
  yearFit=rbind(yearFit,curfits)
}
par(mar=c(5,5,4,3))
meanFit=apply(yearFit,2,mean)
meanFitSum=NULL
for(i.day in 1:365){
  meanFitSum=c(meanFitSum,sum(rep(meanFit,2)[i.day:(i.day+duration-1)]))
}

x11(width=9,height=6)
if(plotExtra==TRUE){
  for(curgen in seq(2,length(years.index),length=10)){
    curgen=round(curgen)
    arheight=rep(max(meanFit)*1.1,N)
    emergeDay=pophistory[[curgen]]$emerge
    plot(meanFit,type='l',ylim=c(0,max(meanFit)*1.2))
    arrows(y0=jitter(arheight,factor=1.5),x0=emergeDay,x1=emergeDay+duration-1,length=.1)
    dev.print(pdf,paste("dailyfit-run",runName,"-gen",curgen,"-meanfit.pdf",sep=""))
    plot(meanFitSum,type='l',ylim=c(0,max(meanFitSum)*1.2),
         main=paste("Mean fitness gained, gen",curgen),
         ylab="Fitness gained",
         xlab="Julian date",
         cex.lab=1.3,
         cex.main=1.3)
    arheight=jitter(rep(max(meanFitSum)*1.05,N),factor=.8)
    arrows(y0=arheight+.05*max(meanFitSum),x0=emergeDay,y1=arheight,length=.1)
    dev.print(pdf,paste("dailyfitSum-run",runName,"-gen",curgen,"-meanfit.pdf",sep=""))
    #now calculate the fitSum for THIS YEAR ONLY
    FitSum=NULL
    for(i.day in 1:365){
      FitSum=c(FitSum,sum(rep(years.list[[years.index[[curgen]]]]$fit.daily,2)[i.day:(i.day+duration-1)]))
    }
    plot(FitSum,type='l',ylim=c(0,max(FitSum)*1.2),
         main=paste("Fitness gained this year, gen",curgen),
         ylab="Fitness gained",
         xlab="Julian date",
         cex.lab=1.3,
         cex.main=1.3)
    arheight=jitter(rep(max(FitSum)*1.05,N),factor=.8)
    arrows(y0=arheight+.05*max(FitSum),x0=emergeDay,y1=arheight,length=.1)
    dev.print(pdf,paste("dailyfitSum-run",runName,"-gen",curgen,"-actualfit.pdf",sep=""))
    #Now plot each of the coefs by emergence day.
    for(coefName in c("b.day","b.temp","b.precip")){
      plot(pophistory[[curgen]][,coefName],pophistory[[curgen]][,"emerge"],
           main=paste(coefName, "by emergence, gen", curgen),
           ylab="emergence",
           xlab=coefName,
           cex.lab=1.3,
           cex.main=1.3)
      dev.print(pdf,paste("Coef_x_emerge-",coefName,"-run",runName,"-gen",curgen,".pdf",sep=""))
    }
    coefList=c("b.day","b.temp","b.precip")
    for(i in 1:3){
      coef1=coefList[i]
      coef2=coefList[(i %% 3)+1]
      curpop=pophistory[[curgen]]
      plot(x=curpop[,coef1],y=curpop[,coef2],type='n',
           main=paste(coef1, "by", coef2, "gen", curgen),
           ylab=coef2,
           xlab=coef1,
           cex.lab=1.3,
           cex.main=1.3)
      points(x=(curpop[curpop[,"emerge"]>364,coef1]),y=(curpop[curpop[,"emerge"]>364,coef2]),pch=3,col='blue')
      points(x=(curpop[curpop[,"emerge"]<365,coef1]),y=(curpop[curpop[,"emerge"]<365,coef2]),pch=1,col='black')
      dev.print(pdf,paste("Coef_x_coef-",coef1,"x",coef2,"-run",runName,"-gen",curgen,".pdf",sep=""))
    }
  }
}
#Calculating changes in mean fitness through time
maxfit=maxActfit=meanfit=emerge.ideal=rep(0,length(years.index))
emerge=matrix(0,ncol=N,nrow=length(years.index))
for(curgen in 1:length(years.index)){
  meanfit[curgen]=mean(pophistory[[curgen]]$Wi)
  maxActfit[curgen]=max(pophistory[[curgen]]$Wi)
  cur.fitness=years.list[[years.index[curgen]]]$fit.daily
  cur.fitness.durated=rollapply(c(cur.fitness,rep(0,duration-1)),duration,by=1,sum)
  maxfit[curgen]=max(cur.fitness.durated)
  emerge[curgen,]=pophistory[[curgen]]$emerge
  emerge.ideal[curgen]=min(which(cur.fitness.durated==maxfit[curgen]))
}
#plot emergence times
maxCount=100 #maximum number of years to count
generations=1:length(years.index)
viewGens=generations
if(length(generations)>maxCount){
  viewGens=floor(seq(min(generations),max(generations),length.out=maxCount))
}
matplot(jitter(viewGens),emerge[viewGens,],type='p',pch=1,col='black',
     main=paste("Emergence days"),
     xlab="Generation",
     ylab="Emergence day",
     cex.lab=1.4,cex.main=1.4)
#Plot the "emerge before last day" indivs
points(viewGens,emerge.ideal[viewGens],col="red",pch=4,lwd=2)
dev.print(pdf,paste("emerge-run",runName,"-gen",curgen,".pdf",sep=""))

#plot mean fitness through time, showing max possible fitness
plot(maxfit,type='l',col='red',
     main=paste("Mean fitness through time for run",runName),
     xlab="generation",
     ylab="Raw mean fitness",
     sub="red is maximum possible",
     ylim=c(0,max(maxfit))
)
points(1:length(meanfit),meanfit,type='l')
dev.print(pdf,paste("meanfit-run",runName,"-gen",curgen,".pdf",sep=""))
#plot mean fitness through time, normalized by max fitness
plot(meanfit/maxfit,type='l',col='black',
     main=paste("Mean fitness / max possible",runName),
     ylab="normalized mean fitness",
     xlab="generation",
     sub="red is maximum possible",
     ylim=c(0,1)
)
abline(h=1,col='red')
dev.print(pdf,paste("meanfitNorm-run",runName,"-gen",curgen,".pdf",sep=""))
#plot max actual fitness through time
plot(maxfit,type='l',col='red',
     main=paste("Max achieved fitness through time for run",runName),
     xlab="generation",
     ylab="Raw max fitness",
     sub="red is maximum possible",
     ylim=c(0,max(maxfit))
)
points(1:length(maxActfit),maxActfit,type='l')
dev.print(pdf,paste("maxfit-run",runName,"-gen",curgen,".pdf",sep=""))

#Looking at coef changes through time
#  doing two things. exp.eff is the "expected" effect size by muliplying the coefficient of each individual by the mean environmental conditions for that generation
#  The act.eff is the actual effect size, found by multiplying the coefficient of each indiv by the environmental conditions of their day of emergence.
#    Those plots use crosses to represent individuals who didn't emerge until the final day, and circles for those that emerged on a normal day (ie their cue
exp.eff=meanTraitEff(years.index,years.list,pophistory,N)
act.eff=actTraitEff(years.index,years.list,pophistory,N)
act.vals=actTraitVals(pophistory,numYears,N)

x11(width=9,height=6)
par(mar=c(5,5,4,4))
if(plotExtra==TRUE){
  traitplot(indivs=act.vals,trait="b.day")
  dev.print(pdf,paste("coefVals-day-run",runName,".pdf",sep=""))
  traitplot(indivs=act.vals,trait="b.temp")
  dev.print(pdf,paste("coefVals-temp-run",runName,".pdf",sep=""))
  traitplot(indivs=act.vals,trait="b.precip")
  dev.print(pdf,paste("coefVals-precip-run",runName,".pdf",sep=""))

  traiteffplot(indivs=exp.eff,trait="b.day")
  dev.print(pdf,paste("coefEffects-day-expected-run",runName,".pdf",sep=""))
  traiteffplot(indivs=exp.eff,trait="b.temp")
  dev.print(pdf,paste("coefEffects-temp-expected-run",runName,".pdf",sep=""))
  traiteffplot(indivs=exp.eff,trait="b.precip")
  dev.print(pdf,paste("coefEffects-precip-expected-run",runName,".pdf",sep=""))

  emergePlot(indivs=act.eff,trait="b.day")
  dev.print(pdf,paste("coefEffects-day-actual-run",runName,".pdf",sep=""))
  emergePlot(indivs=act.eff,trait="b.temp")
  dev.print(pdf,paste("coefEffects-temp-actual-run",runName,".pdf",sep=""))
  emergePlot(indivs=act.eff,trait="b.precip")
  dev.print(pdf,paste("coefEffects-precip-actual-run",runName,".pdf",sep=""))
}
#Plot it all in one
x11()
par(mfrow=c(3,1))
traiteffplot(indivs=exp.eff,trait="b.day")
traiteffplot(indivs=exp.eff,trait="b.temp")
traiteffplot(indivs=exp.eff,trait="b.precip")
dev.print(pdf,paste("coefEffects-all-expected-run",runName,".pdf",sep=""))

traiteffplot(indivs=exp.eff[exp.eff[,"gen"]>burnIn,],trait="b.day")
traiteffplot(indivs=exp.eff[exp.eff[,"gen"]>burnIn,],trait="b.temp")
traiteffplot(indivs=exp.eff[exp.eff[,"gen"]>burnIn,],trait="b.precip")
dev.print(pdf,paste("coefEffects-all-expected-postburn-run",runName,".pdf",sep=""))

emergePlot(indivs=act.eff,trait="b.day")
emergePlot(indivs=act.eff,trait="b.temp")
emergePlot(indivs=act.eff,trait="b.precip")
dev.print(pdf,paste("coefEffects-all-actual-run",runName,".pdf",sep=""))

emergePlot(indivs=act.eff[act.eff[,"gen"]>burnIn,],trait="b.day")
emergePlot(indivs=act.eff[act.eff[,"gen"]>burnIn,],trait="b.temp")
emergePlot(indivs=act.eff[act.eff[,"gen"]>burnIn,],trait="b.precip")
dev.print(pdf,paste("coefEffects-all-actual-postburn-run",runName,".pdf",sep=""))

traitplot(indivs=act.vals,trait="b.day")
traitplot(indivs=act.vals,trait="b.temp")
traitplot(indivs=act.vals,trait="b.precip")
dev.print(pdf,paste("coefVals-all-run",runName,".pdf",sep=""))

traitplot(indivs=act.vals[act.vals[,"gen"]>burnIn,],trait="b.day")
traitplot(indivs=act.vals[act.vals[,"gen"]>burnIn,],trait="b.temp")
traitplot(indivs=act.vals[act.vals[,"gen"]>burnIn,],trait="b.precip")
dev.print(pdf,paste("coefVals-all-postburn-run",runName,".pdf",sep=""))

