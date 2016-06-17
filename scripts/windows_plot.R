require(vegan)
require(scatterplot3d)
setwd("results")
resultsdir=sprintf("resRun%s",runName)
dir.create(resultsdir,showWarnings = FALSE)
setwd(resultsdir)

#Make a matrix of daily fitnesses for all years, cutting out the last day of leap years
yearFit=matrix(0,nrow=length(years.index),ncol=365)
for(i in years.index){
  curfits=years.list[[i]]$fit.daily
  if(length(curfits)==366){curfits=curfits[-366]} #to handle leap years, remove last day
  yearFit[i,]=curfits;
}
#Calculate the mean fitness accrued each day
par(mar=c(5,5,4,3))
meanFit=apply(yearFit,2,mean)
meanFitSum=NULL
#calculate the mean fitness for emerging on day x (for all days) [this is using mean fitness accrued per day]
for(i.day in 1:365){
  meanFitSum=c(meanFitSum,sum(rep(meanFit,2)[i.day:(i.day+duration-1)]))
}

x11(width=9,height=6)
if(plotExtra==TRUE){
  for(curgen in c(1:20,seq(21,length(years.index),length=10))){
    curgen=round(curgen)
    #arheight=rep(max(meanFit)*1.1,N) #upper bound used for making plots look good
    emergeDay=pophistory[[curgen]]$emerge
    #     plot(meanFit,type='l',ylim=c(0,max(meanFit)*1.2))
    #     arrows(y0=jitter(arheight,factor=1.5),x0=emergeDay,x1=emergeDay+duration-1,length=.1)
    #     dev.print(pdf,paste("dailyfit-run",runName,"-gen",curgen,"-meanfit.pdf",sep=""))
    #     plot(meanFitSum,type='l',ylim=c(0,max(meanFitSum)*1.2),
    #          main=paste("Mean fitness gained, gen",curgen),
    #          ylab="Fitness gained",
    #          xlab="Julian date",
    #          cex.lab=1.3,
    #          cex.main=1.3)
    #     arheight=jitter(rep(max(meanFitSum)*1.05,N),factor=.8)
    #     arrows(y0=arheight+.05*max(meanFitSum),x0=emergeDay,y1=arheight,length=.1)
    #     dev.print(pdf,paste("dailyfitSum-run",runName,"-gen",curgen,"-meanfit.pdf",sep=""))
    #now calculate the fitSum for THIS YEAR ONLY
    FitSum=NULL
    for(i.day in 1:365){
      FitSum=c(FitSum,sum(c(years.list[[years.index[[curgen]]]]$fit.daily,rep(0,duration))[i.day:(i.day+duration-1)]))
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
    #     for(coefName in c("b.day","b.temp","b.precip")){
    #       plot(pophistory[[curgen]][,coefName],pophistory[[curgen]][,"emerge"],
    #            main=paste(coefName, "by emergence, gen", curgen),
    #            ylab="emergence",
    #            xlab=coefName,
    #            cex.lab=1.3,
    #            cex.main=1.3)
    #       dev.print(pdf,paste("Coef_x_emerge-",coefName,"-run",runName,"-gen",curgen,".pdf",sep=""))
    #     }
    #     coefList=c("b.day","b.temp","b.precip")
    #     for(i in 1:3){
    #       coef1=coefList[i]
    #       coef2=coefList[(i %% 3)+1]
    #       curpop=pophistory[[curgen]]
    #       plot(x=curpop[,coef1],y=curpop[,coef2],type='n',
    #            main=paste(coef1, "by", coef2, "gen", curgen),
    #            ylab=coef2,
    #            xlab=coef1,
    #            cex.lab=1.3,
    #            cex.main=1.3)
    #       points(x=(curpop[curpop[,"emerge"]>364,coef1]),y=(curpop[curpop[,"emerge"]>364,coef2]),pch=3,col='blue')
    #       points(x=(curpop[curpop[,"emerge"]<365,coef1]),y=(curpop[curpop[,"emerge"]<365,coef2]),pch=1,col='black')
    #       dev.print(pdf,paste("Coef_x_coef-",coef1,"x",coef2,"-run",runName,"-gen",curgen,".pdf",sep=""))
    #     }
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
#Plot when organisms emerge
matplot(jitter(viewGens),emerge[viewGens,],type='p',pch=1,col='black',
        main=paste("Emergence days"),
        xlab="Generation",
        ylab="Emergence day",
        cex.lab=1.4,cex.main=1.4)
#add `optimal emergence day' - note this is a vast oversimplification
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
#plot max potential fitness and max actual fitness through time
#  ie the fittest individual of each generatoin
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
#  The act.eff is the actual effect size, found by multiplying the coefficient of each indiv by the environmental conditions of their day of emergence.
#    Those plots use crosses to represent individuals who didn't emerge until the final day, and circles for those that emerged on a normal day (ie their cue
act.eff=actTraitEff(years.index,years.list,pophistory,N,traits)
act.vals=actTraitVals(pophistory,numYears,N)
traitslist=sprintf("b.%s",traits)

#x11(width=9,height=6)
par(mar=c(5,5,4,4))
#Plot it all in one
par(mfrow=c(3,1))
ind=1
while(ind<=length(traitslist)){
  plotlist=ind:min(ind+2,length(traitslist))
  for(cur.trait in plotlist){
    emergePlot(indivs=act.eff,trait=traitslist[cur.trait])
  }
  dev.print(pdf,paste("coefEffects-",paste(traitslist[plotlist],collapse='-'),"-actual-run",runName,".pdf",sep=""))
  for(cur.trait in plotlist){
    emergePlotYlim(indivs=act.eff,trait=traitslist[cur.trait],ylim=c(-50,150))
  }
  dev.print(pdf,paste("coefEffects-",paste(traitslist[plotlist],collapse='-'),"-actual-run-ylim",runName,".pdf",sep=""))
  for(cur.trait in plotlist){
  emergePlot(indivs=act.eff[act.eff[,"gen"]>burnIn,],trait=traitslist[cur.trait])
  }
  dev.print(pdf,paste("coefEffects-",paste(traitslist[plotlist],collapse='-'),"-actual-postburn-run",runName,".pdf",sep=""))
  for(cur.trait in plotlist){
  traitplot(indivs=act.vals,trait=traitslist[cur.trait])
  }
  dev.print(pdf,paste("coefVals-",paste(traitslist[plotlist],collapse='-'),"-run",runName,".pdf",sep=""))
  for(cur.trait in plotlist){
  traitplot(indivs=act.vals[act.vals[,"gen"]>burnIn,],trait=traitslist[cur.trait])
    }
  dev.print(pdf,paste("coefVals-",paste(traitslist[plotlist],collapse='-'),"-postburn-run",runName,".pdf",sep=""))
  ind=max(plotlist)+1
}
#Plot phenotypes through time

par(mfrow=c(1,1))
pop.temp=do.call(rbind.data.frame,pophistory[c(seq(2,length(years.index),by=100))])
traitmins=NULL
for(i.trait in traitslist){traitmins=c(traitmins,min(pop.temp[,i.trait]))}
traitmaxs=NULL
for(i.trait in traitslist){traitmaxs=c(traitmaxs,max(pop.temp[,i.trait]))}
#x11(width=9,height=6)
if(plotPheno==TRUE){
  traitslist=sprintf("b.%s",traits)
  for(curgen in c(1,seq(2,length(years.index),by=100))){
    curgen=round(curgen)
    cur.pop=pophistory[[curgen]]
    if(length(traitslist)<4){ #can just do 3d plot
      scatterplot3d(jitter(as.matrix(cur.pop[,traitslist])),type='h',
                    xlim=c(traitmins[1],traitmaxs[1]),
                    ylim=c(traitmins[2],traitmaxs[2]),
                    zlim=c(traitmins[3],traitmaxs[3]))
      dev.print(pdf,sprintf("pheno3d-gen%06d-run%s.pdf",curgen,runName))
    }else{
      cur.ndms=metaMDS(cur.pop[,traitslist],k=2,trymax=100,autotransform = TRUE)
      plot(cur.ndms,main=paste("NDMS for generation", curgen))
      dev.print(pdf,sprintf("pheno3d-gen%06d-run%s.pdf",curgen,runName))
    }
  }
}
