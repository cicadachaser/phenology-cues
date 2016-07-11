#For comparing two runs

set_wrkdir()
load(paste("results/",names[1],"/",names[1],"_summary.RData",sep=""))
assign(x="sims1",value=list(means=store.mean,names=store.names,coEff=store.coEff,finalpops))
load(paste("results/",names[2],"/",names[2],"_summary.RData",sep=""))
assign(x="sims2",value=list(means=store.mean,names=store.names,coEff=store.coEff,finalpops))



###########Everything below this needs to be updated for the new modular approach###########

setwd("results")
resultsdir=paste("compare",names)
dir.create(resultsdir,showWarnings = FALSE)
setwd(resultsdir)

col.list=c('red','blue')

matplot(((numYears-viewLength+1):numYears),rbind(sim1[["store.mean"]],sim2[["store.mean"]]),type='l',col=c(rep('red',numsims),rep('blue',numsims)),
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
#}