#For comparing two runs
require(stats)
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
load(paste("results/",runsnames[1],"/",runsnames[1],"_summary.RData",sep=""))
assign(x="sims1",value=list(means=store.mean,max=store.max,names=store.names,coEff=store.coEff,finalpops=finalpops))
load(paste("results/",runsnames[2],"/",runsnames[2],"_summary.RData",sep=""))
assign(x="sims2",value=list(means=store.mean,max=store.max,names=store.names,coEff=store.coEff,finalpops=finalpops))
max1=apply(sims1[["max"]],2,max)
max2=apply(sims2[["max"]],2,max)



setwd("results")
resultsdir=paste("compare-",names[1],"-vs-",names[2],sep="")
unlink(resultsdir,recursive = TRUE)
dir.create(resultsdir,showWarnings = FALSE)
setwd(resultsdir)


#plotting scaled means of all runs
col.list=c('red','blue')
dev.new(width=9,height=6)
matplot(cbind(t(sims1[["means"]][,(length(sims1$max)-viewLength):length(sims1$max)] /
                  sims1[["max"]][,(length(sims1$max)-viewLength):length(sims1$max)]),
              t(sims2[["means"]][,(length(sims2$max)-viewLength):length(sims2$max)] /
                  sims2[["max"]][,(length(sims2$max)-viewLength):length(sims2$max)])),
        type='l',col=c(rep("chocolate",numsims),rep("cornflowerblue",numsims)),
        main=paste("Scaled fitness through time for all runs",names[1],names[2]),
        xlab="generation",
        ylab="Raw mean fitness",
        cex.lab=1.5,
        cex.main=1.8
)
legend(x="bottomright",legend=c(sprintf("max %s",names)),
       fill=c("chocolate","cornflowerblue"),cex=2)
abline(h=1)
dev.print(pdf,paste("compare-allruns-scaled.pdf",sep=""))


#Plotting mean of means of runs
meanmeans.1=apply(sims1[["means"]],2,mean)
meanmax.1=apply(sims1[["max"]],2,mean)
meanmeans.2=apply(sims2[["means"]],2,mean)
meanmax.2=apply(sims2[["max"]],2,mean)

matplot(cbind(meanmax.1[(length(sims1$max)-viewLength):length(sims1$max)],
              meanmax.2[(length(sims2$max)-viewLength):length(sims2$max)]),
        type='l',col=c(rep('red',numsims),rep('blue',numsims)),
        main=paste("Mean fitness through time for mean of runs",names[1],names[2]),
        xlab="generation",
        ylab="Raw mean fitness",
        ylim=c(0,max(c(meanmax.1,meanmax.2)))
)
matpoints(cbind(meanmeans.1[(length(sims1$max)-viewLength):length(sims1$max)],
                meanmeans.2[(length(sims2$max)-viewLength):length(sims2$max)]),
          type='l',col=c("chocolate","cornflowerblue"))
legend(x="bottomright",legend=c(sprintf("max possible %s",names),sprintf("mean %s",names)),
       fill=c(col.list,"chocolate","cornflowerblue"),cex=2)
dev.print(pdf,paste("compare-means.pdf",sep=""))

# #Difference of means through time
# plot((meanmeans.1-meanmeans.2),type='l',col="red",
#      main=paste("Fitness of ",names[1],"minus",names[2]),
#      xlab="generation",
#      ylab="Raw mean fitness"
# )
# abline(h=0)
# dev.print(pdf,paste("difference-of-means.pdf",sep=""))

#histogram of difference of means
hist((sims1[["means"]][,(length(sims1$max)-viewLength):length(sims1$max)]-
        sims2[["means"]][,(length(sims2$max)-viewLength):length(sims2$max)]),
     breaks=30,
     main="histogram of differences",
     xlab=paste("Fitness of ",names[1],"minus",names[2]),
     sub="green is mean")
abline(v=0,col='red',lwd=2)
abline(v=mean(sims1[["means"]]-sims2[["means"]]),lwd=2,col='green')
dev.print(pdf,paste("hist-difference-of-means.pdf",sep=""))

#Scaled mean of means fitness through time
matplot(cbind(meanmeans.1[(length(sims1$max)-viewLength):length(sims1$max)] /
                meanmax.1[(length(sims1$max)-viewLength):length(sims1$max)],
              meanmeans.2[(length(sims2$max)-viewLength):length(sims2$max)] /
                meanmax.2[(length(sims2$max)-viewLength):length(sims2$max)]),
        type='l',col=c("chocolate","cornflowerblue"),
        main=paste("Scaled fitness through time for mean of runs",names[1],names[2]),
        xlab="generation",
        ylab="Raw mean fitness",
        ylim=c(0,1)
)
#Add trend for first param set
xvals=1:viewLength
model=loess((meanmeans.1[(length(sims1$max)-viewLength):length(sims1$max)] /
               meanmax.1[(length(sims1$max)-viewLength):length(sims1$max)])~xvals)
pred=predict(model,xvals)
points(pred,type='l',lwd=2,col="chocolate",lty=2)
#Add trend for second param set
xvals=1:viewLength
model=loess((meanmeans.2[(length(sims2$max)-viewLength):length(sims2$max)] /
               meanmax.2[(length(sims2$max)-viewLength):length(sims2$max)])~xvals)
pred=predict(model,xvals)
points(pred,type='l',lwd=2,col="cornflowerblue",lty=2)

abline(h=1,col='red')
legend(x="bottomright",legend=names,fill=c("chocolate","cornflowerblue"),cex=2)
dev.print(pdf,paste("compare-means-scaled.pdf",sep=""))

#barplot(ish) for the fitness values towards the end of the run
latefit1=apply(sims1[["means"]][,(length(sims1$max)-viewLength):length(sims1$max)],1,sum)
latefit2=apply(sims2[["means"]][,(length(sims1$max)-viewLength):length(sims1$max)],1,sum)
plot(jitter(c(rep(1,numsims),rep(2,numsims)),factor=.1),c(latefit1,latefit2),
     xlim=c(.5,2.5),ylim=c(0,max(c(latefit1,latefit2))),
     xaxt='n',
    main=paste("Comparing sum fitness over final",viewLength, "years"),
    xlab="",
    ylab="Sum fitness"
)
axis(1,at=c(1,2),labels = names)
dev.print(pdf,paste("latefitness.pdf",sep=""))
set_wrkdir()
dev.off()
