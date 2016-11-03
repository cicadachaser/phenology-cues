#Setting up for parallel processing
require(doSNOW); require(parallel); require(doParallel)
nClust=detectCores(all.tests=FALSE,logical=TRUE)
c1<-makeCluster(nClust-1)
registerDoParallel(c1)

#Getting basic information
source(paste("scripts/var_exper/runfiles/",runFile,sep = "")) #call the runFile which as all the variables defined
startTime=proc.time()
#Creating combinations of year and day variations using matrix trick
yrstdmat=matrix(seq(0,yearstdMax,length=numpts),ncol=numpts,nrow=numpts,byrow = TRUE)
daystdmat=matrix(seq(0,daystdMax,length=numpts),ncol=numpts,nrow=numpts,byrow = FALSE)
yearstds=yrstdmat[1:(numpts^2)]
daystds=daystdmat[1:(numpts^2)]

require(zoo) #for use in fitness calculations

#Reading in the appropriate fitness function
set_wrkdir()
source(file = paste("fitcurve/",fitfile,".R",sep=""))

#Reading in mean year, assigning to appropriate variables
load(paste("data-years/",locName,"Dat.Rdata",sep = ""),envir= e <- new.env())
precip.means=e$daily.means$PRCP.means
temp.means=e$daily.means$TMAX.means
remove(list="e")

#smoothing across three repetitions of the year, using middle one (so that boundaries touch)
smooth.tempyrs=loess.smooth(rep(1:(365*3)),rep(temp.means[1:365],3),evaluation = 365*3,span=1/18)
smooth.temp=smooth.tempyrs$y[366:(366+364)]

#Smoothing again, this time across precip
smooth.precipyrs=loess.smooth(rep(1:(365*3)),rep(precip.means[1:365],3),evaluation = 365*3,span=1/18)
smooth.precip=smooth.precipyrs$y[366:(366+364)]

#create mean year dataframe, using either the smoothed precip or all zeros for precip
if(actprecip==FALSE){
  meanYr=as.data.frame(cbind(1:365,smooth.temp,rep(0,365)))
}else{
  meanYr=as.data.frame(cbind(1:365,smooth.temp,smooth.precip))
}
colnames(meanYr)=c("day","temp","precip")

#Calculate the "base temp" for cumulative temp calculations, taking baseTempQth quantile
baseTemp=sort(meanYr$temp)[round(365/baseTempQ)]

#Pull in the appropriate functions
source("scripts/var_exper/prepstuff.R")
mutdist=0 #to avoid error in the part of rate_setup.R that we're not using
source("scripts/rate_setup.R")
source("scripts/windows_subs.R")

#Make a dataframe to store the results of each run
totnum=length(yearstds)*length(traitslist)*slownum #total number of runs

#Iterate through each each combination of stdev'
res=foreach(i.stdev = 1:length(yearstds)) %dopar% {
  #Re-load each library needed
  require(zoo) #for use in choosing initial points to check
  ###################################################
  # Produce a set of years of appropriate variances #
  ###################################################
  yearstd=yearstds[i.stdev]
  daystd=daystds[i.stdev]

  #generate years of climate
  years.list=list()
  count=1
  for(i.year in 1:numYears){
    offset=round(rnorm(1,sd=yearstd))
    newtemp= meanYr[((1:365+offset[1]) %% 365)+1,2]
    newtemp=newtemp+rnorm(365,sd=daystd)
    newYear=as.data.frame(cbind(day=1:365,
                                temp=newtemp,
                                precip=meanYr$precip,
                                cuprecip=cumsum(pmax(0,meanYr$precip)),
                                cutemp=cumsum(pmax(0,newtemp-baseTemp)),
                                daysq=(1:365)^2,
                                tempsq=(newtemp)^2,
                                precipsq=(meanYr$precip)^2,
                                cutempsq=cumsum(pmin(0,newtemp-baseTemp)^2),
                                cuprecipsq=cumsum(pmin(0,meanYr$precip)^2)
    ))
    fit.daily=fit_fn(newYear)
    fit.tot=c(rollapply(c(fit.daily,rep(0,duration-1)),duration,by=1,sum))
    newYear=cbind(newYear,fit.daily=fit.daily,fit.tot=fit.tot)
    years.list[[count]]<-newYear
    count=count+1
  }


  #####################
  #create initial starting points, check them
  resmat=NULL
  #FOR B.DAY
  dayres=opt_day(years.list)
  resmat=rbind(resmat,c(daystd,yearstd,paste(trait,collapse=", "),dayres))



  for(trait in traitslist){ #iterate through each trait
    N=pointcheck
    b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,N)
    randnums=runif(N)*maxcues[[i.trait]]*length(trait)
    randnums[randnums==0]=1/(10^10)
    curvals=randnums
    curname=paste("b.",i.trait,sep="")
    assign(curname,curvals)
    newpop<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
    yrfit=matrix(0,nrow=N,ncol=numYears)
    for(i in 1:numYears){
      pop<-fitness(year=years.list[[i]],newpop=newpop,duration=duration,traits=trait)
      yrfit[,i]=pop$fit
    }
    geofit=apply(yrfit,1,function(x){-sum(log(x))})
    ordpop=newpop[(order(geofit)),]
    b.traits=sprintf("b.%s",trait)
    topop=ordpop[1:fastnum,b.traits]

    #Now, fast-check the best points
    store.fast=list()
    res.fast=matrix(0,ncol=length(trait)+1,nrow=fastnum)
    colnames(res.fast)<-c("geofit",b.traits)
    if(length(b.traits)==1){
      for(i in 1:fastnum){
        temp=store.fast[[i]]=optim(par=ordpop[i,sprintf("b.%s",trait)],fn=obj_fn,duration=duration,
                                   yrs=years.list,traits=trait,method="Brent",lower=1e-5,upper=maxcues[[trait]],
                                   control=list(maxit=1000,abstol=1/10^4,reltol=1/10^4))
        res.fast[i,1]=temp$value
        res.fast[i,2:(length(trait)+1)]=temp$par
      }
    }
    res.fast=res.fast[order(res.fast[,"geofit"]),]

    #finally, slow convergence
    store.slow=list()
    res.slow=matrix(0,ncol=length(trait)+1,nrow=slownum)
    colnames(res.slow)<-c("geofit",b.traits)
    for(i in 1:slownum){
      temp=store.slow[[i]]=optim(par=res.fast[i,2:(1+length(trait))],fn=obj_fn,duration=duration,
                                 yrs=years.list,traits=trait,method="Brent",lower=1e-5,upper=maxcues[[trait]],
                                 control=list(maxit=50000,abstol=1/10^10))
      res.slow[i,1]=temp$value
      res.slow[i,2:(length(trait)+1)]=temp$par
    }
    #Save results, update counter, print progress
    resmat=rbind(resmat,cbind(rep(daystd,slownum),rep(yearstd,slownum),rep(paste(trait,collapse=", "),slownum),res.slow))
  }
  resmat
}
# df <- data.frame(matrix(unlist(res), ncol=5, byrow=F))
overall.res = do.call(rbind.data.frame, res)
names(overall.res)=c("daystd","yearstd","trait","geofit","traitval")

# Save results, make figures
set_wrkdir()
dir.create(paste("results/fig1/",runnum,sep=""))
save.image(file=paste("results/fig1/",runnum,"/fig1dat-version",runnum,".Rdata",sep=""))
source("scripts/var_exper/fig1_plot.R")
runTime=proc.time()-startTime;print(runTime)
