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
fitfile="standgauss"
locName="davis"
traits=traitslist=c("day","temp","cutemp")
yearvarMax=2500 #max year variation to test
dayvarMax=25 #max day variation to test
numpts=10 #number of different levels of day and year variation to test (total of numpts^2 var combinations)
lag=1 #
duration=10
pointcheck=1000
numYears=100
fastnum=20
slownum=10
runnum="test4" #for identifying run number

#Creating combinations of year and day variations using matrix trick
yrvarmat=matrix(seq(0,yearvarMax,length=numpts),ncol=numpts,nrow=numpts,byrow = TRUE)
dayvarmat=matrix(seq(0,dayvarMax,length=numpts),ncol=numpts,nrow=numpts,byrow = FALSE)
yearvars=yrvarmat[1:(numpts^2)]
dayvars=dayvarmat[1:(numpts^2)]

require(zoo) #for use in fitness calculations
require(lhs) #for use in choosing initial points to check

set_wrkdir()
source(file = paste("fitcurve/",fitfile,".R",sep=""))

load(paste("data-years/",locName,"Dat.Rdata",sep = ""),envir= e <- new.env())
precip.means=e$daily.means$PRCP.means
temp.means=e$daily.means$TMAX.means
remove(list="e")

#smoothing across three repretitions of the year, using middle one (so that boundaries touch)
smooth.tempyrs=loess.smooth(rep(1:(365*3)),rep(temp.means[1:365],3),evaluation = 365*3,span=1/18)
# plot(smooth.tempyrs)
# points(rep(temp.means[1:365],3))
smooth.temp=smooth.tempyrs$y[366:(366+364)]

smooth.precipyrs=loess.smooth(rep(1:(365*3)),rep(precip.means[1:365],3),evaluation = 365*3,span=1/18)
# plot(smooth.precipyrs)
# points(rep(precip.means[1:365],3))
smooth.precip=smooth.precipyrs$y[366:(366+364)]

#CURRENTLY USING PRECIP = 0!
meanYr=as.data.frame(cbind(1:365,smooth.temp,rep(0,365)))
colnames(meanYr)=c("day","temp","precip")


source("scripts/var_exper/prepstuff.R")
mutdist=0 #to avoid error in the part of rate_setup.R that we're not using
source("scripts/rate_setup.R")
source("scripts/windows_subs.R")

totnum=length(yearvars)*length(traitslist)*slownum
overall.res=data.frame(dayvar=rep(-99,totnum),
                       yearvar=rep(-99,totnum),
                       trait=rep("init",totnum),
                       geofit=rep(-99,totnum),
                       traitval=rep(-99,totnum),
                       stringsAsFactors = FALSE)
resind=1 #results index for tracking where to put each iteration of results

for(i.var in 1:length(yearvars)){
  ###################################################
  # Produce a set of years of appropriate variances #
  ###################################################
  yearvar=yearvars[i.var]
  dayvar=dayvars[i.var]

  #generate years
  years.list=list()
  count=1
  for(i.year in 1:numYears){
    offset=round(rnorm(1,sd=yearvar))
    newtemp= meanYr[((1:365+offset[1]) %% 365)+1,2]
    newtemp=newtemp+rnorm(365,sd=sqrt(dayvar))
    newYear=as.data.frame(cbind(day=1:365,
                                temp=newtemp,
                                precip=meanYr$precip,
                                cuprecip=cumsum(meanYr$precip),
                                cutemp=cumsum(newtemp),
                                daysq=(1:365)^2,
                                tempsq=(newtemp)^2,
                                precipsq=(meanYr$precip)^2,
                                cutempsq=cumsum((newtemp)^2),
                                cuprecipsq=cumsum((meanYr$precip)^2)
    ))
    fit.daily=fit_fn(newYear)
    newYear=cbind(newYear,fit.daily=fit.daily)
    years.list[[count]]<-newYear
    count=count+1
  }

  #####################
  for(trait in traitslist){ #iterate through each trait
    N=pointcheck
    b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,N)
    lhc=randomLHS(N,length(trait))
    count=1
    for(i.trait in trait){
      if(start[[i.trait]][1]==0 & start[[i.trait]][2]==0){
        curvals=rep(0,N)
      }else{
        randnums=lhc[,count]*maxcues[[i.trait]]*length(trait)
        randnums[randnums==0]=1/(10^10)
        curvals=randnums
      }
      curname=paste("b.",i.trait,sep="")
      assign(curname,curvals)
      count=count+1
    }

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
    overall.res[resind:(resind+slownum-1),]=cbind(rep(dayvar,slownum),rep(yearvar,slownum),rep(paste(trait,collapse=", "),slownum),res.slow)
    resind=resind+slownum
    print(resind/totnum)
  }
}

set_wrkdir()
dir.create(paste("results/fig1/",runnum,sep=""))
save.image(file=paste("results/fig1/",runnum,"/fig1dat-version",runnum,".Rdata",sep=""))
source("scripts/var_exper/fig1_plot.R")
proc.time()
