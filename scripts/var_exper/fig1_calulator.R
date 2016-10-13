fitfile="standgauss"
locName="davis"
lag=1
duration=10
pointcheck=1000
numYears=10

require(zoo) #for use in fitness calculations
require(lhs) #for use in choosing initial points to check

set_wrkdir()
source(file = paste("fitcurve/",fitfile,".R",sep=""))

e=load("data-years/",locName,"Dat.Rdata")
precip.means=daily.means$PRCP.means
temp.means=daily.means$TMAX.means
####RUN LOESS on these!



start<-list(#these are used to generate the starting values of individuals. Starting values will produce
  # individuals who, if they relied solely on a single cue, would emerge uniformly across a range of
  # cues from 1 to x times approximately the max value, where x is the number of traits used.
  # The max value is selected to be approximately the maximum found in the first five years of the davis
  # data set, and the min is set to approximately the minimum found. Note that we can't use zero.
  day=c(1,maxcues$day*length(traits)),
  temp=c(5,maxcues$temp*length(traits)),
  precip=c(.01,maxcues$precip*length(traits)),
  cutemp=c(5,maxcues$cutemp*length(traits)),
  cuprecip=c(.01,maxcues$cuprecip*length(traits)),
  daysq=c(1,maxcues$daysq*length(traits)),
  tempsq=c(5,maxcues$tempsq*length(traits)),
  precipsq=c(.01,maxcues$precipsq*length(traits)),
  cutempsq=c(5,maxcues$cutempsq*length(traits)),
  cuprecipsq=c(.01,maxcues$cuprecipsq*length(traits))
)




obj_fn<-function(x,duration,yrs,traits){
  #objective function for feeding to optim()
  #Parameters:
  # x is an individual (with traits)
  # duration is the duration an individual is emerged
  # yrs is the list of years to be used (ie climate regime)
  # traits is the list of traits being used
  #Returns:
  # the negative sum log of yearly fitness. This maps to geometic fitness,
  #   with small values being better (convenient for optim)
  b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=0
  # start with traits = 0
  indiv<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
  # create "empty" individual
  b.traits=sprintf("b.%s",traits)  # list of traits in b.day format
  indiv[b.traits]=x #fill in the values for the traits we're actually using
  yrfit=rep(0,length(yrs)) #initialize vector for storing yearly fitness
  for(i in 1:length(yrs)){ #iterate through all years, calculate fitness for each
    yrfit[i]=fitness(year=yrs[[i]],newpop=indiv,duration=duration,traits=traits)$fit
  }
  return(-sum(log(yrfit+.Machine$double.eps)))
}
################################################

#Iterate through param values
 yearVar=50
 dayVar=13

 ###################################################

#generate years
years.list=list()
count=1
for(i.year in 1:numYears){
  offset=round(rnorm(1,sd=yearVar))
  newtemp= meanYr[((1:365+offset[1]) %% 365)+1,2]
  newtemp=newtemp+rnorm(365,sd=sqrt(dayVar))
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
traits="day"
N=pointcheck
b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,N)
lhc=randomLHS(N,length(traits))
count=1
for(i.trait in traits){
  if(start[[i.trait]][1]==0 & start[[i.trait]][2]==0){
    curvals=rep(0,N)
  }else{
    randnums=lhc[,count]*maxcues[[i.trait]]*length(traits)
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
  pop<-fitness(year=years.list[[i]],newpop=newpop,duration=duration,traits=traits)
  yrfit[,i]=pop$fit
}
geofit=apply(yrfit,1,function(x){-sum(log(x))})
ordpop=newpop[(order(geofit)),]
b.traits=sprintf("b.%s",traits)
topop=ordpop[1:fastnum,b.traits]

#Now, fast-check the best points
store.fast=list()
res.fast=matrix(0,ncol=length(traits)+1,nrow=fastnum)
colnames(res.fast)<-c("geofit",b.traits)
for(i in 1:fastnum){
  print(i)
  temp=store.fast[[i]]=optim(par=ordpop[i,sprintf("b.%s",traits)],fn=obj_fn,duration=duration,
                             ###################I AM HERE
                             ###################TRYING TO MAKE OPTIM WORK RIGHT
                             ################### MAY NEED TO FIDDLE WITH LOWER AND UPPER BOUNDS
                             yrs=years.list,traits=traits,method="Brent",lower=,upper=,
                             control=list(maxit=1000,abstol=1/10^4,reltol=1/10^4))
  res.fast[i,1]=temp$value
  res.fast[i,2:(length(traits)+1)]=temp$par
}
#NOTE: we are MINIMIZING the objective function, which is the NEGATVE of the log geometric fitness
res.fast=res.fast[order(res.fast[,"geofit"]),]


#finally, slow convergence
store.slow=list()
res.slow=matrix(0,ncol=length(traits)+1,nrow=slownum)
colnames(res.slow)<-c("geofit",b.traits)
for(i in 1:slownum){
  print(i)
  temp=store.slow[[i]]=optim(par=res.fast[i,2:(1+length(traits))],fn=obj_fn,duration=duration,
                             yrs=years.list[years.index],traits=traits,
                             control=list(maxit=50000,abstol=1/10^10))
  res.slow[i,1]=temp$value
  res.slow[i,2:(length(traits)+1)]=temp$par
}
set_wrkdir()
save(list=c("store.fast","res.fast","store.slow","res.slow"),file=paste("results/analytic/",runsname,"_analytic.RData",sep=""))

