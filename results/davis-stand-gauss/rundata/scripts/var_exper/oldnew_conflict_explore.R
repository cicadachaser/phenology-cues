#seeing if the new method gives us the same fitness that the old method does
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
setwd("results/fig1/compare3-2yr")
load("fig1dat-versioncompare3-2yr.Rdata",envir=brute<-new.env())
set_wrkdir()
source("scripts/var_exper/runfiles/compare3-2yr.R") #call the runFile which as all the variables defined
set_wrkdir()
source("scripts/var_exper/runfiles/compare1-10yr.R") #call the runFile which as all the variables defined
mutdist=0 #to avoid error in the part of rate_setup.R that we're not using
source("scripts/rate_setup.R")
source("scripts/windows_subs.R")


stds=unique(cbind(daystd=varcompare$daystd,yearstd=varcompare$yearstd))
yearstds=stds[,"yearstd"]
daystds=stds[,"daystd"]
require(zoo)

oldnewfit=rep(0,dim(varcompare)[1])
for(i.test in 1:dim(varcompare)[1]){
  curline=varcompare[i.test,]
  listlistind=which(yearstds==curline$yearstd &daystds==curline$daystd)
  years.list=brute$yearlistlist[[listlistind]]
  b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=0
  # start with traits = 0
  indiv<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
  # create "empty" individual
  b.traits=sprintf("b.%s",curline$trait)  # list of traits in b.day format
  indiv[b.traits]=curline$traitval.new #fill in the values for the traits we're actually using
  yrfit=rep(0,length(years.list)) #initialize vector for storing yearly fitness
  for(i in 1:length(years.list)){ #iterate through all years, calculate fitness for each
    yrfit[i]=fitness(year=years.list[[i]],newpop=indiv,duration=duration,traits=traits)$fit
  }
  oldnewfit[i.test]=sum(log(yrfit+.Machine$double.eps))
}
fit.tot=c(rollapply(c(years.list[[1]]$fit.daily,rep(0,duration-1)),duration,by=1,sum))
for(i in 2:length(years.list))
  fit.tot=rbind(fit.tot,rollapply(c(years.list[[i]]$fit.daily,rep(0,duration-1)),duration,by=1,sum))
max(apply(log(fit.tot),2,sum))
oldnewfit-varcompare$geofit.old

#testing to see if fitness function behaves as expected
b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=0
# start with traits = 0
indiv<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
indiv$b.day=99.9
i=1
res=fitness(year=years.list[[i]],newpop=indiv,duration=duration,traits=traits)
fit.tot=c(rollapply(c(years.list[[1]]$fit.daily,rep(0,duration-1)),duration,by=1,sum))

yrfit()
