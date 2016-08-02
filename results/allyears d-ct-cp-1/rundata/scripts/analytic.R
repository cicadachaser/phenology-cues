#Script for optimizing traits for geometric fitness
#uses same parameter script files as simulation
rm(list=setdiff(ls(), c("set_wrkdir","runsname")))
##############################################################################
# The following lines are only important if you're running this script alone #
##############################################################################
# rm(list=ls()) #clear workspace
# runsname="parameterexample"
# pointcheck=10000 #number of points to evaluate initially
# fastnum=20 #number of points to test quickly
# slownum=10 #number of points to test slowly
#  set_wrkdir<-function(){
#   #function for setting working directory to the right place given the current computer/user
#   if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin" || Sys.getenv("USERNAME")=="Collin.work"){ #If it's collin
#     if(Sys.info()["nodename"]=="DESKTOP-D6QSU8F"){
#       setwd("G:\\Repos\\phenology-cues") #desktop
#     }else{
#       setwd("C:\\Repos\\phenology-cues") #desktop
#     }
#   }else{
#     if(Sys.getenv("COMPUTERNAME")=="ENT-YANG01"){
#       setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")#desktop
#     }else{
#       setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")} #laptop
#   }
# }

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
  return(-sum(log(yrfit+1/10^10)))
}
################################
# Optimizing geometric fitness #
################################
#Prepping:
#runsname="parameterexample" #Name of the parameter file, minus .R part
set_wrkdir()
require(zoo) #for use in fitness calculations
require(lhs) #for use in choosing initial points to check
source("scripts/windows_subs.R") #read in standard functions
source(paste("parameters/",runsname,".R",sep="")) #read in parameters
source(paste("fitcurve/",fitshape,".R",sep="")) #read in fitness curve function
source("scripts/rate_setup.R") #here using only for the cuesmax
years.list=NULL #initialize list

if(runType=="standard"){
  years.stuff=yeargen.davis(best.temp = best.temp,sd.temp = sd.temp,
                            best.precip = best.precip,sd.precip = sd.precip)
  years.list=years.stuff[[1]]
  years.index=years.stuff[[2]]
} else if(runType=="ithaca"){
  years.stuff=yeargen.ithaca(best.temp = best.temp,sd.temp = sd.temp,
                             best.precip = best.precip,sd.precip = sd.precip)
  years.list=years.stuff[[1]]
  years.index=years.stuff[[2]]
}
set_wrkdir()



#first, point-check
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
yrfit=matrix(0,nrow=N,ncol=length(years.index))
for(i in 1:length(years.index)){
  pop<-fitness(year=years.list[[years.index[i]]],newpop=newpop,duration=duration,traits=traits)
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
  temp=store.fast[[i]]=optim(par=topop[i,],fn=obj_fn,duration=duration,
                             yrs=years.list[years.index],traits=traits,
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

