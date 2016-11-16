#Make sure to source prepstuff.R
load("G:/Repos/phenology-cues/results/fig1/opt3-2yr-comp1listlist/fig1dat-versionopt3-2yr-comp1listlist_new.Rdata")
source('G:/Repos/phenology-cues/scripts/windows_subs.R', encoding = 'UTF-8')
require(zoo)

duration=10;lag=1

years.list=yearlistlist[[5]]
#recalculate fit.tot with lag included
for(i.ylist in 1:length(years.list)){
    years.list[[i.ylist]]=years.list[[i.ylist]][,-which(names(years.list[[i.ylist]])=="fit.tot")]
   fit.tot=c(rollapply(c(years.list[[i.ylist]]$fit.daily,rep(0,duration-1)),duration,by=1,sum))
   fit.tot=c(fit.tot[-(1:lag)],rep(0,lag))#for lag of 1
   years.list[[i.ylist]]=cbind(years.list[[i.ylist]],fit.tot=fit.tot)

  }


numYears=length(years.list)
ptstore=rep(0,365*numYears)
for(i.list in 1:numYears){
  ptstore[((i.list-1)*365+1):(i.list*365)]=cummax(years.list[[i.list]][,"temp"])
}
ptstore=sort(unique(ptstore))
testpts=ptstore[-length(ptstore)]+diff(ptstore)/2 #take midpoint between each ptstore

#New approach:
fitstore=matrix(0,ncol=length(testpts),nrow=numYears)
for(i.list in 1:numYears){
  #fitting piecewise function of day~cutemp, but add in extreme endpoints to avoid extrapolation issues
  #First, though, need to remove duplicates, and create a way to map to the appropriate day
  cmaxtemp=cummax(years.list[[i.list]][,"temp"]) #cumulative maximum temp
  unik <- !duplicated(cmaxtemp)  ## logical vector of unique values (where first appearance of number is "unique"
  mapseq=seq_along(cmaxtemp)[unik]  ## indices of these uniques
  cmaxtemp=cmaxtemp[unik] #version of cmaxtemp
  fcur=approxfun(x=c(-10,cmaxtemp,max(testpts)*1.1),
                 y=c(0,1:(length(cmaxtemp)+1)))
  ind=ceiling(fcur(testpts)) #indices to map to appropriate days of year using mapseq
  days=c(1,mapseq,365)[ind+1]
  fit=years.list[[i.list]][days,"fit.tot"]
  fit[is.na(fit)]=0
  fitstore[i.list,]=fit
}
geofit.new=apply(log(fitstore),2,sum)


#Old approach:
b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,length(testpts))
# start with traits = 0
indiv<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
# create "empty" individual
traits="temp"
b.traits=sprintf("b.%s",traits)  # list of traits in b.day format
indiv[,b.traits]=testpts #fill in the values for the traits we're actually using
yrfit=matrix(0,nrow=length(years.list),ncol=length(testpts)) #initialize vector for storing yearly fitness
for(i in 1:length(years.list)){ #iterate through all years, calculate fitness for each
  yrfit[i,]=fitness(year=years.list[[i]],newpop=indiv,duration=duration,traits=traits)$fit
}
geofit.old=apply(log(yrfit),2,sum)
geofit.old-geofit.new

##################################

source('G:/Repos/phenology-cues/scripts/var_exper/prepstuff.R', encoding = 'UTF-8')
#Now confirming that the prepstuff version works:
official.temp=opt_temp(years.list)
opt.temp2=opt_temp2(years.list)
#Doesn't work!
official.temp-geofit.new
opt.temp2-geofit.new










#########
#Scratch
#########
