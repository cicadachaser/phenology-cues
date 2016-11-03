
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


opt_day<-function(years.list){
  numYears=length(years.list)
  fitstore=matrix(0,nrow=numYears,ncol=365) #matrix to store fitness per day
  for(i.list in 1:numYears){
    fitstore[i.list,]=years.list[[i.list]][,"fit.tot"]
  }
  geofit=apply(log(fitstore),2,sum)
  return(c(geofit=max(geofit),b.day=which(geofit==max(geofit))))
}
