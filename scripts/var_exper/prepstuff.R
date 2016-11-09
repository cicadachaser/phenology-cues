
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
  return(c(geofit=max(geofit),b.day=which(geofit==max(geofit))-.1))
}

opt_temp=function(years.list){
  numyears=length(years.list)
  ptstore=rep(0,365*numYears)
  for(i.list in 1:numYears){
    ptstore[((i.list-1)*365+1):(i.list*365)]=cummax(years.list[[i.list]][,"temp"])
  }
  ptstore=sort(unique(ptstore))
  testpts=ptstore[-length(ptstore)]+diff(ptstore)/2 #take midpoint between each ptstore
  fitstore=matrix(0,ncol=length(testpts),nrow=numYears)
  for(i.list in 1:numYears){
    #fitting piecewise function of day~cutemp, but add in extreme endpoints to avoid extrapolation issues
    fcur=approxfun(x=c(-10,cummax(years.list[[i.list]][,"temp"]),max(testpts)*1.1),
                   y=c(1,years.list[[i.list]][,"day"],365))
    days=ceiling(fcur(testpts)) #mapping test points to emergence days
    fit=years.list[[i.list]][days,"fit.tot"]
    fitstore[i.list,]=fit
  }
  geofit=apply(log(fitstore),2,sum)
  # plot(testpts,geofit,type='l')
  return(c(geofit=max(geofit),b.temp=testpts[geofit==max(geofit)][1]))
}

opt_cutemp<-function(years.list){
  numYears=length(years.list)
  ptstore=rep(0,365*numYears)
  for(i.list in 1:numYears){
    ptstore[((i.list-1)*365+1):(i.list*365)]=years.list[[i.list]][,"cutemp"]
  }
  ptstore=sort(unique(ptstore))
  testpts=ptstore[-length(ptstore)]+diff(ptstore)/2 #take midpoint between each ptstore
  fitstore=matrix(0,ncol=length(testpts),nrow=numYears)
  for(i.list in 1:numYears){
    #fitting piecewise function of day~cutemp, but add in extreme endpoints to avoid extrapolation issues
    fcur=approxfun(x=c(-10,years.list[[i.list]][,"cutemp"],max(testpts)*1.1),
                   y=c(1,years.list[[i.list]][,"day"],365))
    days=ceiling(fcur(testpts)) #mapping test points to emergence days
    fit=years.list[[i.list]][days,"fit.tot"]
    fitstore[i.list,]=fit
  }
  geofit=apply(log(fitstore),2,sum)
  # plot(testpts,geofit,type='l')
  return(c(geofit=max(geofit),b.cutemp=testpts[geofit==max(geofit)][1]))
}
