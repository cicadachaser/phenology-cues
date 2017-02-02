#setting up year
create_yrs=function(meanYr,numYears,yearstd,daystd){
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
    fit.tot=c(fit.tot[-(1:lag)],rep(0,lag))#for lag of 1
    newYear=cbind(newYear,fit.daily=fit.daily,fit.tot=fit.tot)
    years.list[[count]]<-newYear
    count=count+1
  }
  return(years.list)
}


#################################################
#Calculating optima using objective function that encorporates fitness function

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


#####################################
# Explicitly solving for fitness surface, fast methods

opt_day<-function(years.list){
  #Function for explicitly calculating geometric fitness for emerging on each day
  #Fastest approach so far
  numYears=length(years.list)
  fitstore=matrix(0,nrow=numYears,ncol=365) #matrix to store fitness per day
  for(i.list in 1:numYears){
    fitstore[i.list,]=years.list[[i.list]][,"fit.tot"]
  }
  geofit=apply(log(fitstore),2,sum)
  return(c(geofit=max(geofit),b.day=which(geofit==max(geofit))-.1))
}

opt_temp=function(years.list){
  #Function for explicitly calculating geometric fitness for all possible temperature choices
  #Fastest approach so far
  numYears=length(years.list)
  ptstore=rep(0,365*numYears)
  for(i.list in 1:numYears){
    ptstore[((i.list-1)*365+1):(i.list*365)]=cummax(years.list[[i.list]][,"temp"])
  }
  ptstore=sort(unique(ptstore))
  testpts=ptstore[-length(ptstore)]+diff(ptstore)/2 #take midpoint between each ptstore
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
  geofit=apply(log(fitstore),2,sum)
  return(c(geofit=max(geofit),b.temp=testpts[geofit==max(geofit)][1]))
}



opt_cutemp<-function(years.list){
  #Function for explicitly calculating geometric fitness for all possible cumulative temperature choices
  #Fastest approach so far
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
    cmaxcutemp=cummax(years.list[[i.list]][,"cutemp"]) #cumulative maximum temp
    unik <- !duplicated(cmaxcutemp)  ## logical vector of unique values (where first appearance of number is "unique"
    mapseq=seq_along(cmaxcutemp)[unik]  ## indices of these uniques
    cmaxcutemp=cmaxcutemp[unik] #version of cmaxtemp
    fcur=approxfun(x=c(-10,cmaxcutemp,max(testpts)*1.1),
                   y=c(0,1:(length(cmaxcutemp)+1)))
    ind=ceiling(fcur(testpts)) #indices to map to appropriate days of year using mapseq
    days=c(1,mapseq,365)[ind+1]
    fit=years.list[[i.list]][days,"fit.tot"]
    fitstore[i.list,]=fit
  }
  geofit=apply(log(fitstore),2,sum)
  # plot(testpts,geofit,type='l')
  # return(geofit)
  return(c(geofit=max(geofit),b.cutemp=testpts[geofit==max(geofit)][1]))
}

#######################
# Fitness versions of explicit solvers (SLOWER)

fit_temp=function(years.list){
  #Function for explicitly calculating geometric fitness for all possible temperature choices
  #USES FITNESS FUNCTION - SLOWER THAN opt_temp
  numYears=length(years.list)
  ptstore=rep(0,365*numYears)
  for(i.list in 1:numYears){
    ptstore[((i.list-1)*365+1):(i.list*365)]=cummax(years.list[[i.list]][,"temp"])
  }
  ptstore=sort(unique(ptstore))
  testpts=ptstore[-length(ptstore)]+diff(ptstore)/2 #take midpoint between each ptstore
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
  geofit=apply(log(yrfit),2,sum)
  # return(geofit)
  return(c(geofit=max(geofit),b.temp=testpts[geofit==max(geofit)][1]))
}

fit_cutemp<-function(years.list){
  #Function for explicitly calculating geometric fitness for all possible temperature choices
  #USES FITNESS FUNCTION - SLOWER THAN opt_cutemp
  numYears=length(years.list)
  ptstore=rep(0,365*numYears)
  for(i.list in 1:numYears){
    ptstore[((i.list-1)*365+1):(i.list*365)]=years.list[[i.list]][,"cutemp"]
  }
  ptstore=sort(unique(ptstore))
  testpts=ptstore[-length(ptstore)]+diff(ptstore)/2 #take midpoint between each ptstore
  b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,length(testpts))
  # start with traits = 0
  indiv<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
  # create "empty" individual
  traits="cutemp"
  b.traits=sprintf("b.%s",traits)  # list of traits in b.day format
  indiv[,b.traits]=testpts #fill in the values for the traits we're actually using
  yrfit=matrix(0,nrow=length(years.list),ncol=length(testpts)) #initialize vector for storing yearly fitness
  for(i in 1:length(years.list)){ #iterate through all years, calculate fitness for each
    yrfit[i,]=fitness(year=years.list[[i]],newpop=indiv,duration=duration,traits=traits)$fit
  }
  geofit=apply(log(yrfit),2,sum)
  # return(geofit)
  return(c(geofit=max(geofit),b.temp=testpts[geofit==max(geofit)][1]))
}

fit_day<-function(years.list){
  #Function for explicitly calculating geometric fitness for emerging on each day
  #USES FITNESS FUNCTION - SLOWER THAN opt_day!
  testpts=(1:365)-.1
  numYears=length(years.list)
  b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,length(testpts))
  # start with traits = 0
  indiv<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
  # create "empty" individual
  traits="day"
  b.traits=sprintf("b.%s",traits)  # list of traits in b.day format
  indiv[,b.traits]=testpts #fill in the values for the traits we're actually using
  yrfit=matrix(0,nrow=length(years.list),ncol=length(testpts)) #initialize vector for storing yearly fitness
  for(i in 1:length(years.list)){ #iterate through all years, calculate fitness for each
    yrfit[i,]=fitness(year=years.list[[i]],newpop=indiv,duration=duration,traits=traits)$fit
  }
  geofit=apply(log(yrfit),2,sum)
  # return(geofit)
  return(c(geofit=max(geofit),b.temp=testpts[geofit==max(geofit)][1]))
}
