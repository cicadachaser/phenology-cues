
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

#UNCERTAIN HOW THIS VERSION DIFFERS, BUT IT DOESN'T GIVE RIGHT ANSWER
# opt_temp=function(years.list){
#   numyears=length(years.list)
#   ptstore=rep(0,365*numYears)
#   for(i.list in 1:numYears){
#     ptstore[((i.list-1)*365+1):(i.list*365)]=cummax(years.list[[i.list]][,"temp"])
#   }
#   ptstore=sort(unique(ptstore))
#   testpts=ptstore[-length(ptstore)]+diff(ptstore)/2 #take midpoint between each ptstore
#   fitstore=matrix(0,ncol=length(testpts),nrow=numYears)
#   for(i.list in 1:numYears){
#     #fitting piecewise function of day~cutemp, but add in extreme endpoints to avoid extrapolation issues
#     cmaxtemp=cummax(years.list[[i.list]][,"temp"]) #cumulative maximum temp
#     unik <- !duplicated(cmaxtemp)  ## logical vector of unique values (where first appearance of number is "unique"
#     mapseq=seq_along(cmaxtemp)[unik]  ## indices of these uniques
#     cmaxtemp=cmaxtemp[unik] #version of cmaxtemp
#     fcur=approxfun(x=c(-10,cmaxtemp,max(testpts)*1.1),
#                    y=c(0,1:(length(cmaxtemp)+1)))
#     ind=ceiling(fcur(testpts)) #indices to map to appropriate days of year using mapseq
#     days=c(1,mapseq,365)[ind+1]
#     fit=years.list[[i.list]][days+lag,"fit.tot"]
#     fitstore[i.list,]=fit
#   }
#   geofit=apply(log(fitstore),2,sum)
#   # plot(testpts,geofit,type='l')
#   return(geofit)
#   # return(c(geofit=max(geofit),b.temp=testpts[geofit==max(geofit)][1]))
# }

opt_temp=function(years.list){
  #checking my functional code
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
