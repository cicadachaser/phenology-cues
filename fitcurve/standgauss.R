fit.parms=list(best.temp=40, sd.temp=10, #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
               best.other=10, sd.other=30) #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
fit_fn<-function(years, #data frame containing temp and the "other.name"
                 other.name, #name of the "other" variable
                 fit.parms){ #the list of fitness parameters
  #note: using with() to essentially attach fit.parms, but without using attach, which pushes them to the global environment (ewwww)
  with(fit.parms, dnorm(years$temp,mean=best.temp,sd=sd.temp)*dnorm(years[,other.name],mean=best.other,sd=sd.other))
}
