fit.parms=list(best.temp=40, shape.temp=10, #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
               best.other=10, shape.other=30, #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
               corr=.5) #correlation between the two terms. allows us to rotate the ellipse that we get in 2d space.
sigmat=with(fit.parms,matrix(c(shape.temp^2,shape.temp*shape.other*corr,shape.temp*shape.other*corr,shape.other^2),2,2,byrow=TRUE))
fit.parms$sigmat=sigmat

fit_fn<-function(years, #data frame containing temp and the "other.name"
                 other.name, #name of the "other" variable
                 fit.parms){ #the list of fitness parameters
  require(mvtnorm)
  #note: using with() to essentially attach fit.parms, but without using attach, which pushes them to the global environment (ewwww)
  with(fit.parms, dmvnorm(as.matrix(cbind(years$temp,years[,other.name])),mean=c(best.temp,best.other), sigma=sigmat))
}
