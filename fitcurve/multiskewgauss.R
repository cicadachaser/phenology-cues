#Fitness function based on multiplying two skew normal distribution
# Note that there is a complication - the mode (highest point) of the
# skew normal distribution doesn't have an exact answer. To continue using
# our current method of specification (ie "best.temp"), we solve for the
# "location" parameter of the skew normal that would give us a "best temp" or
# "best other" at the given value.

require(sn)
best.temp=30; shape.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.other=10; shape.other=40; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
skew.t=skew.o=-5 #skew values for temp and other skewed gaussians
level.corr=10 #correlation for the two covariates. Cannot exceed the smaller of shape.temp and shape.other!


#Converting best.temp and best.other to location parameters
obj_fun <- function(location,goal,maxrange,omega,alpha){
  #First, make a matrix that is two columns, each corresponding to temp and "other" values for each point to test
  #  Note that this should represent a grid we're testing.
  resolution = 1000 #number of points to test in each direction
  xraw.temp = seq(-(goal[1]), goal[1]*2, length=resolution) #all temperature values to test
  x.temp = as.vector(matrix(xraw.temp,resolution, resolution, byrow = FALSE)) #make a vector of repeated values
  xraw.other = seq(-(goal[2]), goal[2]*2, length=resolution) #all other values to test
  x.other = as.vector(matrix(xraw.other, resolution, resolution, byrow = TRUE)) #make a vector of repeated values
  xmat = cbind(x.temp,x.other)
  z = dmsn(xmat,xi=location,Omega=omega,alpha=alpha)
  return((xmat[which(z==max(z)),1]-goal[1])^2+(xmat[which(z==max(z)),2]-goal[2])^2)
}
omega=matrix(c(shape.temp,level.corr,level.corr,shape.other),2,2)
goal=c(best.temp,best.other)
out=optim(par=c(best.temp,best.other),
            fn=obj_fun,
            goal=goal,
            maxrange=80,
            omega=omega,
            alpha=c(skew.t,skew.o))
loc.t=out$par[1] #location value for temp that gives the target peak
loc.o=out$par[2] #location value for other that gives the target peak


#Turning our parameters into a list for feeding into the fit_fn
fit.parms=list(loc.t,loc.o,
               shape.temp,shape.other,
               skew.t,skew.o,omega)

# ## Use the lines below to visualize what the skew normal curves look like.
# x=seq(0,80,.1)
# plot(x,dsn(x,xi=loc.t,omega=shape.temp,alpha=skew.t)^2,type='l')
# abline(v=best.temp)
#
fit_fn<-function(years,other.name,fit.parms){
  #note: using with() to essentially attach fit.parms, but without using attach, which pushes them to the global environment (ewwww)
  with(fit.parms,
       dmsn(cbind(years$temp,years[,other.name]),xi=c(loc.t,loc.o),Omega=omega,alpha=c(skew.t,skew.o)))
}
