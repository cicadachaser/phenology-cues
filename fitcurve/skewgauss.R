#Fitness function based on multiplying two skew normal distribution
# Note that there is a complication - the mode (highest point) of the
# skew normal distribution doesn't have an exact answer. To continue using
# our current method of specification (ie "best.temp"), we solve for the
# "location" parameter of the skew normal that would give us a "best temp" or
# "best other" at the given value.

require(sn)
best.temp=40; shape.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.other=10; shape.other=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
skew.t=skew.o=-10 #skew values for temp and other skewed gaussians

#Converting best.temp and best.other to location parameters
obj_fun <- function(location,goal,maxrange,omega,alpha){
  x=seq(-(maxrange),maxrange,by=.1)
  z=dsn(x,xi=location,omega=omega,alpha=alpha)
  return((x[which(z==max(z))]-goal)^2)
}
goal=best.temp
loc.t=optimize(obj_fun,interval = c(goal/2,goal*2),
               goal=goal,maxrange=80,omega=shape.temp,alpha=skew.t)
loc.t=loc.t$minimum #Location parameter (xi) that gives target shape
goal=best.other
loc.o=optimize(obj_fun,interval = c(goal/2,goal*2),goal=goal,maxrange=80,omega=shape.other,alpha=skew.o)
loc.o=loc.o$minimum #Location parameter (xi) that gives target shape

#Turning our parameters into a list for feeding into the fit_fn
fit.parms=list(loc.t,loc.o,
               shape.temp,shape.other,
               skew.t,skew.o)

## Use the lines below to visualize what the skew normal curves look like.
# x=seq(0,80,.1)
# plot(x,dsn(x,xi=loc.t,omega=shape.temp,alpha=skew.t),type='l')
# abline(v=best.temp)

fit_fn<-function(years,other.name,fit.parms){
  #note: using with() to essentially attach fit.parms, but without using attach, which pushes them to the global environment (ewwww)
  with(fit.parms,
       dsn(years$temp,xi=loc.t,omega=shape.temp,alpha=skew.t)*
         dsn(years[,other.name],xi=loc.o,omega=shape.other,alpha=skew.o))
}
