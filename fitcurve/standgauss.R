best.temp=45; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.precip=10; sd.precip=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
fit_fn<-function(years){
 dnorm(years$temp,mean=best.temp,sd=sd.temp)*dnorm(years$precip,mean=best.precip,sd=sd.precip)
}
