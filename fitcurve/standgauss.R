best.temp=40; sd.temp=10; #The optimal temp and the sd for the temp-by-fitness curve (which is gaussian)
best.other=10; sd.other=30; #The optimal precip and the sd for the precip-by-fitness curve (which is gaussian)
other.name="precip"
fit_fn<-function(years,other.name){
 dnorm(years$temp,mean=best.temp,sd=sd.temp)*dnorm(years[,other.name],mean=best.other,sd=sd.other)
}
