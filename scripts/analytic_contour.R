#contour plotting
#For using with analytic script - requires parameter values from analytic.
require(graphics)
pointdense=100
plot.traits=c("day","temp")
fix.traits=c("precip")
fix.vals=cbind(rep(80,pointdense^2))
x=rep(seq(1/10^10,maxcues[[plot.traits[1]]]*length(traits),length=pointdense),pointdense)
ytemp=seq(1/10^10,maxcues[[plot.traits[2]]]*length(traits),length=pointdense)
y=x*0
for(i in 1:length(ytemp)){
  y[(1+(i-1)*pointdense):(i*pointdense)]=rep(ytemp[i],pointdense)
}

calc.traits=c(plot.traits,fix.traits)
calc.vals=cbind(x,y,fix.vals)

N=pointdense^2
b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,N)
count=1
for(i.trait in calc.traits){
  randnums=calc.vals[,count]
  randnums[randnums==0]=1/(10^10)
  curvals=randnums
  curname=paste("b.",i.trait,sep="")
  assign(curname,curvals)
  count=count+1
}
newpop<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
yrfit=matrix(0,nrow=N,ncol=length(years.index))
for(i in 1:length(years.index)){
  pop<-fitness(year=years.list[[years.index[i]]],newpop=newpop,duration=duration,traits=traits)
  yrfit[,i]=pop$fit
}
geofit=apply(yrfit,1,function(x){-sum(log(x))})
z=matrix(geofit,pointdense,pointdense,byrow=TRUE)

x11()
par(mar=c(7,5,5,3))
contour(x=seq(1/10^10,maxcues[[plot.traits[1]]]*length(traits),length=pointdense),
        y=seq(1/10^10,maxcues[[plot.traits[2]]]*length(traits),length=pointdense),
        z=z,
        xlab=plot.traits[1],
        ylab=plot.traits[2],
        sub=paste(fix.traits,"set to",fix.vals),
        main=paste("Contour of", runsname,"\n small is good"),
        cex.lab=1.6,cex.main=1.6,cex.sub=1.4)
#marking out infinities
inf.geo=which(is.infinite(geofit))
points(x[inf.geo],y[inf.geo],pch='.')