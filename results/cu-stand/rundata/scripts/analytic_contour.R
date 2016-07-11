#contour plotting
#For using with analytic script - requires results and parameter values from analytic.
require(graphics)
# The below is needed only if you don't run the parameter file
# pointdense=100
# num.plots=10
# plot.traits=c("day","temp") #traits to vary
fix.traits=c("precip") #traits to hold constant
fix.b=sprintf("b.%s",fix.traits)
figpath=paste("results/analytic/",runsname,sep="")
dir.create(figpath,showWarnings = FALSE)
if(length(fix.traits)>1){ #if there are more than one fixed traits
  fix.max=apply(res.slow[,fix.b],2,max)
  fix.min=apply(res.slow[,fix.b],2,min)
  fix.valmat=apply(rbind(fix.min,fix.max),2,function(x){seq(x[1],x[2],length=num.plots)})
}else{
  fix.valmat=t(t(seq(min(res.slow[,fix.b]),max(res.slow[,fix.b]),length=num.plots)))
}

dev.new(width=9,height=6)
par(mar=c(7,5,5,3))
for(i.vals in 1:num.plots){
  fix.traits.values=fix.valmat[i.vals,] #values to be held constant, in same order as fix.traits (if there are multiple traits)
  #####################################
  fix.vals=matrix(fix.traits.values,nrow=pointdense^2,ncol=length(fix.traits.values),byrow=TRUE)
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
  contour(x=seq(1/10^10,maxcues[[plot.traits[1]]]*length(traits),length=pointdense),
          y=seq(1/10^10,maxcues[[plot.traits[2]]]*length(traits),length=pointdense),
          z=z,
          xlab=plot.traits[1],
          ylab=plot.traits[2],
          sub=paste(fix.traits,"set to",round(fix.traits.values)),
          main=paste("Contour of", runsname,"\n small is good"),
          cex.lab=1.6,cex.main=1.6,cex.sub=1.4)
  #marking out infinities
  inf.geo=which(is.infinite(geofit))
  points(x[inf.geo],y[inf.geo],pch='.')
  b.traits=sprintf("b.%s",plot.traits)
  points(res.slow[,b.traits[1]],res.slow[,b.traits[2]],col='red')
  #points(res.slow[,b.traits[1]][1],res.slow[,b.traits[2]][1],col='blue')
  dev.print(pdf,paste(figpath,sprintf("/fix-%s-at-%04.0f.pdf",fix.traits,fix.traits.values),sep=""))
}
dev.off()