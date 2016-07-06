#Does all the setup of population, runs simulation

#######################
# initializing population
#######################
##intialize a population of N individuals
# Their min and max values are determined by the start$ parameters
#Start by setting them all equal to zero. Fill in with the working traits from the traits variable
b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,N)
for(i.trait in traits){
  if(start[[i.trait]][1]==0 & start[[i.trait]][2]==0){
    curvals=rep(0,N)
  }else{
    randnums=runif(n=N,min=start[[i.trait]][1],max=start[[i.trait]][2])
    randnums[randnums==0]=1/(10^10)
    curvals=randnums
  }
  curname=paste("b.",i.trait,sep="")
  assign(curname,curvals)
}
newpop<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
pop<-selection(newpop,duration,year=years.list[[years.index[1]]],N,traits=traits)
###########################
## Running the Simulation #
###########################
pophistory=runSim(startpop=pop,years.list=years.list,
                  years.ind=years.index,N=N,duration=duration,
                  sds=sds,mutrate=mutrate,generations=numYears,traits=traits)
#Note: we've already used year 1 in initiating the pop
#Based on the value of "runType", generate the appropriate type of data.