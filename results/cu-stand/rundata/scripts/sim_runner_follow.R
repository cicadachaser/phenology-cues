
#Does all the setup of population, runs simulation
#######################
# initializing population
#######################
#This is the "follow" version, which reads in the old pophistory and uses that to initialize populations
name.old=paste("resRun",runsname.old,simcount,sep="")
#need old runsname
load(paste("results/",runsname.old,"/",name.old,"/dat.Rdata",sep=""),envir=e<-new.env())
pop.old=e[["pophistory"]]
pop=pop.old[[length(pop.old)]]
e=NULL;pop.old=NULL
###########################
## Running the Simulation #
###########################
pophistory=runSim(startpop=pop,years.list=years.list,
                  years.ind=years.index,N=N,duration=duration,
                  sds=sds,mutrate=mutrate,generations=numYears,traits=traits)
#Note: we've already used year 1 in initiating the pop
#Based on the value of "runType", generate the appropriate type of data.