#Saving results from the evolution simulations

#Turn results from list to data frame
pophist.table <-do.call(rbind.data.frame, pophistory)

setwd("results")
resultsdir=sprintf("%s/resRun%s",runsname,runName)
unlink(resultsdir,recursive = TRUE)
dir.create(resultsdir,showWarnings = FALSE)
setwd(resultsdir)
write.table(pophist.table,file=paste("pophist_run",runName,".csv",sep=""),row.names=FALSE,sep=",")
parnames=c(
  "best.precip",
  "best.temp",
  "duration",
  "N",
  "numYears",
  "runType"
)
parvals=get(parnames)
meta=sprintf("%s has value %f",parnames,parvals)
sink(paste("par_values_run",runName,".txt",sep=""))
cat("Parameters for simulation. Weather from davis data. \n")
cat(" Run type = ")
cat(runType)
cat("\n")
cat("traits = ")
cat(traits)
cat("\n")
for(i in meta){cat(paste(i,"\n"))}
cat(" sds= \n")
print(sds)
cat("\n mutrate=\n")
print(mutrate)
cat("\n start=\n")
print(start)
cat("\n years.index=\n")
print(years.index)
sink()

save(list=ls(all.names=TRUE),file="dat.RData")
