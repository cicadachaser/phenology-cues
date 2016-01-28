#Saving results from the evolution simulations

#Turn results from list to data frame
pophist.table <-do.call(rbind.data.frame, pophistory)

setwd("results")
resultsdir=sprintf("resRun%d",runNumber)
dir.create(resultsdir,showWarnings = FALSE)
setwd(resultsdir)
write.table(pophist.table,file=paste("pophist_run",runNumber,".csv",sep=""),sep=",")
parnames=c(
  "best.precip",
  "best.temp",
  "duration",
  "N"
)
parvals=get(parnames)
meta=sprintf("%s has value %f",parnames,parvals)
sink(paste("par_values_run",runNumber,".txt",sep=""))
cat("Parameters for simulation. Weather from davis data. \n")
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

save(list=c("pophistory","years.list","years.index"),file="dat.RData")
