
  require(doSNOW); require(parallel); require(doParallel)
  nClust=detectCores(all.tests=FALSE,logical=TRUE)
  c1<-makeCluster(min(nClust-1,5))
  registerDoParallel(c1)

  testfun=function(x){x^3}
  foreach(i=1:3) %dopar% testfun(i)
