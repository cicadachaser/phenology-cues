
  require(doSNOW); require(parallel); require(doParallel)
  nClust=detectCores(all.tests=FALSE,logical=TRUE)
  c1<-makeCluster(nClust-1)
  registerDoParallel(c1)

  testfun=function(x){x^3}
  res=foreach(i=1:50) %dopar% {
    testfun(i);
    (c(i,testfun(i)));
  }

