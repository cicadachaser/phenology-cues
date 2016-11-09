setwd("G:/Repos/phenology-cues/results/fig1/opt3-2yr-comp1listlist")
load("fig1dat-versionopt3-2yr-comp1listlist.Rdata",envir=opt<-new.env())
opt.varmeans=opt$overall.res[order(opt$overall.res$yearstd,opt$overall.res$daystd,opt$overall.res$trait),]
setwd("G:/Repos/phenology-cues/results/fig1/compare3-2yr")
load("fig1dat-versioncompare3-2yr.Rdata",envir=brute<-new.env())
overall.res=brute$overall.res
varmeans=aggregate(as.numeric(cbind(overall.res$geofit,overall.res$traitval)),
                   by=list(daystd=as.factor(overall.res$daystd),
                           yearstd=as.factor(overall.res$yearstd),
                           trait=overall.res$trait),
                   FUN=mean)
varmeans$yearstd=as.numeric(as.character(varmeans$yearstd))
varmeans$daystd=as.numeric(as.character(varmeans$daystd))
varmeans$x=as.numeric(as.character(varmeans$x))
varpts=unique(varmeans[,c("daystd","yearstd")])
varmeans=varmeans[order(varmeans$yearstd,varmeans$daystd,varmeans$trait),]
brute$varmeans=varmeans

varcompare=cbind(varmeans, opt.varmeans[,c("geofit","traitval")])
names(varcompare)[4:6]=c("geofit.old","geofit.new","traitval.new")
varcompare=cbind(varcompare,res.diff=varcompare$geofit.old+varcompare$geofit.new)
varcompare$geofit.old=-varcompare$geofit.old
varcompare=cbind(varcompare,frac.diff=varcompare$res.diff/mean(varcompare$geofit.old+varcompare$geofit.new))
names(varcompare)[8]="frac.diff"
set_wrkdir()
save(varcompare,file="varCompare_2yr.Rdata")


#Are the two lists of years identical
identical(opt$yearlistlist,brute$yearlistlist)

for(i in 1:length(opt$yearlistlist)){
  for(j in 1:length(opt$yearlistlist[[1]])){
    testopt=opt$yearlistlist[[i]][[j]]
    testbrute=brute$yearlistlist[[i]][[j]]
    print(identical(testopt,testbrute))
  }
}
