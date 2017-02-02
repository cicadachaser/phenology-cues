set_wrkdir()
setwd("results/fig1/opt4-comp")
load("fig1dat-versionopt4-comp.Rdata",envir=opt<-new.env())
opt.varmeans=opt$overall.res[order(opt$overall.res$yearstd,opt$overall.res$daystd,opt$overall.res$trait),]
set_wrkdir()
setwd("results/fig1/fit4-comp")
load("fig1dat-versionfit4-comp.Rdata",envir=brute<-new.env())
overall.res=brute$overall.res
oldtrait=aggregate(as.numeric(overall.res$traitval),
                   by=list(daystd=as.factor(overall.res$daystd),
                           yearstd=as.factor(overall.res$yearstd),
                           trait=overall.res$trait),
                   FUN=mean)
oldtrait=oldtrait$x
# varmeans=aggregate(as.numeric(overall.res$geofit),
# setwd("G:/Repos/phenology-cues/results/fig1/opt3-2yr-comp1listlist")
# load("fig1dat-versionopt3-2yr-comp1listlist.Rdata",envir=opt<-new.env())
# opt.varmeans=opt$overall.res[order(opt$overall.res$yearstd,opt$overall.res$daystd,opt$overall.res$trait),]
# setwd("G:/Repos/phenology-cues/results/fig1/compare3-2yr")
# load("fig1dat-versioncompare3-2yr.Rdata",envir=brute<-new.env())
overall.res=brute$overall.res
varmeans=aggregate(as.numeric(overall.res$geofit),
                   by=list(daystd=as.factor(overall.res$daystd),
                           yearstd=as.factor(overall.res$yearstd),
                           trait=overall.res$trait),
                   FUN=mean)
varmeans$yearstd=as.numeric(as.character(varmeans$yearstd))
varmeans$daystd=as.numeric(as.character(varmeans$daystd))
varmeans$x=as.numeric(as.character(varmeans$x))
varpts=unique(varmeans[,c("daystd","yearstd")])
varmeans=cbind(varmeans,oldtrait)
varmeans=varmeans[order(varmeans$yearstd,varmeans$daystd,varmeans$trait),]
brute$varmeans=varmeans

varcompare=cbind(varmeans, opt.varmeans[,c("geofit","traitval")])
names(varcompare)[4:7]=c("geofit.old","traital.old","geofit.new","traitval.new")
varcompare=cbind(varcompare,res.diff=varcompare$geofit.old-varcompare$geofit.new)
varcompare$geofit.old=varcompare$geofit.old
varcompare=cbind(varcompare,frac.diff=varcompare$res.diff/mean(varcompare$geofit.old+varcompare$geofit.new))
names(varcompare)[9]="frac.diff"
set_wrkdir()
save(varcompare,file="varCompare_2yr.Rdata")

#Let's look at cases where they didn't match
varcompare[,c(1:3,7,8)]
varcompare[which(abs(varcompare[,"res.diff"])>10^-10),]


=======
names(varcompare)[4:6]=c("geofit.old","geofit.new","traitval.new")
varcompare=cbind(varcompare,res.diff=varcompare$geofit.old+varcompare$geofit.new)
varcompare$geofit.old=-varcompare$geofit.old
varcompare=cbind(varcompare,frac.diff=varcompare$res.diff/mean(varcompare$geofit.old+varcompare$geofit.new))
names(varcompare)[8]="frac.diff"
set_wrkdir()
save(varcompare,file="varCompare_2yr.Rdata")

>>>>>>> 8f45bd1e6a5c8f2a0cfce45534c0605b131e8b45

#Are the two lists of years identical
identical(opt$yearlistlist,brute$yearlistlist)

for(i in 1:length(opt$yearlistlist)){
  for(j in 1:length(opt$yearlistlist[[1]])){
    testopt=opt$yearlistlist[[i]][[j]]
    testbrute=brute$yearlistlist[[i]][[j]]
    print(identical(testopt,testbrute))
  }
}
