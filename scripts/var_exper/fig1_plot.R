set_wrkdir()
load(file=paste("results/fig1/",runnum,"/fig1dat-version",runnum,".Rdata",sep=""))
varmeans=aggregate(as.numeric(overall.res$geofit),
                   by=list(dayvar=as.factor(overall.res$dayvar),
                           yearvar=as.factor(overall.res$yearvar),
                           trait=overall.res$trait),
                   FUN=mean)
varmeans$yearvar=as.numeric(as.character(varmeans$yearvar))
varmeans$dayvar=as.numeric(as.character(varmeans$dayvar))
varpts=unique(varmeans[,c("dayvar","yearvar")])
winner=rep("init",nrow(varpts))
for(i.pt in 1:nrow(varpts)){
  cur.ind=which(varmeans$dayvar==varpts$dayvar[i.pt] & varmeans$yearvar==varpts$yearvar[i.pt])
  curvar=varmeans[cur.ind,]
  if(sum(curvar$x==min(curvar$x))>1){
    winner[i.pt]="tie"
  }else{
    winner[i.pt]=curvar[which.min(curvar$x),]$trait
  }
}
# x11()
# m=matrix(c(1,1,1,2),nrow=1)
# layout(m)
# win.levels=unique(winner)
# col=c("grey","red","blue","green")[1:length(win.levels)]
# image(matrix(winner==win.levels[1],ncol=numpts),
#       x=sort(as.numeric(levels(varpts$dayvar))),
#       y=sort(as.numeric(levels(varpts$yearvar))),
#       col=c(NA,col=col[1]),
#       xlab="Day-to-day variance",
#       ylab="year-to-year variance")
#
# for(i in 2:length(win.levels)){
#   image(matrix(winner==win.levels[i],ncol=numpts),
#       x=sort(as.numeric(levels(varpts$dayvar))),
#       y=sort(as.numeric(levels(varpts$yearvar))),
#         col=c(NA,col[i]),add=TRUE)
# }
# par(xpd=TRUE)
# legend(legend=win.levels, col=col, title="Best Trait",xpd = TRUE)

########################
#Making the plots
########################
#First, make heat map of winners
set_wrkdir()
setwd(paste("results/fig1/",runnum,sep=""))

require(ggplot2)

df=data.frame(x=varmeans$dayvar,
              y=varmeans$yearvar,
              z=winner)
x11()
ggplot(df, aes(x,y))+
  geom_tile(aes(fill=z))+
  labs(title="Winning trait",
       x="day to day variance",
       y="year to year variance")
dev.print(pdf,paste("winners-heatmap-good-run-",runnum,".pdf",sep=""))

for(i.yearvar in unique(varmeans$yearvar)){
  cur.df=varmeans[which(varmeans$yearvar==i.yearvar),]
  print(ggplot(cur.df,aes(dayvar,x))+
    geom_line(aes(color=factor(trait)))+
    labs(title=paste("fitness by variance, yearvar at",i.yearvar),
         x="day to day variance",
         y="geometric fitness"))
  dev.print(pdf,paste("day-by-fit-yrvar",round(i.yearvar),"-run-",runnum,".pdf",sep=""))
}


for(i.dayvar in unique(varmeans$dayvar)){
  cur.df=varmeans[which(varmeans$dayvar==i.dayvar),]
  print(ggplot(cur.df,aes(yearvar,x))+
    geom_line(aes(color=factor(trait)))+
    labs(title=paste("fitness by variance, dayvar at",i.dayvar),
         x="year to year variance",
         y="geometric fitness"))
  dev.print(pdf,paste("year-by-fit-dayvar",round(i.dayvar),"-run-",runnum,".pdf",sep=""))
}
