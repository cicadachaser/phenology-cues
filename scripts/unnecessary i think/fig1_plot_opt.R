set_wrkdir()
load(file=paste("results/fig1/",runnum,"/fig1dat-version",runnum,".Rdata",sep=""))
varmeans=overall.res
varpts=unique(varmeans[,c("daystd","yearstd")])
winner=rep("init",nrow(varpts))
yearstd=daystd=rep(-10,nrow(varpts))
for(i.pt in 1:nrow(varpts)){
  cur.ind=which(varmeans$daystd==varpts$daystd[i.pt] & varmeans$yearstd==varpts$yearstd[i.pt])
  curvar=varmeans[cur.ind,]
  yearstd=
  if(sum(curvar$geofit==min(curvar$geofit))>1){
    winner[i.pt]="tie"
  }else{
    winner[i.pt]=curvar[which.max(curvar$geofit),]$trait
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

df=data.frame(x=varpts$daystd,
              y=varpts$yearstd,
              z=winner)
x11()
print(ggplot(df, aes(x,y))+
  geom_tile(aes(fill=z))+
  labs(title="Winning trait",
       x="day to day variance",
       y="year to year variance"))
dev.print(pdf,paste("winners-heatmap-run-",runnum,".pdf",sep=""))

for(i.yearstd in unique(varmeans$yearstd)){
  cur.df=varmeans[which(varmeans$yearstd==i.yearstd),]
  print(ggplot(cur.df,aes(daystd,geofit))+
    geom_line(aes(color=factor(cur.df$trait)))+
    labs(title=paste("fitness by stdev, yearstdev at",i.yearstd),
         x="day to day stdev",
         y="geometric fitness"))
  dev.print(pdf,paste("day-by-fit-yrstd",round(i.yearstd),"-run-",runnum,".pdf",sep=""))
}


for(i.daystd in unique(varmeans$daystd)){
  cur.df=varmeans[which(varmeans$daystd==i.daystd),]
  print(ggplot(cur.df,aes(yearstd,geofit))+
    geom_line(aes(color=factor(trait)))+
    labs(title=paste("fitness by stdev, daystd at",i.daystd),
         x="year to year std",
         y="geometric fitness"))
  dev.print(pdf,paste("year-by-fit-daystd",round(i.daystd),"-run-",runnum,".pdf",sep=""))
}
