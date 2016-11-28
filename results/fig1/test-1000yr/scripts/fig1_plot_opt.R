set_wrkdir()
load(file=paste("results/fig1/",runnum,"/fig1dat-version",runnum,".Rdata",sep=""))
varmeans=overall.res
varpts=unique(varmeans[,c("daystd","yearstd")])
winner=rep("init",nrow(varpts))
yearstds=daystds=rep(-10,nrow(varpts))
for(i.pt in 1:nrow(varpts)){
  cur.ind=which(varmeans$daystd==varpts$daystd[i.pt] & varmeans$yearstd==varpts$yearstd[i.pt])
  curvar=varmeans[cur.ind,]
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
#Now print in black and white
print(ggplot(df, aes(x,y))+
  geom_tile(aes(fill=z))+
    scale_fill_grey()+
  labs(title="Winning trait",
       x="day to day variance",
       y="year to year variance"))
dev.print(pdf,paste("winners-heatmap-bw-run-",runnum,".pdf",sep=""))

for(i.yearstds in unique(varmeans$yearstd)){
  cur.df=varmeans[which(varmeans$yearstd==i.yearstds),]
  print(ggplot(cur.df,aes(daystd,geofit))+
    geom_line(aes(color=factor(cur.df$trait)))+
    labs(title=paste("fitness by stdev, yearstdsev at",i.yearstds),
         x="day to day stdev",
         y="geometric fitness"))
  dev.print(pdf,paste("day-by-fit-yrstd",round(i.yearstds),"-run-",runnum,".pdf",sep=""))
}


for(i.daystds in unique(varmeans$daystd)){
  cur.df=varmeans[which(varmeans$daystd==i.daystds),]
  print(ggplot(cur.df,aes(yearstd,geofit))+
    geom_line(aes(color=factor(trait)))+
    labs(title=paste("fitness by stdev, daystds at",i.daystds),
         x="year to year std",
         y="geometric fitness"))
  dev.print(pdf,paste("year-by-fit-daystds",round(i.daystds),"-run-",runnum,".pdf",sep=""))
}
