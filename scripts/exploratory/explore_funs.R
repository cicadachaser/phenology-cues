set_wrkdir<-function(){
  #function for setting working directory to the right place given the current computer/user
  if(Sys.getenv("USERNAME")=="Collin" || Sys.getenv("USERNAME")=="collin" || Sys.getenv("USERNAME")=="Collin.work"){ #If it's collin
    if(Sys.info()["nodename"]=="DESKTOP-D6QSU8F"){
      setwd("G:\\Repos\\phenology-cues") #desktop
    }else{
      setwd("C:\\Repos\\phenology-cues") #desktop
    }
  }else{
    if(Sys.getenv("COMPUTERNAME")=="ENT-YANG01"){
      setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")#desktop
    }else{
      setwd("C:\\Users\\louie\\Documents\\GitHub\\phenology-cues")} #laptop
  }
}

exp_fitfun <- function(fit.file="standgauss.R",
                       min.temp=0, #min temp value to plot
                       max.temp=50, #max temp value to plot
                       min.other=0, #min "other" value to plot
                       max.other=50, #max "other" value to plot
                       numpts=500, #Resolution of image - number points in each X and Y dimension
                       climate.fun="yeargen.davis", #for plotting some random actual data.
                       n.plotyears=3 #number of random years to plot data from
                       #  set to FALSE if you want to not include points on the plot
){
  #Function for producing plot of fitness function
  # Note that it uses temp and "other" as the two traits used by the fitness function
  # 'other' might be precip, or it might be our upcoming "moisture" or any other climate factor
  #COMMENT THIS FUNCTION!

  #Make case file, fitfuns.R
  require(fields) #for use in image.plot
  set_wrkdir()
  source(paste("fitcurve/",fit.file,sep=""))
  tempmat=matrix(seq(min.temp,max.temp,length=numpts),ncol=numpts,nrow=numpts,byrow = TRUE)
  #matrix of temp values
  othermat=matrix(seq(min.temp,max.temp,length=numpts),ncol=numpts,nrow=numpts,byrow = FALSE)
  #vector of other values
  tempv=tempmat[1:(numpts^2)] #make vector of temp values to test
  otherv=othermat[1:(numpts^2)] #make vector of other values to test
  pseudo.year=data.frame(temp=tempv,placeholder=otherv) #Turn our vectors into a fake year to feed the fit_fn
  names(pseudo.year)[2] = other.name #replaces "placeholder" with the specified other.name
  fit=fit_fn(pseudo.year,other.name=other.name) #calculate fitness
  fit.mat=matrix(fit,numpts,numpts) #convert fitness to matrix
  image.plot(x=seq(min.other,max.other,length=numpts), #make a heatmap of fitness in the given range of temp and other
             y=seq(min.temp,max.temp,length=numpts),
             z=fit.mat,
             xlab=other.name,
             ylab="Temperature",
             main="Fitness gained per day",
             cex.lab=1.6,cex.main=1.6,cex.axis=1.5
  )
  if(climate.fun!=FALSE){ #if climate.File is a string
    source("scripts/windows_subs.R")
    yearfun=get(climate.fun) #find the function represented tby the climate.fun string
    list.yr=yearfun(best.other,sd.other,best.temp,sd.temp,baseTemp=0,other.name=other.name)[[1]]
    #Pull in the appropriate year data
    samp.years=sample(1:length(list.yr),n.plotyears) #choose n.plotyears random years
    oneyr=do.call("rbind", list.yr[samp.years]) #create a single data frame
    points(jitter((oneyr[,other.name])),jitter(oneyr$temp),lwd=2) #plot the points
  }
}

exp_moist <- function(dat.file="davisDat.Rdata",decay=.2,numyears=3){
  ##IN PROGRESS
  #NEED TO:
  #  FINISH
  #  COMMENT
  #  WRITE CORRESPONDING USE SCRIPT
  set_wrkdir()
  source(paste("fitcurve/",fit.file,sep=""))
  source("scripts/windows_subs.R")
  yearstuff=yeargen(dat.file=dat.file,best.other,sd.other,best.temp,sd.temp,baseTemp=0,other.name=other.name,decay=decay)
  list.yr=yearstuff[[1]]
  moiststore=NULL
  for(i.year in 1:length(list.yr)){moiststore=cbind(moiststore,list.yr[[i.year]]$moist[1:365])}
  meanmoist=apply(moiststore,1,mean)
  plot(1:365,meanmoist,col="black",lwd=3,lty=2,type='l',ylim=c(0,max(meanmoist)*2.5))
  for(i in 1:numyears){
    points(1:365,list.yr[[i]]$moist[1:365],type='l',col=i)
  }
}

