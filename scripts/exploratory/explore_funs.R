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


fittot_fun <- function(years.list,duration=10,lag=1){
  #For adding fit.tot to a years.list file
  require(zoo)
  for(i.year in 1:length(years.list)){
    fit.tot=c(rollapply(c(years.list[[i.year]]$fit.daily,rep(0,duration-1)),duration,by=1,sum))
    fit.tot=c(fit.tot[-(1:lag)],rep(0,lag))#for lag of 1
    years.list[[i.year]]=cbind(years.list[[i.year]],fit.tot)
  }
  return(years.list)
}

exp_fitfun <- function(fit.file="standgauss.R",
                       min.temp=0, #min temp value to plot
                       max.temp=50, #max temp value to plot
                       min.other=0, #min "other" value to plot
                       max.other=50, #max "other" value to plot
                       numpts=500, #Resolution of image - number points in each X and Y dimension
                       dat.file="davisDat.R", #for plotting some random actual data.
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
    list.yr=yeargen(dat.file=dat.file, fit.parms=fit.parms)[[1]]
    #Pull in the appropriate year data
    samp.years=sample(1:length(list.yr),n.plotyears) #choose n.plotyears random years
    oneyr=do.call("rbind", list.yr[samp.years]) #create a single data frame
    points(jitter((oneyr[,other.name])),jitter(oneyr$temp),lwd=2) #plot the points
  }
}

exp_moist <- function(dat.file="davisDat.Rdata", #Whiche climate file to use
                      fit.file="standgauss.R", #Which fitness file to use
                      decay=.2, #The decay parameter for calculating moisture
                      baseTemp=0, #The base temperature to use for calucating cumulative temperature
                      numyears=3, #the number of years to plot for the individual plots
                      other.name="moist", #The second climate variable to use for fitness (in addition to temp)
                      duration=10, #duration used in calculating total fitness
                      lag=1){ #lag used in calculating total fitness
  set_wrkdir()
  source(paste("fitcurve/",fit.file,sep=""))
  source("scripts/windows_subs.R")
  yearstuff=yeargen(dat.file=dat.file, fit.parms=fit.parms,decay=decay,baseTemp=baseTemp,other.name=other.name)
  list.yr=yearstuff[[1]] #grab the years.list part of yearstuff
  list.yr=fittot_fun(years.list=list.yr,duration=duration,lag=lag) #add total fitness
  #Create average moisture measure
  moiststore=NULL
  fitstore=NULL
  for(i.year in 1:length(list.yr)){
    moiststore=cbind(moiststore,list.yr[[i.year]]$moist[1:365])
    fitstore=cbind(fitstore,list.yr[[i.year]]$fit.tot[1:365])
  }
  meanmoist=apply(moiststore,1,mean)
  meanfit=apply(fitstore,1,mean)
  #Plotting
  par(mfcol=c(2,2),oma=c(0,0,2,0))
  plot(1:365,meanmoist,col="black",lwd=2,lty=1,type='l',
       main="Mean moisture across all years",
       ylab="Moisture",
       xlab="day")
  plot(1:365,meanfit,col="black",lwd=2,lty=1,type='l',
       main="Mean fitness across all years",
       ylab="Fitness",
       xlab="day")
  plot(1:365,list.yr[[1]]$moist[1:365],col="black",lwd=1,lty=1,type='l',
       main=paste("Moisture for first",numyears,"years"),
       ylab="Moisture",
       xlab="day")
  if(numyears>1){
    for(i in 2:numyears){
      points(1:365,list.yr[[i]]$moist[1:365],type='l',col=i)
    }
  }
  plot(1:365,list.yr[[1]]$fit.tot[1:365],col="black",lwd=1,lty=1,type='l',
       main=paste("fitness for first",numyears,"years"),
       ylab="fitness",
       xlab="day")
  if(numyears>1){
    for(i in 2:numyears){
      points(1:365,list.yr[[i]]$fit.tot[1:365],type='l',col=i)
    }
  }
  mtext(paste("Using: ", dat.file,", ", fit.file,", decay = ",decay,sep=""), outer = TRUE, cex = 1.5)
}



exp_corr_plot <- function(df.yr,interest,textpos="top"){
  #function for plotting within exp_corr
  vals.interest=df.yr[,interest]
  out=lm(df.yr$fit.tot ~ vals.interest)
  xjit.amount=sort(abs(unique(vals.interest)))[which(sort(abs(unique(vals.interest)))>1*10^-10)[1]] #ugly, but finds smallest non-zero difference between x values
  plot(jitter(vals.interest,amount=.4*xjit.amount),jitter(df.yr$fit.tot),pch=".",
       main=paste("Correlation of fitness with", interest),
       ylab="Total fitness",
       xlab=interest,
       cex.lab=1.5)
  abline(out,col='blue',lwd=2)
  #Add information about our line
  legend(x=textpos,legend=paste("r-squared =", round(summary(out)$r.squared,digits=4)),
         bty="n",cex=1.5)
}

exp_corr <- function(dat.file="davisDat.Rdata", #Whiche climate file to use
                     traits=c("temp","moist","day","cutemp"),
                     text.pos=FALSE, #For custom position of legends, replace with a vector of position strings
                     fit.file="standgauss.R", #Which fitness file to use
                     decay=.2, #The decay parameter for calculating moisture
                     baseTemp=0, #The base temperature to use for calucating cumulative temperature
                     numyears=3, #the number of years to plot for the individual plots
                     other.name="moist", #The second climate variable to use for fitness (in addition to temp)
                     duration=10, #duration used in calculating total fitness
                     lag=1){ #lag used in calculating total fitness
  set_wrkdir()
  source(paste("fitcurve/",fit.file,sep=""))
  source("scripts/windows_subs.R")
  yearstuff=yeargen(dat.file=dat.file, fit.parms=fit.parms,decay=decay,baseTemp=baseTemp,other.name=other.name)
  list.yr=fittot_fun(years.list=yearstuff[[1]],duration=duration,lag=lag) #grab the years.list part of yearstuff
  df.yr=do.call(rbind.data.frame,list.yr)
  numrows=floor(sqrt(length(traits)))
  numcols=ceiling(sqrt(length(traits)))
  par(mfrow=c(numrows,numcols),oma=c(0,0,2,0))
  text.position=rep("top",length(traits))
  if(text.pos[1]!=FALSE) text.position=text.pos
  for(i in 1:length(traits)){
      exp_corr_plot(df.yr=df.yr,interest=traits[i], textpos=text.position[i])
  }
  mtext(paste("Using: ", dat.file,", ", fit.file,", decay = ",decay,sep=""), outer = TRUE, cex = 1.5)
}
