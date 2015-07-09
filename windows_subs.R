fitness<-function(t.start,t.duration,W){ # fitness is the sum of W over the lifespan
  #Function for giving fitness of individuals based on their start time, duration, and the W.
  #Inputs:
  #  t.start: vector of start times, one per individual
  #  t.duration: vector of duration times, one per individual
  #  W: vector of the `goodness of environment'
  #Returns:
  #  res: vector of the fitnesses of each individual
  #
  ind.start<-round((t.start*length(W))/12)
  ind.duration<-round((t.duration*length(W))/12)
  ind.end<-ind.start+ind.duration
  wtemp=rep(W,2) #! make a double-long W to simplify summing fitness.
  res=rep(-1000,times=length(ind.start)); #create a place to store results, fill with arbitrary identifable number
  for(i in 1:length(ind.start)){
    res[i]=sum(wtemp[ind.start[i]:ind.end[i]])
  }
  return(res)
}

selection1<-function(newpop,W){
  #Function for carrying out `hard selection' - individuals that have fitness lower than the mid-range do not reproduce.
  #Inputs:
  #  newpop: data frame with two columns (t.start and t.duration), and a row for every individual
  #  W: vector of goodness of environment, to pass to fitness() function
  #Returns:
  #  newpop: a matrix with the traits and fitness of each individual. Has 6 columns - t.start,
  #     t.duration, Wi, Ws, Wp, and Wnum
  #
  newWi<-mapply(fitness,newpop$t.start,newpop$t.duration,MoreArgs=list(W=W)) #MoreArgs argument lets us pass extra bits to mapply
  newWs<-2*((newWi)-mean(range(newWi)))/(max(newWi)-min(newWi)) #rescaled between -1 and 1, centered on the mid-range
  newWsurv<-newWs*(newWs>0) #Wsurv: all individuals with Ws<0 have zero fitness #!vectorized this: ~10x faster, maybe a bit less readible.
  newWp<-newWsurv/sum(newWsurv) #Wp is the proportional fitness after mortality
  #newWnum<-round(N*newWp) #Wnum is the integer number of offspring for each individual, population maintained at N
  newWnum=t(rmultinom(1,size=N,prob=newWp)) #To avoid potential rounding weirdness, had individuals assigned via the multinomial distribution
  newpop<-cbind(newpop,newWi,newWs,newWp,newWnum)
  colnames(newpop)<-c("t.start","t.duration","Wi","Ws","Wp","Wnum")
  return(newpop)
}

#selection2 "soft selection" all individuals reproduce, with variable fitness
#function expands a two-col dataframe with t. start and t. duration into a 6-col with fitness measures
selection2<-function(newpop,W){
  #Function for carrying out `soft selection' - all individuals reproduce, with variable fitness.
  #Inputs:
  #  newpop: data frame with two columns (t.start and t.duration), and a row for every individual
  #  W: vector of goodness of environment, to pass to fitness() function
  #Returns:
  #  newpop: a matrix with 6 columns - t.start, t.duration, Wi, Ws, Wp, and Wnum
  #
  newWi<-mapply(fitness,newpop$t.start,newpop$t.duration,MoreArgs=list(W=W))
  newWs<-(newWi-min(newWi))/(max(newWi)-min(newWi)) #rescaled between 0 and 1, centered on the mid-range
  newWsurv<-newWs*(newWs>0) #newWsurv: all individuals survive (some may have zero fitness, none have neg fitness)
  newWp<-newWsurv/sum(newWsurv) #Wp is the proportional fitness after mortality
#   newWnum<-round(N*newWp) #Wnum is the integer number of offspring for each individual, population maintained at N
  newWnum=(rmultinom(1,size=N,prob=newWp)) #To avoid potential rounding weirdness, had individuals assigned via the multinomial distribution
  newpop<-cbind(newpop,newWi,newWs,newWp,newWnum)
  colnames(newpop)<-c("t.start","t.duration","Wi","Ws","Wp","Wnum")
  return(newpop)
}

#selection3 "tunable selection" survival cut-off can be determined by the user
#function expands a two-col dataframe with t. start and t. duration into a 6-col with fitness measures
selection3<-function(newpop,cutoff,W){
  #Function for carrying out selection with a custom cutoff - all individuals reproduce, with variable fitness.
  #Inputs:
  #  newpop: data frame with two columns (t.start and t.duration), and a row for every individual
  #  cutoff: fitness below which individuals don't reproduce.
  #  W: vector of goodness of environment, to pass to fitness() function
  #Returns:
  #  newpop: a matrix with 6 columns - t.start, t.duration, Wi, Ws, Wp, and Wnum
  #
  newWi<-mapply(fitness,newpop$t.start,newpop$t.duration,MoreArgs=list(W=W))
  newWs<-qunif((newWi-min(newWi))/(max(newWi)-min(newWi)),min=-1,max=1) #rescaled between -1 and 1, centered on the median
  newWsurv<-newWs*(newWs>max(0,cutoff)) #Wsurv: all individuals with Ws<cutoff (or zero) have zero fitness
  newWp<-newWsurv/sum(newWsurv) #Wp is the proportional fitness after mortality
#   newWnum<-round(N*newWp) #Wnum is the integer number of offspring for each individual, population maintained at N
  newWnum=(rmultinom(1,size=N,prob=newWp)) #To avoid potential rounding weirdness, had individuals assigned via the multinomial distribution
  newpop<-cbind(newpop,newWi,newWs,newWp,newWnum)
  colnames(newpop)<-c("t.start","t.duration","Wi","Ws","Wp","Wnum")
  return(newpop)
}


#! Added tuning parameter
#! Made function take vectors, return a matrix of new start and new duration
mutation<-function(t.start,t.duration,sd.start=0.5,sd.dur=0.5){
  #Function that creates random offspring with variable t.start and t.duration values
  #Inputs:
  #  t.start: vector of starting times, one for each individual
  #  t.duration: vector of durations, one for each individual
  #  sd.start: standard deviation for mutation of start time
  #  sd.dur: standard deviation for mutation of duration
  #Returns
  #  2-dimensional matrix of new starting times and new durations.
  #
  new.t.start<-t.start+round(rnorm(length(t.start),mean=0,sd=sd.start),1) #Gaussian random mutation; adjust sd for "mutation rate"
  #consider not rounding? The initial values aren't integers
  new.t.start<-new.t.start*(new.t.start>0) #return negative values to zero
  new.t.start<-((new.t.start-1) %% 12)+1 #wrapping t.start - using modulo operator for speed. It's a cool function, but gives zeros unless you offset with -1 and then +1
  new.t.duration<-t.duration+round(rnorm(length(t.start),mean=0,sd=sd.dur),1) #Gaussian random mutation; adjust sd for "mutation rate"
  #consider not rounding? The initial values aren't integers
  new.t.duration<-new.t.duration*(new.t.duration>0) #return near zero or negative values to 0.1
  return(cbind(new.t.start,new.t.duration))
}

reproduction<-function(pop){
  #Function for handling reproduction
  #Inputs:
  #  pop: population data frame generated by selection() function
  #Returns:
  #  Next generation as a data frame.
  repop<-pop[pop$Wnum>0,]
  expandpop<-data.frame(t.start=rep(0,N),t.duration=rep(0,N))
  ind=1
  for(i in 1:nrow(repop)){
    for (j in 1:repop[i,6]) {
      expandpop[ind,]<-c(repop[i,1],repop[i,2])
      ind=ind+1
    }
  }
  newpop = mutation(expandpop[,1],expandpop[,2]) #remove initial placeholder row
  colnames(newpop)<-c("t.start","t.duration")
  return(data.frame(newpop))
}

#environmental variation function
envar<-function(y1,y2,month,y1.opt,y2.opt,sd=.5){
  #Function for adding noise to environment and generating a fine-grain resolution environmental `goodness' metric W
  #Inputs:
  #  y1: vector of monthly averages of temp
  #  y2: vector of monthly average precip
  #  y1.opt: optimal temp value
  #  y2.opt: optimal precip value
  #  sd:standard deviation defining the amount of variation between years, defaults to 0.5
  #Returns
  #  W: vector of immediate fitness-gains for all points in time.
  #
  newy1<-y1+round(rnorm(1,mean=0,sd=sd),0) # add or subtract 1 from each month (as in `wet yr')
  newy2<-y2+round(rnorm(1,mean=0,sd=sd),0)
  #! Some alternative ways to add variation: 
  #!   Don't need round(), as temp and precip are continuous variables.
  #!   Can use sample() with a vector of probabilities to specify likelihood of 
  #!     various different increases/decreases (if want discrete variations).
  #!   Could have offset years, where monthly averages occur late or early
  model1<-loess(newy1~month, span=.35); #create model for smoothing the year across months
  xv1<-seq(0,12,0.0001)
  yv1<-predict(model1,data.frame(month=xv1))# create smooth year
  model2<-loess(newy2~month, span=.35);
  xv2<-seq(0,12,0.0001)
  yv2<-predict(model2,data.frame(month=xv2))
  W1r<-dnorm(yv1,mean=y1.opt,sd=y1.opt) #Fitness as a function of temp - distance from optimal temp over time.
  W1<-(W1r-min(W1r))/(max(W1r)-min(W1r)) #rescaled between 0 to 1
  W2r<-dnorm(yv2,mean=y2.opt,sd=y2.opt) #Fitness as a function of precip
  W2<-(W2r-min(W2r))/(max(W2r)-min(W2r)) #rescaled between 0 to 1
  Wr<-W1*W2 #raw fitness - multiply fitness by precip and temp
  W<-2*(Wr-mean(range(Wr)))/(max(Wr)-min(Wr)) #defining the combined fitness landscale, rescaled between -1 and 1 to prevent long-lived strategies
  return(W)
}


runSim<-function(pop,y1,y2,month,y1.opt,y2.opt,generations=24,graphics=FALSE){
  #Function that actually runs the simulation (calling the other functions above)
  #Inputs:
  #  pop: initial population
  #  y1: baseline monthly mean temperatures
  #  y2: baseline monthly mean precip
  #  month: vector of months corresponding to y1, y2
  #  y1.opt: optimal temperature
  #  y2.opt: optimal precip
  #  generations: number of generations to simulate. Defaults to 24
  #  graphics: boolean, defaults to false. If true, carry out some plotting operations
  #Returns
  #  pophistory: list where each element represents the full population data frame for each generation
  #
  pophistory<-list(pop) #initialize the population history
  for(g in 1:generations){
    #reproduction
    newpop<-reproduction(pop)
    #inter-annual variation
    W<-envar(y1=y1,y2=y2,month=month,y1.opt=y1.opt,y2.opt=y2.opt)
    #selection
    newpop<-selection2(newpop,W)
    pophistory[[g+1]]<-newpop
    pop<-newpop
    #add arrows for each individual
    xlines<-pop$t.start
    ylines<-2*(newpop$Ws-mean(c(min(pop$Ws),max(pop$Ws))))/(max(pop$Ws)-min(pop$Ws)) #rescale height of arrow between -1 and 1
    if(graphics==TRUE){
      dev.new(height=4, width=8)
      par(mar=c(4, 4, 1, 1) + 0.1)
      plot(xv2,W,type="l",xaxp=c(0, 12, 12))
      abline(h=0,lty=3)
      arrows(xlines,ylines,x1=(pop$t.start+pop$t.duration),lty=1,length=.1)
    }
  }
  return(pophistory)
}