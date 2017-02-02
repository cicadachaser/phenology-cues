#Updated through "fitness" function

emerge_sub<-function(x){min(c(which(x), 366))} #function for use in "apply" within emergence()

emergence<-function(year, newpop, traits){
  #Function for calculating the emergence day of the given individual in the given year.
  #  Calculates emergence value as E= b.day*day+b.temp*temp+b.precip*precip
  # Then finds the first day when the calculated E is greater than 100 (100 chosen for arbitrary convenience)
  # Inputs:
  #  Indiv: (individual) has three important attributes: $b.day, $b.temp, $b.precip
  #  year: current year data. Includes columns $day, $tmax, $precip
  # Output:
  #  day (Julian) of emergence
  # newpop[newpop==0] = 10^10
  trans = 0*newpop #trans will hold the tranformed trait values for linear combination and 0s for unused traits
  b.traits = sprintf("b.%s", traits)
  trans[, b.traits] = 100/newpop[, b.traits]
  E = trans$b.day %*% t(year$day) +
    trans$b.temp %*% t(year$temp) +
    trans$b.precip %*% t(year$precip) +
    trans$b.cutemp %*% t(year$cutemp) +
    trans$b.cuprecip %*% t(year$cuprecip) +
    trans$b.daysq %*% t(year$daysq) +
    trans$b.tempsq %*% t(year$tempsq) +
    trans$b.precipsq %*% t(year$precipsq) +
    trans$b.cutempsq %*% (t(year$cutempsq)^2) +
    trans$b.cuprecipsq %*% t(year$cuprecipsq)
  emerge = apply(E>100, 1, emerge_sub)
  return(emerge) #find the first day where emergence value is greater than 100 (or the last day of the year)
}
fitness<-function(year, newpop, duration, traits){
  #Function for giving fitness of individuals based on their start time, duration, and the W.
  # fitness is the sum of W over the lifespan
  #FOR SIMPLICITY, ASSUMING END OF YEAR MEANS DEATH. CHANGE IF APPROPRIATE.
  #Inputs:
  #  year: data frame of climate and fitness information for current year. Includes $fit.daily column
  #  newpop: matrix of individuals, with each row corresponding to an individual, rows $b.const, $b.day, $b.temp, $b.precip
  #  duration: number of days organism is emerged
  #Returns:
  #  fit: vector of the fitnesses of each individual
  #
  evect = emergence(year = year, newpop = newpop, traits = traits)
  #Create a vector with the "total fitness you experience if you emerge on this day" values.
  # (using rollapply function from "zoo" library to make this fast)
  # Adding a 0 at the end of the vector: if you didn't emerge in the normal year, you "emerge" on day 366 which has zero fitness. ie if you don't emerge you die.
  fitVals = c(rollapply(c(year$fit.daily, rep(0, duration-1)), duration, by = 1, sum), rep(0, lag+1))
  fit = fitVals[evect+lag] #pmin is "parallel min"
  return(data.frame(fit = fit, emerge = evect))
}

selection<-function(newpop, duration, year, N, traits){
  #Function for carrying out `soft selection' - all individuals reproduce, with variable fitness.
  #Inputs:
  #  newpop: data frame with the current population trait values; each row is an individual
  #  year: data on the year, including a column for daily fitness
  #  duration: lifespan of organisms
  #  N: number of individuals in the population
  #Returns:
  #  newpop: a matrix with the current population traits, plus raw fitness (Wi),
  #     rescaled fitness(Ws), proportional fitness after mortality(Wp), and number of offspring (Wnum)
  #
  out = fitness(year = year, newpop = newpop, duration = duration, traits = traits)
  newWi = out$fit
  newWnum = (rmultinom(1, size = N, prob = newWi)) #To avoid potential rounding weirdness, had individuals assigned via the multinomial distribution
  init.colnames = colnames(newpop)
  newpop<-cbind(newpop, out$emerge, newWi, newWnum)
  colnames(newpop)<-c(init.colnames, "emerge", "Wi", "Wnum")
  return(newpop)
}

mutation<-function(poptraits, sds, mutrate, N, fattail){
  #Function that creates random offspring with variable t.start and t.duration values
  #Inputs:
  #  poptraits: a matrix of just the traits for the current population
  #  sds: 1-d data frame of standard deviations for mutations of the various traits. Has values $const, $day, $temp, $precip
  #  mutrate: 1-d data frame for PROBABILITY of mutation for each of the traits. Has values $const, $day, $temp, $precip
  #Returns
  #  2-dimensional matrix of the new (post-mutation) traits of the population
  #
  mat.runif = matrix(runif(length(sds)*N), nrow = N, ncol = length(sds)) #generate matrix of random uniform numbers for testing
  #make a data frame with N rows, each row being a duplicate of the mutrate list
  mutframe = data.frame(t(unlist(mutrate))) #one row
  test.mutate = mutframe[rep(1, N), ] #sneaky way to make N rows
  mat.mutate = mat.runif<test.mutate # which traits of which individuals mutated?
  #This is a little dicey, but when executed correctly (like here), allows me to make our entire rnorm matrix in one shot.
  #It relies on the INCREDIBLY STUPID FACT that R recycles vectors when needed. I don't like it, but the code works and is quick.
  #Using cauchy now for fat tails
  vals.mutate = matrix(rnorm(n = N*length(sds), mean = 0, sd = unlist(sds)), N, length(sds), byrow = TRUE)
  if(fattail==TRUE)
  vals.mutate = matrix(rcauchy(n = N*length(sds), location = 0, scale = unlist(sds)), N, length(sds), byrow = TRUE)
  colnames(vals.mutate)<-sprintf("b.%s", names(sds))
  poptraits = poptraits+vals.mutate*mat.mutate #Take current population, add mutations only for individuals and traits that mutated.
  return(poptraits)
}

reproduction<-function(pop){
  #THIS FUNCTION CAN PROBABLY BE MADE MUCH FASTER
  #Function for handling reproduction
  #Inputs:
  #  pop: population data frame generated by selection() function. Includes Wnum
  #Returns:
  #  Next generation as a data frame.
  repop<-pop[pop$Wnum>0, ]
  nameslist = sprintf("b.%s", names(sds))
  expandpop<-data.frame(t(rep(0, length(sds))))[rep(1, N), ]
  colnames(expandpop)<-nameslist
  ind = 1
  for(i in 1:nrow(repop)){
    for (j in 1:repop$Wnum[i]) {
      expandpop[ind, ]<-repop[i, !names(repop) %in% c("emerge", "Wi", "Wnum")]
      ind = ind+1
    }
  }
  return(expandpop)
}
runSim<-function(startpop, years.list, years.ind, N, duration, sds, mutrate, generations, traits, fattail, graphics = FALSE){
  #Function that actually runs the simulation (calling the other functions above)
  #Inputs:
  #  startpop: initial population
  #  years.list: List of dataframes for daily information on each year (MUST INCLUDE DAILY FITNESS)
  #  years.ind: vector of indices for the year to use for each generation.
  #  N: number of individuals in the population
  #  duration: number of days all individuals is in an emerged state
  #  sds: standard deviation for the distribution of mutation sizes for each trait
  #  mutrate: probability of mutation for each trait (per individual).
  #  generations: number of generations to simulate.
  #  graphics: boolean, defaults to false. If true, carry out some plotting operations
  #Returns
  #  pophistory: list where each element represents the full population data frame for each generation
  #
  pop = startpop
  pophistory<-list(cbind(startpop, gen = rep(1, N))) #initialize the population history
  for(g in 2:generations){
    #reproduction
    cur.year = years.list[[years.ind[g]]]
    newpop<-reproduction(pop = pop)
    newpop<-mutation(poptraits = newpop, sds = sds, mutrate = mutrate, N = N, fattail = fattail)
    newpop<-selection(newpop = newpop, duration = duration, year = cur.year, N = N, traits = traits)
    pophistory[[g]]<-cbind(newpop, gen = rep(g, N))
    pop<-newpop
  }
  return(pophistory)
}

######################
# Analysis Functions #
######################

#####
actTraitVals<-function(pophistory, numYears, N){
  traitslist = sprintf("b.%s", traits)
  coef.indiv = matrix(data = 0, ncol = (3+length(traitslist)), nrow = N*numYears,
                    dimnames = list(NULL, c("gen", traitslist, "relfit", "emerge")))
  ind = 1
  for(i.gen in 1:numYears){
    curhist = pophistory[[i.gen]]
    coef.indiv[ind:(ind+N-1), "gen"] = rep(i.gen, N)
    coef.indiv[ind:(ind+N-1), "relfit"] = curhist$Wi
    coef.indiv[ind:(ind+N-1), "emerge"] = curhist$emerge
    coef.indiv[ind:(ind+N-1), traitslist] = as.matrix(curhist[, traitslist])
    ind = ind+N
  }
  return(coef.indiv)
}

##NEXT DO THIS ONE!
# actTraitEff<-function(years.ind, years.list, pophistory, N){
#   #  Function for calculating the actual effect size of each coefficient for each indiv
#   #    This is done by finding the conditions when each individual emerged, and calculating the effect of each coefficient on that day.
#   #  Inputs:
#   traitslist = sprintf("b.%s", traits)
#   coef.indiv = matrix(data = 0, ncol = (3+length(traitslist)), nrow = N*numYears,
#                     dimnames = list(NULL, c("gen", traitslist, "relfit", "emerge")))
#   ind = 1
#   for(i.gen in 1:length(years.ind)){
#     curhist = pophistory[[i.gen]] #store the current year of population date
#     curyear = years.list[[years.ind[[i.gen]]]] #store the current year of envi conditions
#     coef.indiv[ind:(ind+N-1), "gen"] = rep(i.gen, N)
#     coef.indiv[ind:(ind+N-1), "relfit"] = curhist$Wi
#     coef.indiv[ind:(ind+N-1), "emerge"] = curhist$emerge
#     for(i.indiv in 1:N){
#       cur.econd = curyear[curhist$emerge[i.indiv], ] #grab the envi conditions of the day of emergence of current indiv
#       coef.indiv[ind, "b.day"] = cur.econd$day*curhist[i.indiv, "b.day"]
#       coef.indiv[ind, "b.temp"] = cur.econd$temp*curhist[i.indiv, "b.temp"]
#       coef.indiv[ind, "b.precip"] = cur.econd$precip*curhist[i.indiv, "b.precip"]
#       ind = ind+1
#     }
#   }
#   return(coef.indiv)
# }

actTraitEff<-function(years.ind, years.list, pophistory, N, traits){
  #  Function for calculating the actual effect size of each coefficient for each indiv
  #    This is done by finding the conditions when each individual emerged, and calculating the effect of each coefficient on that day.
  #  Inputs:
  numYears = length(pophistory)
  traitslist = sprintf("b.%s", traits)
  coef.indiv = matrix(data = 0, ncol = (3+length(traitslist)), nrow = N*numYears,
                    dimnames = list(NULL, c("gen", traitslist, "relfit", "emerge")))
  ind = 1
  for(i.gen in 1:length(years.ind)){
    curhist = pophistory[[i.gen]] #store the current year of population date
    curyear = years.list[[years.ind[[i.gen]]]] #store the current year of envi conditions
    coef.indiv[ind:(ind+N-1), "gen"] = rep(i.gen, N)
    coef.indiv[ind:(ind+N-1), "relfit"] = curhist$Wi
    coef.indiv[ind:(ind+N-1), "emerge"] = curhist$emerge
    for(i.indiv in 1:N){
      cur.econd = curyear[curhist$emerge[i.indiv], ] #grab the envi conditions of the day of emergence of current indiv
      coef.indiv[ind, traitslist] = as.matrix(cur.econd[, unlist(traits)]*100/curhist[i.indiv, traitslist])
      ind = ind+1
    }
  }
  return(coef.indiv)
}


#Does all the setup of population, runs simulation
#######################
# initializing population
#######################
##intialize a population of N individuals
# Their min and max values are determined by the start$ parameters
#Start by setting them all equal to zero. Fill in with the working traits from the traits variable
simrunner = function(N,traits,start,years.list,years.index,duration,sds,mutrate,numyears,fattail){
  b.day=b.temp=b.precip=b.cutemp=b.cuprecip=b.daysq=b.tempsq=b.precipsq=b.cutempsq=b.cuprecipsq=rep(0,N)
  for(i.trait in traits){
    if(start[[i.trait]][1]==0 & start[[i.trait]][2]==0){
      curvals=rep(0,N)
    }else{
      randnums=runif(n=N,min=start[[i.trait]][1],max=start[[i.trait]][2])
      randnums[randnums==0]=1/(10^10)
      curvals=randnums
    }
    curname=paste("b.",i.trait,sep="")
    assign(curname,curvals)
  }
  newpop<-data.frame(b.day,b.temp,b.precip,b.cutemp,b.cuprecip,b.daysq,b.tempsq,b.precipsq,b.cutempsq,b.cuprecipsq)
  pop<-selection(newpop,duration,year=years.list[[years.index[1]]],N,traits=traits)
  ###########################
  ## Running the Simulation #
  ###########################
  pophistory=runSim(startpop=pop,years.list=years.list,
                    years.ind=years.index,N=N,duration=duration,
                    sds=sds,mutrate=mutrate,generations=numYears,traits=traits, fattail=fattail)
  #Note: we've already used year 1 in initiating the pop
  #Based on the value of "runType", generate the appropriate type of data.
  return(pophistory)
}



####################################
# Year generation functions
####################################

yeargen<-function(dat.file, #the climate data file to use
                  fit.parms, #LIST of fitness parameters for the given fitness file
                  baseTemp = 0,
                  other.name = "moist",
                  decay = .2,
                  latitude = 38.5,
                  moist.norm = FALSE){
  #I have recently changed this!
  #We now use a single yeargen function for any climate data - the dat.file argument lets us specify the climate
  #And it uses "temp" and "other" as the two climate covariates to feed into the fitness function. Other might refer
  # to precip or moisture, but could be anything we want. Specified with the "other.name" variable
  # The decay parameter is used when creating the moisture coefficient
  #the datFile should refer to a pre-existing imputed data faile, like datFile = "davisDat.Rdata"
  require(geosphere)
  set_wrkdir()
  fileName = dat.file
  envdat = new.env()
  load(paste("data-years/", fileName, sep = ""), envir = envdat)
  years.list = envdat$yearlist
  years.ind = envdat$yearnames
  years.temp = do.call(rbind.data.frame, years.list)
  naminds = which(colnames(years.temp) %in% c("DAY.OF.YEAR", "YEAR", "PRCP", "TMAX"))
  colnames(years.temp)[naminds] = c("day", "year", "precip", "temp")
  years.temp = years.temp[, c("day", "year", "precip", "temp")]
  years.temp = years.temp[years.temp$day<366, ]#REMOVE LEAP DAY BY TAKING OFF LAST DAY OF YEAR FOR THOSE YEARS
  years.temp = cbind(years.temp,
                   moist = 0*years.temp$temp, #environmental moisture
                   photo = 0*years.temp$temp, #photoperiod
                   dprecip = 0*years.temp$temp,
                   dtemp = 0*years.temp$temp,
                   dmoist = 0*years.temp$temp,
                   cutemp = 0*years.temp$temp,
                   dcutemp = 0*years.temp$temp,
                   cuprecip = 0*(years.temp$precip),
                   daysq = 0*(years.temp$day),
                   tempsq = 0*(years.temp$temp),
                   precipsq = 0*(years.temp$precip),
                   cutempsq = 0*(cumsum(years.temp$temp)),
                   cuprecipsq = 0*(cumsum(years.temp$precip))
  )
  years.list = split(years.temp, f = years.temp$year)
  for(i.year in 1:length(years.list)){
    #calculating moisture
    rain = years.list[[i.year]]$precip
    moist.vec = rain*0
    moist = 0
    for(i in 1:length(rain)){
      moist = moist*exp(-decay)+rain[i]
      moist.vec[i] = moist
    }
    years.list[[i.year]]$moist = moist.vec;
    #calculating everything else
    years.list[[i.year]]$photo = daylength(latitude, 1:365)
    years.list[[i.year]]$dprecip = c(0, diff(years.list[[i.year]]$precip));
    years.list[[i.year]]$dtemp = c(0, diff(years.list[[i.year]]$temp));
    years.list[[i.year]]$dmoist = c(0, diff(years.list[[i.year]]$moist));
    years.list[[i.year]]$cutemp = cumsum(pmax(0, years.list[[i.year]]$temp-baseTemp));
    years.list[[i.year]]$dcutemp = c(0, diff(years.list[[i.year]]$cutemp));
    years.list[[i.year]]$cuprecip = cumsum(years.list[[i.year]]$precip);
    years.list[[i.year]]$daysq = (years.list[[i.year]]$day)^2;
    years.list[[i.year]]$tempsq = (years.list[[i.year]]$temp)^2;
    years.list[[i.year]]$precipsq = (years.list[[i.year]]$precip)^2;
    years.list[[i.year]]$cutempsq = cumsum((years.list[[i.year]]$temp)^2);
    years.list[[i.year]]$cuprecipsq = cumsum((years.list[[i.year]]$precip)^2);
  }
  yrsdf = do.call(rbind.data.frame, years.list)
  if(moist.norm != FALSE){
    yrsdf$moist = yrsdf$moist/max(yrsdf$moist)*moist.norm
  }
  fit.daily = fit_fn(years = yrsdf, other.name = other.name, fit.parms = fit.parms)
  yrsdf = cbind(yrsdf, fit.daily)
  years.list = split(yrsdf, yrsdf$year)
  return(list(years.list, years.ind))
}


yeargen.const<-function(numYears){
  #generate a sequence of years with identical, gaussian fitness curves, and constant envi conditions.
  #In this test, fitness is a gauss function centered on day 150
  modelYear = data.frame(day = 1:365, temp = rep(20, 365), precip = rep(.5, 365), fit.daily = dnorm(1:365, mean = 150, sd = 30))
  years.list = list(modelYear)
  # Each year data frame has $day, $precip, $tmean, $tmax, $tmin
  # This will be the same list for all configurations of years - this is essentially just our year database
  years.index = rep(1, numYears) # This is the list of which year.list data to use for each generation of the model
  return(list("years.list" =years.list, "years.index" = years.index))
}

yeargen.rand<-function(numYears){
  #generate a sequence of years with identical, gaussian fitness curves, and randomly fluctuating envi conditions.
  #In this test, fitness is a gauss function centered on day 150
  modelYear = data.frame(day = 1:365, temp = runif(n = 365, min = 0, max = 40), precip = rexp(n = 365, rate = 10), fit.daily = dnorm(1:365, mean = 150, sd = 30))
  years.list = list(modelYear)
  for(i in 2:max(numYears, 100)){
    years.list[[i]] = data.frame(day = 1:365, temp = runif(n = 365, min = 0, max = 40), precip = rexp(n = 365, rate = 10), fit.daily = dnorm(1:365, mean = 150, sd = 30))
  }
  # Each year data frame has $day, $precip, $tmean, $tmax, $tmin
  # This will be the same list for all configurations of years - this is essentially just our year database
  return(list("years.list" = years.list, "years.index" = years.index))
}

yeargen.davistest<-function(numYears, best.temp, sd.temp, best.precip, sd.precip){
  #Then we use the davis input data, and take some of the "good years" - ie no NANs - for our populations
  #input data
  davis.daily<-read.csv("davis-data/626713.csv", header = T, na.strings = "-9999")
  davis.daily$PRCP<-davis.daily$PRCP/10 #precips are reported in tenths of mm
  davis.daily$TMAX<-davis.daily$TMAX/10 #temps are reported in tenths of degree C
  davis.daily$TMIN<-davis.daily$TMIN/10 #temps are reported in tenths of degree C
  davis.daily$DATE2<-as.Date(as.character(davis.daily$DATE), format = "%Y %m %d") #DATE2 is date formatted
  davis.daily$JULIAN<-julian(davis.daily$DATE2, origin = as.Date("1892-12-31")) #1893-01-01 is day 1...
  davis.daily$YEAR<-as.numeric(substr(davis.daily$DATE, 1, 4)) #simple field for year
  davis.daily$MONTH<-as.numeric(substr(davis.daily$DATE, 5, 6)) #simple field for month
  davis.daily$DAY<-as.numeric(substr(davis.daily$DATE, 7, 8)) #simple field for day
  davis.daily<-davis.daily[, c("DATE2", "JULIAN", "YEAR", "MONTH", "DAY", "PRCP", "TMAX", "TMIN")] #simplified dataframe
  davis.yearlist<-split(davis.daily, davis.daily$YEAR) #list of each year separated
  #calculates the "day of year", i.e. Jan 1 is 1, and 12/31 is 365
  #adds a DAY.OF.YEAR column to each dataframe in the year list
  davis.yearnames<-unique(davis.daily$YEAR)
  for (i in 1:length(davis.yearnames)){
    davis.yearlist[[i]]$DAY.OF.YEAR<-julian(davis.yearlist[[i]]$DATE2, origin = as.Date(paste(davis.yearnames[i], "01", "01", sep = "-")))+1 #add +1 so that the first day of the year is 1, not zero.
  }
  yearlist.store = davis.yearlist
  goodyears = NULL
  for(iyear in davis.yearnames){
    nacount = sum(sum(is.na(davis.yearlist[[as.character(iyear)]])))
    daycount = dim(davis.yearlist[[as.character(iyear)]])[1]
    if(nacount==0 & daycount>364){goodyears = c(goodyears, iyear)}
  }
  davis.yearlist = davis.yearlist[as.character(goodyears)]
  davis.yearnames<-goodyears #gives a list of all the years in the data
  davis.daily<-unsplit(yearlist.store, davis.daily$YEAR) #using legacy "yearlist.store" to make unsplit happy
  # DAY.OF.YEAR = rep(0, dim(davis.daily)[1])
  # for(i in 1:length(DAY.OF.YEAR)){
  #   DAY.OF.YEAR[i] = sprintf("%02d%02d", davis.daily[i, "MONTH"], davis.daily[i, "DAY"])
  #
  # }
  # davis.daily = cbind(davis.daily, DAY.OF.YEAR)
  # davis.daily<-unsplit(davis.daily, davis.daily$YEAR)
  davis.daily.means<-aggregate(cbind(TMAX, TMIN, PRCP)~DAY.OF.YEAR, data = davis.daily[davis.daily$YEAR %in% goodyears, ], mean)
  davis.yearvar<-data.frame(row.names = davis.yearnames) #dataframe to hold environmental variability
  for (i in 1:length(davis.yearnames)){
    #temporary dataframe to compare with mean conditions
    #this creates a VAR.x for each year and a VAR.y for the daily means
    comparison<-merge(davis.yearlist[[i]], davis.daily.means, by = "DAY.OF.YEAR")
    #number of complete cases (is.na = F) for each year
    davis.yearvar[i, "TMAX.N"]<-sum(complete.cases(davis.yearlist[[i]]$TMAX))
    davis.yearvar[i, "TMIN.N"]<-sum(complete.cases(davis.yearlist[[i]]$TMIN))
    davis.yearvar[i, "PRCP.N"]<-sum(complete.cases(davis.yearlist[[i]]$PRCP))
    #sum of squared differences with an average year - how weird is each year?
    #some years have incomplete data, so this is the mean SS per observed day
    davis.yearvar[i, "TMAX.SS"]<-(sum(comparison$TMAX.x-comparison$TMAX.y, na.rm = T)^2)/davis.yearvar[i, "TMAX.N"]
    davis.yearvar[i, "TMIN.SS"]<-(sum(comparison$TMIN.x-comparison$TMIN.y, na.rm = T)^2)/davis.yearvar[i, "TMIN.N"]
    davis.yearvar[i, "PRCP.SS"]<-(sum(comparison$PRCP.x-comparison$PRCP.y, na.rm = T)^2)/davis.yearvar[i, "PRCP.N"]
    #CV within years - how variable is each year?
    davis.yearvar[i, "TMAX.CV"]<-sd(comparison$TMAX.x, na.rm = T)/mean(comparison$TMAX.x, na.rm = T)
    davis.yearvar[i, "TMIN.CV"]<-sd(comparison$TMIN.x, na.rm = T)/mean(comparison$TMIN.x, na.rm = T)
    davis.yearvar[i, "PRCP.CV"]<-sd(comparison$PRCP.x, na.rm = T)/mean(comparison$PRCP.x, na.rm = T)
    #sum of differences (not squared) with an average year - how hot/wet is each year?
    #some years have incomplete data, so this is the mean difference per observed day
    davis.yearvar[i, "TMAX.DEL"]<-sum(comparison$TMAX.x-comparison$TMAX.y, na.rm = T)/davis.yearvar[i, "TMAX.N"]
    davis.yearvar[i, "TMIN.DEL"]<-sum(comparison$TMIN.x-comparison$TMIN.y, na.rm = T)/davis.yearvar[i, "TMIN.N"]
    davis.yearvar[i, "PRCP.DEL"]<-sum(comparison$PRCP.x-comparison$PRCP.y, na.rm = T)/davis.yearvar[i, "PRCP.N"]
    ######################
    # Fitness generation #
    ######################
    #For now, daily incremental fitness will be found by multiplying two gaussian functions together:
    #  one for temp, that's maximized at best.temp with sd tempsd
    #  the other for precip that's maximized at best.precip with sd precipsd
    # We will then normalize the results to vary from 0 to 1
    years.list = davis.yearlist
    for(i.year in 1:length(years.list)){
      newyear = years.list[[i.year]]
      newyear = newyear[, c("DAY.OF.YEAR", "TMAX", "PRCP")]
      colnames(newyear)<-c("day", "tmax", "precip")
      daily.fit = dnorm(newyear$tmax, mean = best.temp, sd = sd.temp)*dnorm(newyear$precip, mean = best.precip, sd = sd.precip)
      daily.fit = (daily.fit-min(daily.fit))/(max(daily.fit)-min(daily.fit))
      years.list[[i.year]] = cbind(newyear, fit.daily = daily.fit)
    }
    # Each year data frame has $day, $precip, $tmean, $tmax, $tmin
    # This will be the same list for all configurations of years - this is essentially just our year database
  }
  years.index = rep(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), length.out = numYears) # This is the list of which year.list data to use for each generation of the model
  years.list = davis.yearlist #Replace this with code to grab a list of data frames. Each data frame is a year.
  return(list("years.index" = years.index, "years.list" = years.list))
}
######################
# Plotting Functions #
######################

emergePlot<-function(indivs, traitName){
  #Function for plotting trait values through time
  #Inputs:
  #  generations: vector of the generation of each individual to be plotted
  #  traivals: vector of the trait value of interest of each individ to be plotted
  #  mainlabel: label for the main graph
  #  ylabel: label for Y axis
  maxCount = 100 #maximum number of years to count
  generations = indivs[, "gen"]
  if(length(unique(generations))>maxCount){
    viewGens = floor(seq(min(generations), max(generations), length.out = maxCount))
    goodInd = generations %in% viewGens
    generations = generations[goodInd]
    indivs = indivs[goodInd, ]
  }
  plot(jitter(generations), indivs[, traitName], type = 'n',
       main = paste("Actual effect size of", traitName),
       xlab = "Generation",
       ylab = paste(traitName, "effect size"),
       cex.lab = 1.4, cex.main = 1.4)
  #Plot the "emerge before last day" indivs
  points(jitter(generations[indivs[, "emerge"]>364]), indivs[indivs[, "emerge"]>364, traitName], pch = 3, col = 'blue')
  points(jitter(generations[indivs[, "emerge"]<365]), indivs[indivs[, "emerge"]<365, traitName], pch = 1)
}
emergePlotYlim<-function(indivs, traitName, ylim){
  #Function for plotting trait values through time
  #Inputs:
  #  generations: vector of the generation of each individual to be plotted
  #  traivals: vector of the trait value of interest of each individ to be plotted
  #  mainlabel: label for the main graph
  #  ylabel: label for Y axis
  maxCount = 100 #maximum number of years to count
  generations = indivs[, "gen"]
  if(length(unique(generations))>maxCount){
    viewGens = floor(seq(min(generations), max(generations), length.out = maxCount))
    goodInd = generations %in% viewGens
    generations = generations[goodInd]
    indivs = indivs[goodInd, ]
  }
  plot(jitter(generations), indivs[, traitName], type = 'n',
       main = paste("Actual effect size of", traitName),
       xlab = "Generation",
       ylab = paste(traitName, "effect size"),
       cex.lab = 1.4, cex.main = 1.4,
       ylim = ylim)
  #Plot the "emerge before last day" indivs
  points(jitter(generations[indivs[, "emerge"]>364]), indivs[indivs[, "emerge"]>364, traitName], pch = 3, col = 'blue')
  points(jitter(generations[indivs[, "emerge"]<365]), indivs[indivs[, "emerge"]<365, traitName], pch = 1)
}
traiteffplot<-function(indivs, traitName){
  #Function for plotting trait values through time
  #Inputs:
  #  generations: vector of the generation of each individual to be plotted
  #  traivals: vector of the trait value of interest of each individ to be plotted
  #  mainlabel: label for the main graph
  maxCount = 100 #maximum number of years to count
  generations = indivs[, "gen"]
  if(length(unique(generations))>maxCount){
    viewGens = floor(seq(min(generations), max(generations), length.out = maxCount))
    goodInd = generations %in% viewGens
    generations = generations[goodInd]
    indivs = indivs[goodInd, ]
  } #  ylabel: label for Y axis
  plot(jitter(indivs[, "gen"]), indivs[, traitName], pch = 1,
       main = paste("Expected effect size of", traitName),
       xlab = "Generation",
       ylab = paste(traitName, "effect size"),
       cex.lab = 1.4, cex.main = 1.4)
}
traitplot<-function(indivs, traitName){
  #Function for plotting trait values through time
  #Inputs:
  #  generations: vector of the generation of each individual to be plotted
  #  traivals: vector of the trait value of interest of each individ to be plotted
  #  mainlabel: label for the main graph
  #  ylabel: label for Y axis
  maxCount = 100 #maximum number of years to count
  generations = indivs[, "gen"]
  if(length(unique(generations))>maxCount){
    viewGens = floor(seq(min(generations), max(generations), length.out = maxCount))
    goodInd = generations %in% viewGens
    generations = generations[goodInd]
    indivs = indivs[goodInd, ]
  }
  plot(jitter(generations), indivs[, traitName], type = 'n',
       main = paste("Coefficient values of", traitName),
       xlab = "Generation",
       ylab = paste(traitName, "values"),
       cex.lab = 1.4, cex.main = 1.4)
  #Plot the "emerge before last day" indivs
  points(jitter(generations[indivs[, "emerge"]>364]), indivs[indivs[, "emerge"]>364, traitName], pch = 3, col = 'blue')
  points(jitter(generations[indivs[, "emerge"]<365]), indivs[indivs[, "emerge"]<365, traitName], pch = 1)
}


######################################################
# Function versions of windows_plot and windows_save #
######################################################
windows_plot = function(list.all){
  allres=within(list.all, {
    setwd("results")
    resultsdir=sprintf("%s/resRun%s",runsname,runName)
    setwd(resultsdir)

    #Make a matrix of daily fitnesses for all years, cutting out the last day of leap years
    yearFit=matrix(0,nrow=length(years.index),ncol=365)
    count=1
    for(i in years.index){
      curfits=years.list[[i]]$fit.daily
      if(length(curfits)==366){curfits=curfits[-366]} #to handle leap years, remove last day
      yearFit[count,]=curfits;
      count=count+1
    }
    #Calculate the mean fitness accrued each day
    par(mar=c(5,5,4,3))
    meanFit=apply(yearFit,2,mean)
    #calculate the mean fitness for emerging on day x (for all days) [this is using mean fitness accrued per day]
    meanFitSum=rollapply(c(meanFit,rep(0,duration-1)),duration,by=1,sum)
    for(i.day in 1:365){
      meanFitSum=c(meanFitSum,sum(rep(meanFit)[i.day:(i.day+duration-1)]))
    }

    x11(width=9,height=6)
    if(plotExtra==TRUE){
      for(curgen in c(1:20,seq(21,length(years.index),length=10))){
        curgen=round(curgen)
        #arheight=rep(max(meanFit)*1.1,N) #upper bound used for making plots look good
        emergeDay=pophistory[[curgen]]$emerge
        #     plot(meanFit,type='l',ylim=c(0,max(meanFit)*1.2))
        #     arrows(y0=jitter(arheight,factor=1.5),x0=emergeDay,x1=emergeDay+duration-1,length=.1)
        #     dev.print(pdf,paste("dailyfit-run",runName,"-gen",curgen,"-meanfit.pdf",sep=""))
        #     plot(meanFitSum,type='l',ylim=c(0,max(meanFitSum)*1.2),
        #          main=paste("Mean fitness gained, gen",curgen),
        #          ylab="Fitness gained",
        #          xlab="Julian date",
        #          cex.lab=1.3,
        #          cex.main=1.3)
        #     arheight=jitter(rep(max(meanFitSum)*1.05,N),factor=.8)
        #     arrows(y0=arheight+.05*max(meanFitSum),x0=emergeDay,y1=arheight,length=.1)
        #     dev.print(pdf,paste("dailyfitSum-run",runName,"-gen",curgen,"-meanfit.pdf",sep=""))
        #now calculate the fitSum for THIS YEAR ONLY
        FitSum=rollapply(c(years.list[[years.index[[curgen]]]]$fit.daily,rep(0,duration-1)),duration,by=1,sum)
        plot(FitSum,type='l',ylim=c(0,max(FitSum)*1.2),
             main=paste("Fitness gained this year, gen",curgen),
             ylab="Fitness gained",
             xlab="Julian date",
             cex.lab=1.3,
             cex.main=1.3)
        arheight=jitter(rep(max(FitSum)*1.05,N),factor=.8)
        arrows(y0=arheight+.05*max(FitSum),x0=emergeDay+lag,y1=arheight,length=.1)
        points(meanFitSum,type='l',lty=3,lwd=2)
        dev.print(pdf,paste("dailyfitSum-run",runName,"-gen",curgen,"-actualfit.pdf",sep=""))
        #Now plot each of the coefs by emergence day.
        #     for(coefName in c("b.day","b.temp","b.precip")){
        #       plot(pophistory[[curgen]][,coefName],pophistory[[curgen]][,"emerge"],
        #            main=paste(coefName, "by emergence, gen", curgen),
        #            ylab="emergence",
        #            xlab=coefName,
        #            cex.lab=1.3,
        #            cex.main=1.3)
        #       dev.print(pdf,paste("Coef_x_emerge-",coefName,"-run",runName,"-gen",curgen,".pdf",sep=""))
        #     }
        #     coefList=c("b.day","b.temp","b.precip")
        #     for(i in 1:3){
        #       coef1=coefList[i]
        #       coef2=coefList[(i %% 3)+1]
        #       curpop=pophistory[[curgen]]
        #       plot(x=curpop[,coef1],y=curpop[,coef2],type='n',
        #            main=paste(coef1, "by", coef2, "gen", curgen),
        #            ylab=coef2,
        #            xlab=coef1,
        #            cex.lab=1.3,
        #            cex.main=1.3)
        #       points(x=(curpop[curpop[,"emerge"]>364,coef1]),y=(curpop[curpop[,"emerge"]>364,coef2]),pch=3,col='blue')
        #       points(x=(curpop[curpop[,"emerge"]<365,coef1]),y=(curpop[curpop[,"emerge"]<365,coef2]),pch=1,col='black')
        #       dev.print(pdf,paste("Coef_x_coef-",coef1,"x",coef2,"-run",runName,"-gen",curgen,".pdf",sep=""))
        #     }
      }
    }
    #Calculating changes in mean fitness through time
    maxfit=maxActfit=meanfit=emerge.ideal=quant50=rep(0,length(years.index))
    emerge=matrix(0,ncol=N,nrow=length(years.index))
    for(curgen in 1:numYears){
      meanfit[curgen]=mean(pophistory[[curgen]]$Wi)
      maxActfit[curgen]=max(pophistory[[curgen]]$Wi)
      cur.fitness=years.list[[years.index[curgen]]]$fit.daily
      cur.fitness.durated=rollapply(c(cur.fitness,rep(0,duration-1)),duration,by=1,sum)
      maxfit[curgen]=max(cur.fitness.durated)
      emerge[curgen,]=pophistory[[curgen]]$emerge
      emerge.ideal[curgen]=min(which(cur.fitness.durated==maxfit[curgen]))-lag
      quant50[curgen]=which(cumsum(cur.fitness.durated)/sum(cur.fitness.durated)>=.5)[1]
    }
    #plot emergence times
    maxCount=100 #maximum number of years to count
    generations=1:length(years.index)
    viewGens=generations
    if(length(generations)>maxCount){
      viewGens=floor(seq(min(generations),max(generations),length.out=maxCount))
    }
    #Plot when organisms emerge
    matplot(jitter(viewGens),emerge[viewGens,]+lag,type='p',pch=1,col='black',
            main=paste("Emergence days"),
            xlab="Generation",
            ylab="Emergence day",
            cex.lab=1.4,cex.main=1.4)
    #add `optimal emergence day' - note this is a vast oversimplification
    points(viewGens,emerge.ideal[viewGens],col="red",pch=4,lwd=2)
    points(viewGens,quant50[viewGens],col="blue",pch=4,lwd=2)
    dev.print(pdf,paste("emerge-run",runName,"-gen",curgen,".pdf",sep=""))

    #plot mean fitness through time, showing max possible fitness
    plot(maxfit,type='l',col='red',
         main=paste("Mean fitness through time for run",runName),
         xlab="generation",
         ylab="Raw mean fitness",
         sub="red is maximum possible",
         ylim=c(0,max(maxfit))
    )
    points(1:length(meanfit),meanfit,type='l')
    dev.print(pdf,paste("meanfit-run",runName,"-gen",curgen,".pdf",sep=""))
    #plot mean fitness through time, normalized by max fitness
    plot(meanfit/maxfit,type='l',col='black',
         main=paste("Mean fitness / max possible",runName),
         ylab="normalized mean fitness",
         xlab="generation",
         sub="red is maximum possible",
         ylim=c(0,1)
    )
    abline(h=1,col='red')
    dev.print(pdf,paste("meanfitNorm-run",runName,"-gen",curgen,".pdf",sep=""))
    #plot max potential fitness and max actual fitness through time
    #  ie the fittest individual of each generatoin
    plot(maxfit,type='l',col='red',
         main=paste("Max achieved fitness through time for run",runName),
         xlab="generation",
         ylab="Raw max fitness",
         sub="red is maximum possible",
         ylim=c(0,max(maxfit))
    )
    points(1:length(maxActfit),maxActfit,type='l')
    dev.print(pdf,paste("maxfit-run",runName,"-gen",curgen,".pdf",sep=""))

    #Looking at coef changes through time
    #  The act.eff is the actual effect size, found by multiplying the coefficient of each indiv by the environmental conditions of their day of emergence.
    #    Those plots use crosses to represent individuals who didn't emerge until the final day, and circles for those that emerged on a normal day (ie their cue
    act.eff=actTraitEff(years.index,years.list,pophistory,N,traits)
    act.vals=actTraitVals(pophistory,numYears,N)
    traitslist=sprintf("b.%s",traits)

    #x11(width=9,height=6)
    par(mar=c(5,5,4,4))
    #Plot it all in one
    par(mfrow=c(3,1))
    ind=1
    while(ind<=length(traitslist)){
      plotlist=ind:min(ind+2,length(traitslist))
      for(cur.trait in plotlist){
        emergePlot(indivs=act.eff,trait=traitslist[cur.trait])
      }
      dev.print(pdf,paste("coefEffects-",paste(traitslist[plotlist],collapse='-'),"-actual-run",runName,".pdf",sep=""))
      for(cur.trait in plotlist){
        emergePlotYlim(indivs=act.eff,trait=traitslist[cur.trait],ylim=c(0,100))
      }
      dev.print(pdf,paste("coefEffects-",paste(traitslist[plotlist],collapse='-'),"-actual-run-ylim",runName,".pdf",sep=""))
      for(cur.trait in plotlist){
        emergePlot(indivs=act.eff[act.eff[,"gen"]>burnIn,],trait=traitslist[cur.trait])
      }
      dev.print(pdf,paste("coefEffects-",paste(traitslist[plotlist],collapse='-'),"-actual-postburn-run",runName,".pdf",sep=""))
      for(cur.trait in plotlist){
        traitplot(indivs=act.vals,trait=traitslist[cur.trait])
      }
      dev.print(pdf,paste("coefVals-",paste(traitslist[plotlist],collapse='-'),"-run",runName,".pdf",sep=""))
      for(cur.trait in plotlist){
        traitplot(indivs=act.vals[act.vals[,"gen"]>burnIn,],trait=traitslist[cur.trait])
      }
      dev.print(pdf,paste("coefVals-",paste(traitslist[plotlist],collapse='-'),"-postburn-run",runName,".pdf",sep=""))
      ind=max(plotlist)+1
    }
    #Plot phenotypes through time

    par(mfrow=c(1,1))
    pop.temp=do.call(rbind.data.frame,pophistory[c(seq(2,length(years.index),by=100))])
    traitmins=NULL
    for(i.trait in traitslist){traitmins=c(traitmins,min(pop.temp[,i.trait]))}
    traitmaxs=NULL
    for(i.trait in traitslist){traitmaxs=c(traitmaxs,max(pop.temp[,i.trait]))}
    #x11(width=9,height=6)
    if(plotPheno==TRUE){
      traitslist=sprintf("b.%s",traits)
      if(length(traitslist)==1){ #only one trait, do histogram
        for(curgen in c(1,seq(2,numYears,length=20))){
          curgen=round(curgen)
          cur.pop=pophistory[[curgen]]
          hist(cur.pop[,traitslist])
          dev.print(pdf,sprintf("pheno-gen%06d-run%s.pdf",curgen,runName))
        }
      }else if(length(traitslist)==2){ #two traits, do 2d plot
        for(curgen in c(1,seq(2,numYears,length=20))){
          curgen=round(curgen)
          cur.pop=pophistory[[curgen]]
          plot(cur.pop[,traitslist])
          dev.print(pdf,sprintf("pheno-gen%06d-run%s.pdf",curgen,runName))
        }
      }else if(length(traitslist)==3){ #do 3d plot
        for(curgen in c(1,seq(2,numYears,length=20))){
          curgen=round(curgen)
          cur.pop=pophistory[[curgen]]
          scatterplot3d(jitter(as.matrix(cur.pop[,traitslist])),type='h',
                        xlim=c(traitmins[1],traitmaxs[1]),
                        ylim=c(traitmins[2],traitmaxs[2]),
                        zlim=c(traitmins[3],traitmaxs[3]))
          dev.print(pdf,sprintf("pheno-gen%06d-run%s.pdf",curgen,runName))
        }
      }else{ #more than 3 traits, do ndms
        for(curgen in c(1,seq(2,numYears,length=20))){
          curgen=round(curgen)
          cur.pop=pophistory[[curgen]]
          cur.ndms=metaMDS(cur.pop[,traitslist],k=2,trymax=100,autotransform = TRUE)
          plot(cur.ndms,main=paste("NDMS for generation", curgen))
          dev.print(pdf,sprintf("pheno-gen%06d-run%s.pdf",curgen,runName))
        }
      }
    }


    # a few new plotting ideas ------------------------------------------------

    library(reshape2)
    library(ggplot2)
    library(gridExtra)

    #coeff.eff.sum<-aggregate(cbind(b.day,b.cutemp,b.cuprecip)~gen,data=act.eff,mean)

    coeff.eff.sum<-aggregate(act.eff[,traitslist]~gen,data=act.eff,mean)
    # coeff.vals.sum<-aggregate(cbind(b.day,b.temp,b.precip)~gen,data=act.vals,mean)
    #to look at a smaller subset of the data, here up to 50 generations
    #coeff.eff.sum<-coeff.eff.sum[coeff.eff.sum$gen<50,]

    coeff.eff.sum.melt<- melt(coeff.eff.sum, id.var="gen")

    p1<-ggplot(coeff.eff.sum.melt,aes(x=gen,y=value,fill=variable))+geom_smooth()

    p2<-ggplot(coeff.eff.sum.melt,aes(x=gen,y=value,color=variable))+geom_line()+geom_point()+geom_smooth()

    grid.arrange(p1,p2,ncol=1)
    dev.print(pdf,sprintf("ggtraits.eff-run%s.pdf",curgen,runName))
    })
  return(list(act.eff=allres$act.eff,
              meanfit=allres$meanfit,
              maxfit=allres$maxfit,
              runName=allres$runName)) #for later use - need it to make summary info
}


windows_save = function(list.all){
  with(list.all, {
    setwd("results")
    resultsdir=sprintf("%s/resRun%s",runsname,runName)
    unlink(resultsdir,recursive = TRUE)
    dir.create(resultsdir,showWarnings = FALSE)
    setwd(resultsdir)
    save(list=ls(all.names=TRUE),file="dat.RData")
  })
}

