maxcues<-list(#these are approximately the max values of the first five years of the Davis data.
  # use them for scaling the starting values, and for
  day=365,
  temp=45,
  precip=100,
  moist=100,
  cutemp=9000,
  cuprecip=600,
  daysq=133225,
  tempsq=1971.36,
  precipsq=10000,
  moistsq=10000,
  cutempsq=250000,
  cuprecipsq=20000
)
start<-list(#these are used to generate the starting values of individuals. Starting values will produce
  # individuals who, if they relied solely on a single cue, would emerge uniformly across a range of
  # cues from 1 to x times approximately the max value, where x is the number of traits used.
  # The max value is selected to be approximately the maximum found in the first five years of the davis
  # data set, and the min is set to approximately the minimum found. Note that we can't use zero.
  day=c(1,maxcues$day*(length(traits))),
  temp=c(5,maxcues$temp*(length(traits))),
  precip=c(.01,maxcues$precip*(length(traits))),
  moist=c(.01,maxcues$precip*(length(traits))),
  cutemp=c(5,maxcues$cutemp*(length(traits))),
#   day=c(1,maxcues$day),
#   temp=c(5,maxcues$temp),
#   precip=c(.01,maxcues$precip),
#   moist=c(.01,maxcues$precip),
#   cutemp=c(5,maxcues$cutemp),
  cuprecip=c(.01,maxcues$cuprecip*sqrt(length(traits))),
  daysq=c(-100*maxcues$daysq*(length(traits)),-99*maxcues$daysq*(length(traits))),
  tempsq=c(-100*maxcues$tempsq*(length(traits)),-99*maxcues$tempsq*(length(traits))),
  precipsq=c(-100*maxcues$precipsq*(length(traits)),-99*maxcues$precipsq*(length(traits))),
  moistsq=c(-100*maxcues$moistsq*(length(traits)),-99*maxcues$moistsq*(length(traits))),
  cutempsq=c(-100*maxcues$cutempsq*(length(traits)),-99*maxcues$cutempsq*(length(traits))),
  cuprecipsq=c(.01,maxcues$cuprecipsq*(length(traits)))
)
temporary<-list( #temporary object for masking unused traits
  day=0,
  temp=0,
  moist=0,
  precip=0,
  cutemp=0,
  cuprecip=0,
  daysq=0,
  tempsq=0,
  moistsq=0,
  precipsq=0,
  cutempsq=0,
  cuprecipsq=0
)
for(i.trait in traits){
  temporary[i.trait]=start[i.trait]
}
start=temporary
####THE SDS BELOW NEED TO BE RESCALED FOR FAIRNESS.
sds<-list( #standard deviations for trait mutations.
  day=mutdist*maxcues$day,
  temp=mutdist*maxcues$temp,
  moist=mutdist*maxcues$moist,
  precip=mutdist*maxcues$precip,
  cutemp=mutdist*maxcues$cutemp,
  cuprecip=mutdist*maxcues$cuprecip, #max
  daysq=mutdist*maxcues$daysq,
  tempsq=mutdist*maxcues$tempsq,
  moistsq=mutdist*maxcues$moistsq,
  precipsq=mutdist*maxcues$precipsq,
  cutempsq=mutdist*maxcues$cutempsq,
  cuprecipsq=mutdist*maxcues$cuprecipsq)
mutrate<-list( #probability of each trait mutating in an individual. Mutations are independent of one another
  day=1,
  temp=1,
  moist=1,
  precip=1,
  cutemp=1,
  cuprecip=1,
  daysq=1,
  tempsq=1,
  moistsq=1,
  precipsq=1,
  cutempsq=1,
  cuprecipsq=1)
#Now ensure that mutrate is zero for any trait that ISN'T in the trait list
temporary<-list( #probability of each trait mutating in an individual. Mutations are independent of one another
  day=0,
  temp=0,
  precip=0,
  moist=0,
  cutemp=0,
  cuprecip=0,
  daysq=0,
  tempsq=0,
  precipsq=0,
  moistsq=0,
  cutempsq=0,
  cuprecipsq=0
)
for(i.trait in traits){
  temporary[i.trait]=mutrate[i.trait]
}
mutrate=temporary