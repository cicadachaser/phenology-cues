maxcues<-list(#these are approximately the max values of the first five years of the Davis data.
  # use them for scaling the starting values, and for
  day=365,
  temp=45,
  precip=100,
  cutemp=9000,
  cuprecip=600,
  daysq=133225,
  tempsq=1971.36,
  precipsq=10000,
  cutempsq=250000,
  cuprecipsq=20000
)
start<-list(#these are used to generate the starting values of individuals. Starting values will produce
  # individuals who, if they relied solely on a single cue, would emerge uniformly across a range of
  # cues from 1 to x times approximately the max value, where x is the number of traits used.
  # The max value is selected to be approximately the maximum found in the first five years of the davis
  # data set, and the min is set to approximately the minimum found. Note that we can't use zero.
  day=c(1,maxcues$day*length(traits)),
  temp=c(5,maxcues$temp*length(traits)),
  precip=c(.01,maxcues$precip*length(traits)),
  cutemp=c(5,maxcues$cutemp*length(traits)),
  cuprecip=c(.01,maxcues$cuprecip*length(traits)),
  daysq=c(1,maxcues$daysq*length(traits)),
  tempsq=c(5,maxcues$tempsq*length(traits)),
  precipsq=c(.01,maxcues$precipsq*length(traits)),
  cutempsq=c(5,maxcues$cutempsq*length(traits)),
  cuprecipsq=c(.01,maxcues$cuprecipsq*length(traits))
)
temporary<-list( #temporary object for masking unused traits
  day=0,
  temp=0,
  precip=0,
  cutemp=0,
  cuprecip=0,
  daysq=0,
  tempsq=0,
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
  precip=mutdist*maxcues$precip,
  cutemp=mutdist*maxcues$cutemp,
  cuprecip=mutdist*maxcues$cuprecip, #max
  daysq=mutdist*maxcues$daysq,
  tempsq=mutdist*maxcues$tempsq,
  precipsq=mutdist*maxcues$precipsq,
  cutempsq=mutdist*maxcues$cutempsq,
  cuprecipsq=mutdist*maxcues$cuprecipsq)
mutrate<-list( #probability of each trait mutating in an individual. Mutations are independent of one another
  day=.01,
  temp=.01,
  precip=.01,
  cutemp=.01,
  cuprecip=.01,
  daysq=.01,
  tempsq=.01,
  precipsq=.01,
  cutempsq=.01,
  cuprecipsq=.01)
#Now ensure that mutrate is zero for any trait that ISN'T in the trait list
temporary<-list( #probability of each trait mutating in an individual. Mutations are independent of one another
  day=0,
  temp=0,
  precip=0,
  cutemp=0,
  cuprecip=0,
  daysq=0,
  tempsq=0,
  precipsq=0,
  cutempsq=0,
  cuprecipsq=0
)
for(i.trait in traits){
  temporary[i.trait]=mutrate[i.trait]
}
mutrate=temporary