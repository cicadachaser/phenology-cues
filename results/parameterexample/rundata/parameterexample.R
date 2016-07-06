#########################
# Simulation parameters #
#########################
runType="standard" ##THIS DETERMINES WHAT KIND OF YEARS WE'RE USING!
#unitTestConst is for running the population through a unit test with the same gaussian fitness every year
#and constant environmental conditions
#unitTestRand will be for running the populations through a
#unit test with the same gaussian fitness every year and random envi conditions
#standard is for running the populations through a set of replications of the first 10 good years of the davis data
traits=c("day","temp","precip")
numsims=5 # number of simulations of each type to do
#traits=c("day","temp","precip","cutemp","cuprecip","daysq","tempsq","precipsq","cutempsq","cuprecipsq")
plotExtra=TRUE # do we plot snapshots of emergence through time?
plotPheno=FALSE # do we plot snapshots of phenotype through time?
viewLength=500 #for comparisons of simulation types,
#  how many generations (starting from the final and working backwards) to plot/compare
duration=10 #number of days organizm is emerged.
N=100 #number of individuals
numYears=500 #number of years to simulate
burnIn=200 #number of years to not plot (to avoid scale issues from broad initial population traits)
mutdist=.01 #What fraction of the total "cue space" should mutations (on average) traverse (kinda).
yearLabel="A" #For deciding which set of year randomization to use
fitshape="standgauss"


