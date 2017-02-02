#########################
# Simulation parameters #
#########################
runType="standard" ##THIS DETERMINES WHAT KIND OF YEARS WE'RE USING!
#unitTestConst is for running the population through a unit test with the same gaussian fitness every year
#and constant environmental conditions
#unitTestRand will be for running the populations through a
#unit test with the same gaussian fitness every year and random envi conditions
#standard is for running the populations through a set of replications
traits=c("day","temp","precip")
#traits=c("day","temp","precip","cutemp","cuprecip","daysq","tempsq","precipsq","cutempsq","cuprecipsq")
numsims=10 # number of simulations of each type to do
plotExtra=TRUE # do we plot snapshots of emergence through time?
plotPheno=TRUE # do we plot snapshots of phenotype through time? 3-D scatterplot
viewLength=500 #for comparisons of simulation types,
#  how many generations (starting from the final and working backwards) to plot/compare
duration=10 #number of days organizm is emerged.
lag=1 #number of days between decision and obtaining fitness
N=1000 #number of individuals
numYears=500 #number of years to simulate
burnIn=200 #number of years to not plot (to avoid scale issues from broad initial population traits)
mutdist=.01 #What fraction of the total "cue space" should mutations (on average) traverse (kinda).
yearLabel="A" #For deciding which set of year randomization to use
fitshape="standgauss"#shape of the fitness function
yearSet="all.years"
#######################
# For analytic script #
#######################
pointcheck=10000 #number of points to evaluate initially
fastnum=20 #number of points to test quickly
slownum=10 #number of points to test slowly

################################
# For analytic script plotting #
################################
pointdense=100
num.plots=10
plot.traits=c("day","cutemp") #Traits to plot on countours. Must be 2 traits from the 'traits' vector
fix.traits=setdiff(traits, plot.traits)
