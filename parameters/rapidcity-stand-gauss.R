#########################
# Simulation parameters #
#########################
runType="ithacaDat.Rdata" ##THIS DETERMINES WHAT KIND OF YEARS WE'RE USING!
#unitTestConst is for running the population through a unit test with the same gaussian fitness every year
#and constant environmental conditions
#unitTestRand will be for running the populations through a
#unit test with the same gaussian fitness every year and random envi conditions
#standard is for running the populations through a set of replications of the first 10 good years of the davis data
traits=c("day","temp","moist","cutemp","daysq","tempsq","moistsq","cutempsq")
# traits=c("day","temp","moist","cutemp")
numsims=10 # number of simulations of each type to do
#traits=c("day","temp","precip","cutemp","cuprecip","daysq","tempsq","precipsq","cutempsq","cuprecipsq")
plotExtra=TRUE # do we plot snapshots of emergence through time?
plotPheno=FALSE # do we plot snapshots of phenotype through time?
viewLength=500 #for comparisons of simulation types
#  how many generations (starting from the final and working backwards) to plot/compare
duration=10 #number of days organizm is emerged.
lag=1 #number of days between decision and obtaining fitness
N=50 #number of individuals
numYears=1000 #number of years to simulate
burnIn=100 #number of years to not plot (to avoid scale issues from broad initial population traits)
mutdist=.005 #What fraction of the total "cue space" should mutations (on average) traverse (kinda).
yearLabel="A" #For deciding which set of year randomization to use
fitshape="standgauss" #which fitness shape to use. use the name of file in fitcurve folder, minuse the .R
yearSet="all.years" #which set of years do we want to use?
baseTemp=0 #cutoff when calculating degree days/cumulative temperature
other.name="moist" #can change this to base fitness of temp and something other than moisture
decay=.1 #decay parameter for moisture calculation
fattail=TRUE#decision on whether or not to use fattail distribution (cauchy) for generating mutation distances
latitude=38.5#for use in photoperiod - chose value near davis
moist.norm=FALSE
##NOTE: specific fitness parameters can be found in the fitness scripts in fitcurve folder.
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
