#Run parameters
fitfile="standgauss" #Fitness function to use. Currently only one defined is standgauss
locName="davis" #Which location's data to use. CUrrently only tested on davis
traits=traitslist=c("day","temp","cutemp") #which traits to compare
runnum="test4" #Name for saving the results

#Climate parameters
yearvarMax=2500 #max year variation to test
dayvarMax=25 #max day variation to test
numpts=10 #number of year and day variation levels to test.
numYears=100 #number of years of climate to generate
actprecip=FALSE #If false, set all precips to zero. If true, use the smoothed precip values

#Organism parameters
lag=1 #number of days of lag between orgs deciding to emerge and them collecting fitness.
duration=10 #number of days organisms are emerged and collecting fitness

#Optimization
pointcheck=1000 #total number of points to check in initial pass
fastnum=20 #number of the best points from pointcheck pass to do a fast optimization on
slownum=10 #number of the best points from fastnum pass to do a slow optimization on
