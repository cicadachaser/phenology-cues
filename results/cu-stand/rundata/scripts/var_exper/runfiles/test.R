#Run parameters
fitfile="standgauss" #Fitness function to use. Currently only one defined is standgauss
locName="davis" #Which location's data to use. CUrrently only tested on davis
traits=traitslist=c("day","temp","cutemp") #which traits to compare
runnum="test" #Name for saving the results

#Climate parameters
yearstdMax=14 #max year variation to test (SD of shifts in days)
daystdMax=7 #max day temp variation to test
numpts=6 #number of year and day variation levels to test (computationally intensive with square).
numYears=300 #number of years of climate to generate 
actprecip=FALSE #If false, set all precips to zero. If true, use the smoothed precip values

#Organism parameters
lag=2 #number of days of lag between orgs deciding to emerge and them collecting fitness.
duration=10 #number of days organisms are emerged and collecting fitness
baseTempQ=30 #quantile of temperature to use for "base temp". Applied to the MEAN YEAR temperatures.

#Optimization
pointcheck=1000 #total number of points to check in initial pass
fastnum=20 #number of the best points from pointcheck pass to do a fast optimization on
slownum=10 #number of the best points from fastnum pass to do a slow optimization on
