obtainBestInitial <- function(x, y, method){	
	source("multivariate-regression.R")
	# Obtain best initial variable for feature selection process.
	# method can be Hotelling (for Hotelling-Lawley trace) or random.
	# x is matrix/dataframe of 'predictors' or variables that vary per class
	# y is matrix/vector of 'response' or class variable(s)
	# returns column number of x of 'best' initial explanatory variable for class
	maxGain = 0.0
	maxVar    = NULL
	if(method=="Hotelling"){
		for(i in 1:NCOL(x)){
			currentTest = mvar(x[,i],y)
			if(currentTest$Hotelling > maxGain){			
				maxGain   = currentTest$Hotelling
				maxVar    = i
			}
		}
	}
	if(method=="random"){
		# Just pick a random column (returns Hotelling-Lawley trace stat)
		maxVar = sample(1:NCOL(x),1)
		currentTest = mvar(x[,maxVar],y)
	         maxGain = currentTest$Hotelling		
	}

	return(list(maxVar=maxVar,maxGain=maxGain))	
}
