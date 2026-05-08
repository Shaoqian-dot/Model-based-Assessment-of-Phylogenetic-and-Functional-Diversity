#Create a function for parameter truncation
trun <- function(parameters, K){#PT: parameter truncation
  #  K=5 # try different values (bigger is better as long as no NaN!)
  whichBig=which(parameters>K)
  whichSmall=which(parameters< -K)
  newPar=parameters
  newPar[whichBig] = K #truncating parameters at K
  newPar[whichSmall] = -K #truncating parameters so no smaller than -K
  names(newPar)=names(parameters)
  return(list(par=newPar,which=c(whichBig,whichSmall)))
}
#Create a function for parameter truncation