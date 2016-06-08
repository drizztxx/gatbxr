## SUS.R          (Stochastic Universal Sampling)
##
## This function performs selection with STOCHASTIC UNIVERSAL SAMPLING.
##
## Syntax:  NewChrIx = sus(FitnV, Nsel)
##
## Input parameters:
##    FitnV     - Vector containing the fitness values of the
##                individuals in the population.
##    Nsel      - number of individuals to be selected
##
## Output:
##               Vector containing the indexes of the selected
##                individuals relative to the original population, shuffled.
##                The new population, ready for mating, can be obtained
##                by calculating OldChrom[NewChrIx,].
##
## Author:     Hartmut Pohlheim (Carlos Fonseca)
##             David Zhao (Modified for R)
##
## Date: 12May2016

sus <- function(FitnV,
                Nsel){
  ## Identify the population size (Nind)
  Nind <- NROW(FitnV)
  
  ## Perform stochastic universal sampling
  cumfit <- cumsum(FitnV)
  trials <- cumfit[Nind]/Nsel * (runif(1) + 0:(Nsel-1))
  NewChrIx <- rep(0,Nsel)
  j = 1
  for(i in 1:Nind){
    repeat{
      if (j > Nsel) break
      if (trials[j] < cumfit[i]) {
        NewChrIx[j] = i
        j = j + 1
      } else break
    }
  }
  ## Shuffle new population
  NewChrIx <- NewChrIx[order(runif(Nsel))]
  return(NewChrIx)
}