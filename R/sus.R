#' @title Stochastic Universal Sampling
#'
#' @description
#' This function performs selection with STOCHASTIC UNIVERSAL SAMPLING.
#'
#' @usage
#' sus(FitnV, Nsel)
#'
#' @param FitnV a vector containing the fitness values of the
#' individuals in the population.
#' @param Nsel  a number indicating individuals to be selected
#'
#' @return
#' a vector containing the indexes of the selected
#' individuals relative to the original population, shuffled.
#' The new population, ready for mating, can be obtained
#' by calculating OldChrom[NewChrIx,].
#' @export
#' @author 
#' The original matlab implementation of mutate was written by Hartmut Pohlheim and
#' Carlos Fonseca.
#' The R implementation was written by David Zhao. 
#' @examples
#' 
#' Chrom = crtbp(40,10)$Chrom
#' 
#' #Ojbective function is the sum of individuals
#' FitnV = ranking(apply(Chrom,1,sum))
#' 
#' Selch = sus(FitnV=FitnV,Nsel=2)

sus <- function(FitnV,
                Nsel){
  ## Identify the population size (Nind)
  Nind <- NROW(FitnV)
  
  ## Perform stochastic universal sampling
  cumfit <- cumsum(FitnV)
  trials <- cumfit[Nind]/Nsel * (runif(1) + 0:(Nsel-1))
  NewChrIx <- rep(0,Nsel)
  for(i in 1:Nsel){
    for (j in 1:Nind){
      if (trials[i] < cumfit[j]){
        NewChrIx[i] <- j
        break
      }
    }
  }
  ## Shuffle new population
  NewChrIx <- NewChrIx[order(runif(Nsel))]
  return(NewChrIx)
}