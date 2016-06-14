#' @title Roulette Wheel Selection
#' 
#' @description
#' This function selects a given number of individuals from a population.
#'
#' @usage
#' rws(FitnV, Nsel)
#'
#' @param FitnV a vector containing the fitness values of the
#' individuals in the population.
#' @param Nsel  a number indicating individuals to be selected
#' 
#' @details
#' \code{rws} probabilistically select \code{Nsel} individuals for reproduction according
#' to their fitness, \code{FitnV}, in the current population.
#' 
#' \code{NewChrIx = rws(FitnV,Nsel)} selects \code{Nsel} individuals from a
#' population using roulette wheel selection. \code{FitnV} is a vector containing a
#' performance measure for each individual in the population. This can be achieved
#' by using the function \code{\link{ranking}} or \code{\link{scaling}} to assign a
#' fitness level to each individual.
#' 
#' \code{rws} is a low-level selection function normally called by \code{\link{select}}.
#' 
#' @section Algorithm:
#' A form of roulette wheel selection is implemented by obtaining a cumulative sum
#' of the fitness vector, \code{FitnV}, and generating \code{Nsel} uniformly at random
#' distributed numbers between \code{0} and \code{sum(FitnV)}. The index of the individuals
#' selected is determined by comparing the generated numbers with the cumulative
#' sum vector. The probability of an individual being selected is the given by:
#' \deqn{F(x_i) = \frac{f(x_i)}{\sum_{i=1}^{Nind} f(x_i)}}{% 
#' F(xi) = f(xi)/sum(f(xi))}
#' where \eqn{f(x_i)}{f(xi)} is the fitness of individual \eqn{x_i}{xi} and 
#' \eqn{F(x_i)}{F(xi)} is the probability of that individual being selected.
#' 
#' @return
#' a vector containing the indexes of the selected
#' individuals relative to the original population, shuffled.
#' The new population, ready for mating, can be obtained
#' by calculating OldChrom[NewChrIx,].
#' @export
#' @author 
#' The original matlab implementation of rws was written by Carlos Fonseca and
#' Andrew Chipperfield.
#' The R implementation was written by David Zhao. 
#' @seealso
#' \code{\link{select}}, \code{\link{sus}}, \code{\link{reins}}, \code{\link{ranking}}, 
#' \code{\link{scaling}}
#' @examples
#' ## Consider a population of 8 individuals with the assigned fitness values, FitnV:
#' FitnV = c(1.5,1.35,1.21,1.07,0.92,0.78,0.64,0.5)
#' ## Select the indices of 6 individuals:
#' NewChrIx = rws(FitnV,6)

rws <- function(FitnV, Nsel){
  ## Identify the population size (Nind)
  Nind <- NROW(FitnV)
  
  ## Perform Stochastic Sampling with Replacement
  cumfit <- cumsum(FitnV)
  trials <- cumfit[Nind] * runif(Nsel)
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