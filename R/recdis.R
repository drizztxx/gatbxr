#' @title Discrete recombination
#' 
#' @description
#' This function performs discrete recombination between
#' pairs of individuals and returns the new individuals after mating.
#'
#' @usage
#' recdis(OldChrom)
#'
#' @param OldChrom  a matrix containing the chromosomes of the old
#' population. Each line corresponds to one individual.
#' @details
#' \code{recdis} is a function applicable to both populations of real-value variables 
#' and binary or integer variables. The pairs are mated in order, odd row with the next even row.
#' If the number of rows in the matrix \code{OldChrom} is odd the the last row is not 
#' mated and added at the end of \code{NewChrom}. The population should therefore be organised 
#' into contiguous pairs that require mating. This can be achieved by using the function
#' \code{\link{ranking}} to assign a fitness level to each individual and a selection function
#' (e.g. \code{\link{select}}) to select individuals with a probability related to their
#' fitness in the current population.
#' @section Algorithm:
#' Descrete recombination exchanges variable valules between the individuals. For 
#' each variable the parent who contributes its variable value to the offspring is 
#' chosen randomly with equal probability.
#' \cr\cr
#' Discrete recombination can generate the corners of the hypercube defined by the parents.
#' @return
#' a matrix containing the chromosomes of the population
#' after mating, ready to be mutated and/or evaluated,
#' in the same format as OldChrom.
#' @seealso
#' \code{\link{recombin}}, \code{\link{recint}}, \code{\link{reclin}}, \code{\link{ranking}}, 
#' \code{\link{sus}}, \code{\link{rws}}
#' @export
#' @author 
#' The original matlab implementation was written by Hartmut Pohlheim.
#' The R implementation was written by David Zhao. 
#' @examples
#' ## Initial a real-valued population
#' FieldDR = matrix(c(-100,-50,-30,-20,100,50,30,20),2,4,byrow=T) 
#' Chrom = crtrp(6,FieldDR)
#' 
#' ## Perform discrete recombination
#' NewChrom = recdis(Chrom)

recdis <- function(OldChrom){
  ## Identify the population size (Nind) and the number of variables (Nvar)
  Nind <- NROW(OldChrom); Nvar <- NCOL(OldChrom)
  
  ## Identify the number of matings
  Xops <- floor(Nind/2)
  
  ## Define which parent gives the value
  Mask1 <- matrix(runif(Xops * Nvar),Xops,Nvar) < 0.5
  Mask2 <- matrix(runif(Xops * Nvar),Xops,Nvar) < 0.5
  
  ## Performs recombination
  NewChrom <- matrix(NA,Nind,Nvar)
  odd <- (1:(Nind-1))[as.logical(1:(Nind-1) %% 2)]
  even <- (1:Nind)[!as.logical(1:Nind %% 2)]
  NewChrom[odd,] <- (OldChrom[odd,] * Mask1) + (OldChrom[even,] * !Mask1)
  NewChrom[even,] <- (OldChrom[odd,] * Mask2) + (OldChrom[even,] * !Mask2)
 
  ## If the number of individuals is odd, the last individual cannot be mated
  ## but must be included in the new population
  if(Nind %% 2){
    NewChrom[Nind,] <- OldChrom[Nind,]
  }
  return(NewChrom)
}