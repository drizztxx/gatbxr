#' @title Intermediate recombination
#' 
#' @description
#' This function performs extended intermediate recombination between
#' pairs of individuals and returns the new individuals after mating.
#'
#' @usage
#' recint(OldChrom)
#'
#' @param OldChrom  a matrix containing the chromosomes of the old
#' population. Each line corresponds to one individual.
#' @details
#' \code{recint} is a function only applicable to populations of real-value variables (
#' and not binary or integer). The pairs are mated in order, odd row with the next even row.
#' If the number of rows in the matrix \code{OldChrom} is odd the the last row is not 
#' mated and added at the end of \code{NewChrom}. The population should therefore be organised 
#' into contiguous pairs that require mating. This can be achieved by using the function
#' \code{\link{ranking}} to assign a fitness level to each individual and a selection function
#' (e.g. \code{\link{select}}) to select individuals with a probability related to their
#' fitness in the current population.
#' @section Algorithm:
#' Intermediate recombination combines parent values using the following rule:
#' \deqn{offspring = parent1 + Alpha * (parent2 - parent1)}
#' where \eqn{Alpha} is a scaling factor chosen uniformly at random in the interval
#' \eqn{[-0.25, 1.25]}. \code{recint} produces a new Alpha for each pair of values to be 
#' combined.
#' \cr\cr
#' Intermediate recombination can generate any point within a hypercube slightly
#' larger than that defined by the parents.
#' \cr\cr
#' Intermediate recombination is similar to line recombination \code{\link{reclin}}. Whereas 
#' \code{recint} uses a new Alpha factor for each pair of values combined together,
#' \code{\link{reclin}} uses one Alpha factor for each pair of parents.
#' @return
#' a matrix containing the chromosomes of the population
#' after mating, ready to be mutated and/or evaluated,
#' in the same format as OldChrom.
#' @seealso
#' \code{\link{recombin}}, \code{\link{recdis}}, \code{\link{reclin}}, \code{\link{ranking}}, 
#' \code{\link{sus}}, \code{\link{rws}}
#' @export
#' @author 
#' The original matlab implementation was written by Hartmut Pohlheim.
#' The R implementation was written by David Zhao. 
#' @examples
#' ## Initial a real-valued population
#' FieldDR = matrix(c(-100,-50,-30,-20,100,50,30,20),2,4,byrow=TRUE) 
#' Chrom = crtrp(6,FieldDR)
#' 
#' ## Perform intermeidate recombination
#' NewChrom = recint(Chrom)

recint <- function(OldChrom){
  ## Identify the population size (Nind) and the number of variables (Nvar)
  Nind <- NROW(OldChrom); Nvar <- NCOL(OldChrom)
  
  ## Identify the number of matings
  Xops <- floor(Nind/2)
  
  ## Performs recombination
  NewChrom <- matrix(NA,Nind,Nvar)
  odd <- (1:(Nind-1))[as.logical(1:(Nind-1) %% 2)]
  even <- (1:Nind)[!as.logical(1:Nind %% 2)]
  
  ## position of value of offspring compared to parents
  Alpha = matrix(-0.25 + 1.5 * runif(Xops * Nvar),Xops,Nvar)
  NewChrom[odd,] <- OldChrom[odd,] + Alpha * (OldChrom[even,] - OldChrom[odd,])
  
  ## the same ones more for second half of offspring
  Alpha = matrix(-0.25 + 1.5 * runif(Xops * Nvar),Xops,Nvar)
  NewChrom[even,] <- OldChrom[odd,] + Alpha * (OldChrom[even,] - OldChrom[odd,])  
  
  ## If the number of individuals is odd, the last individual cannot be mated
  ## but must be included in the new population
  if(Nind %% 2){
    NewChrom[Nind,] <- OldChrom[Nind,]
  }
  return(NewChrom)
}