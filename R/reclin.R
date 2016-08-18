#' @title Line recombination
#' 
#' @description
#' This function performs line recombination between
#' pairs of individuals and returns the new individuals after mating.
#'
#' @usage
#' reclin(OldChrom)
#'
#' @param OldChrom  a matrix containing the chromosomes of the old
#' population. Each line corresponds to one individual.
#' @details
#' \code{reclin} is a function only applicable to populations of real-value variables (
#' and not binary or integer). The pairs are mated in order, odd row with the next even row.
#' If the number of rows in the matrix \code{OldChrom} is odd the the last row is not 
#' mated and added at the end of \code{NewChrom}. The population should therefore be organised 
#' into contiguous pairs that require mating. This can be achieved by using the function
#' \code{\link{ranking}} to assign a fitness level to each individual and a selection function
#' (e.g. \code{\link{select}}) to select individuals with a probability related to their
#' fitness in the current population.
#' @section Algorithm:
#' Line recombination combines parent values using the following rule:
#' \deqn{offspring = parent1 + Alpha * (parent2 - parent1)}
#' where \eqn{Alpha} is a scaling factor chosen uniformly at random in the interval
#' \eqn{[-0.25, 1.25]}. \code{reclin} produces one Alpha factor for each pair of parents 
#' combined.
#' \cr\cr
#' Line recombination can generate any point on a slghtly longer 
#' line than that defined by the parents.
#' \cr\cr
#' Line recombination is similar to intermediate recombination \code{\link{recint}}. Whereas 
#' \code{reclin} uses one Alpha factor for each pair of parents combined together,
#' \code{\link{recint}} uses a new Alpha factor for each pair of values.
#' @return
#' a matrix containing the chromosomes of the population
#' after mating, ready to be mutated and/or evaluated,
#' in the same format as OldChrom.
#' @seealso
#' \code{\link{recombin}}, \code{\link{recdis}}, \code{\link{recint}}, \code{\link{ranking}}, 
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
#' ## Perform line recombination
#' NewChrom = reclin(Chrom)

reclin <- function(OldChrom){
  ## Identify the population size (Nind) and the number of variables (Nvar)
  Nind <- NROW(OldChrom); Nvar <- NCOL(OldChrom)
  
  ## Identify the number of matings
  Xops <- floor(Nind/2)
  
  ## Performs recombination
  NewChrom <- matrix(NA,Nind,Nvar)
  odd <- (1:(Nind-1))[as.logical(1:(Nind-1) %% 2)]
  even <- (1:Nind)[!as.logical(1:Nind %% 2)]
  
  ## position of value of offspring compared to parents
  Alpha = matrix(rep(-0.25 + 1.5 * runif(Xops),Nvar),Xops,Nvar)
  NewChrom[odd,] <- OldChrom[odd,] + Alpha * (OldChrom[even,] - OldChrom[odd,])
  
  ## the same ones more for second half of offspring
  Alpha = matrix(rep(-0.25 + 1.5 * runif(Xops),Nvar),Xops,Nvar)
  NewChrom[even,] <- OldChrom[odd,] + Alpha * (OldChrom[even,] - OldChrom[odd,])  
  
  ## If the number of individuals is odd, the last individual cannot be mated
  ## but must be included in the new population
  if(Nind %% 2){
    NewChrom[Nind,] <- OldChrom[Nind,]
  }
  return(NewChrom)
}