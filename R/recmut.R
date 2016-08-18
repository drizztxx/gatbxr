#' @title Line recombination with mutation features
#' 
#' @description
#' This function performs line recombination with mutation features between
#' pairs of individuals.
#' @usage
#' recmut(OldChrom, FieldDR, MutOpt = NULL)
#' @inheritParams mutate
#' @param MutOpt  an optional vector containing mutation rate and shrink value:
#' \enumerate{
#'   \item MutOpt[1] a number containing the recombination rate in the range \eqn{[0, 1]}.
#' By default this value is set to 1.
#'   \item MutOpt[2] a number for shrinking the
#' mutation range in the range \eqn{[0, 1]}. 
#' By default this vaule is assumed to 1 (no shrinking).
#' }
#' @details
#' \code{recmut} is a function only applicable to populations of real-value variables (
#' and not binary or integer). The pairs are mated in order, odd row with the next even row.
#' If the number of rows in the matrix \code{OldChrom} is odd the the last row is not 
#' mated and added at the end of \code{NewChrom}. The population should therefore be organised 
#' into contiguous pairs that require mating. This can be achieved by using the function
#' \code{\link{ranking}} to assign a fitness level to each individual and a selection function
#' (e.g. \code{\link{select}}) to select individuals with a probability related to their
#' fitness in the current population.
#' \cr\cr
#' The offsprings of a pair of two parents are computed as follows:\cr
#' \eqn{offsping1 = parent1 + RecMx * range * Mutopt[2] * delta *Diff}\cr
#' \eqn{offsping2 = parent2 + RecMx * range * Mutopt[2] * delta * -Diff}\cr
#' where,\cr \eqn{MutMx} = 1 or -1 with probability \eqn{MutOpt[1]}, (- with 
#' probability 0.9) else 0. \cr
#' \eqn{range = 0.5 *} domain of variable (search interval defined
#' by \code{FieldDR}). \cr
#' \eqn{delta = sum(0 to m-1)pi*2^-i, pi = 1} wiht probability 
#' \eqn{1/m}, else 0, m = 20.\cr
#' \eqn{Diff = (parent2 - parent1)/absolute(parent1 - parent2)}.
#' \cr\cr
#' The recombiantion operator \code{recmut} generates offspring in a direction defined by
#' the parents (line recombination). It tests more often outside the area defined by the
#' parents and in the direction of parent1. The point for the offspring is defined by
#' features of the mutation operator. The probability of small step sizes is greater than
#' that of bigger steps (see \code{\link{mutbga}}).
#' @seealso
#' \code{\link{mutate}}, \code{\link{mutbga}}, \code{\link{reclin}}
#' @return
#' a matrix containing the chromosomes of the population
#' after mating, ready to be mutated and/or evaluated,
#' in the same format as OldChrom.
#' @export
#' @note
#' This function doesn't work with high level recombination function \code{\link{recombin}}
#' @author 
#' The original matlab implementation of crtrp was written by Hartmut Pohlheim. 
#' The R implementation was written by David Zhao. 
#' @examples
#' ## Initial a real-valued population
#' FieldDR = matrix(c(-100,-50,-30,-20,100,50,30,20),2,4,byrow=TRUE) 
#' Chrom = crtrp(6,FieldDR)
#' 
#' ## Perform line recombination with mutation
#' NewChrom = recmut(Chrom,FieldDR)

recmut <- function(OldChrom,FieldDR,MutOpt = NULL){
  
  ## Parameter consistency check
  Nind <- NROW(OldChrom); Nvar <- NCOL(OldChrom)
  
  mF <- NROW(FieldDR); nF <- NCOL(FieldDR)
  if (mF != 2) stop("FieldDR must be a matrix with 2 rows")
  if (Nvar != nF) stop("FieldDR and OldChrom disagree")
  
  if (is.null(MutOpt)) {
    MutR <- 1; MutShrink <- 1
  } else if (length(MutOpt) == 1) {
    MutR <- MutOpt; MutShrink <- 1
  } else if (length(MutOpt) == 2){
    MutR <- MutOpt[1]; MutShrink <- MutOpt[2]
  } else stop("Too many parameters in MutOpt")
  
  if (MutR < 0 || MutR > 1) {
    MutR <- 1
    warning("Parameter for recombination rate must be in [0, 1], this parameter has
            been temporarily set to 1")
  }
  
  if (MutShrink < 0 || MutShrink > 1) {
    MutShrink <- 1
    warning("Parameter for shrinking mutation range must be in [0, 1], this parameter has
            been temporarily set to 1, which means no shrinking")
  }  
  
  ## Identify the number of matings
  Xops <- floor(Nind/2)
  
  ## perform mutation
  ## get range matrix for every variable
  Range <- as.vector(0.5 * MutShrink * diff(FieldDR))
  upper <- matrix(rep(FieldDR[2,],Nind),Nind,Nvar,byrow = T)
  lower <- matrix(rep(FieldDR[1,],Nind),Nind,Nvar,byrow = T)
  
  ## zeros and ones for mutate or not this variable, together with Range
  Range <- do.call(rbind,lapply(1:Xops,function(x) Range * (runif(Nvar) < MutR)))
  
  ## compute, if + or - sign
  Range <- Range * (1 - 2 * (matrix(runif(Xops*Nvar),Xops,Nvar) < 0.9))  
  
  ## compute distance between mating pairs
  NormO <- sapply(1:Xops,function(x) abs(norm(as.matrix(OldChrom[2*x,]),type="2") -
                    norm(as.matrix(OldChrom[2*x - 1,]),type="2")))
  
  ## compute difference between variables divided by distance
  ChromDiff <- do.call(rbind,lapply(1:Xops,function(x) diff(OldChrom[(2*x-1):(2*x),])/NormO[x]))
  
  ## compute delta value for all individuals
  ACCUR <- 20
  Vect <- 2 ^ (-1 * (0:(ACCUR-1)))
  Delta <- (matrix(runif(Xops*ACCUR),Xops,ACCUR) < 1/ACCUR) %*% Vect
  Delta <- matrix(rep(Delta,Nvar),Xops,Nvar)
  
  ## perform recombination
  NewChrom <- matrix(NA,Nind,Nvar)
  odd <- (1:(Nind-1))[as.logical(1:(Nind-1) %% 2)]
  even <- (1:Nind)[!as.logical(1:Nind %% 2)]
  NewChrom[odd,] <- OldChrom[odd,] + Range * Delta * ChromDiff
  NewChrom[even,] <- OldChrom[even,] + Range * Delta * ChromDiff
  
  ## If the number of individuals is odd, the last individual cannot be mated
  ## but must be included in the new population
  if(Nind %% 2){
    NewChrom[Nind,] <- OldChrom[Nind,]
  }
  
  ## Ensure variables boundaries, compare with lower and upper boundaries
  upper.logic <- NewChrom > upper
  lower.logic <- NewChrom < lower
  if (any(upper.logic)){
    NewChrom <- NewChrom * !upper.logic + upper * upper.logic
  }
  if (any(lower.logic)){
    NewChrom <- NewChrom * !lower.logic + lower * lower.logic
  }  
  return(NewChrom)
}