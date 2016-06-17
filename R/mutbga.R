#' @title Mutation of real-valued population
#' @description
#' This function implements the mutation operator of the Breeder Genetic
#' Algorithm.
#' @details
#' \code{mutbga} takes a matrix \code{OldChrom} containing the real
#' representation of the individuals in the current population,
#' mutates the individuals with given probability and returns
#' the resulting population.
#' \cr\cr
#' The mutataion of a variable is computed as follows:\cr
#' mutated variable = variable + \eqn{MutMx * range * MutOpt[2] * delta}, 
#' where \eqn{MutMx} = 1 or -1 with probability \eqn{MutOpt[1]}, (+ or - with equal
#' probability) else 0. \eqn{range = 0.5 *} domain of variable (search interval defined
#' by \code{FieldDR}). \eqn{delta = sum(0 to m-1)pi*2^-i, pi = 1} wiht probability 
#' \eqn{1/m}, else 0, m = 20.\cr
#' \cr
#' With \eqn{m = 20}, the mutation operator is able to locate the optimum up to a 
#' precision of \eqn{range * MutOpt[2] * 2^-19}
#' \cr\cr
#' \code{mutbga} is able to generate most points in the hypercube defined by the 
#' variable of the individual and the range of the mutation. However, it test more
#' often near the variable, that is, the probability of small step size is greater
#' than that of larger step sizes.
#' @usage
#' mutbga(OldChrom, FieldDR, MutOpt)
#' @param OldChrom  a matrix containing the chromosomes of the old
#' population. Each line corresponds to one individual.
#' @inheritParams crtrp
#' @inheritParams mutate
#' @return
#' a matrix containing the chromosomes of the population
#' after mutation in the same format as OldChrom.
#' @export
#' @seealso
#' \code{\link{mutate}}, \code{\link{recdis}}, \code{\link{recint}}, \code{\link{recmut}}, 
#' \code{\link{reclin}}
#' @references
#' Muhlenbein, H., Schlierkamp-Voosen, D. (1993) \emph{Predictive Models for the Breeder 
#' Genetic Algorithm: I. Continuous Parameter Optimization}. Evolutionary Computation.
#' @author 
#' The original matlab implementation of crtrp was written by Hartmut Pohlheim. 
#' The R implementation was written by David Zhao. 
#' @examples
#' ## Initial a real-valued population
#' FieldDR = matrix(c(-100,-50,-30,-20,100,50,30,20),2,4,byrow=T) 
#' Chrom = crtrp(6,FieldDR)
#' 
#' ## Perform mutation
#' NewChrom = mutgba(Chrom,FieldDR)

mutbga <- function(OldChrom, FieldDR, MutOpt=NULL){
  
  ## Parameter consistency check
  Nind <- NROW(OldChrom); Nvar <- NCOL(OldChrom)
  
  mF <- NROW(FieldDR); nF <- NCOL(FieldDR)
  if (mF != 2) stop("FieldDR must be a matrix with 2 rows")
  if (Nvar != nF) stop("FieldDR and OldChrom disagree")
  
  if (is.null(MutOpt)) {
    MutR <- 1/Nvar; MutShrink <- 1;
  } else if (length(MutOpt) == 1) {
    MutR <- MutOpt; MutShrink <- 1
  } else if (length(MutOpt) == 2){
    MutR <- MutOpt[1]; MutShrink <- MutOpt[2]
  } else stop("Too many parameters in MutOpt")
  
  if (MutR < 0 || MutR > 1) {
    MutR <- 1/Nvar
    warning("Parameter for mutation rate must be in [0, 1], this parameter has
            been temporarily set to 1/NCOL(FieldDR)")
  }
  
  if (MutShrink < 0 || MutShrink > 1) {
    MutShrink <- 1
    warning("Parameter for shrinking mutation range must be in [0, 1], this parameter has
            been temporarily set to 1, which means no shrinking")
  }  
  
  ## perform mutation
  ## get range matrix for every variable
  Range <- as.vector(0.5 * MutShrink * diff(FieldDR))
  upper <- matrix(rep(FieldDR[2,],Nind),Nind,Nvar,byrow = T)
  lower <- matrix(rep(FieldDR[1,],Nind),Nind,Nvar,byrow = T)
  ## zeros and ones for mutate or not this variable, together with Range
  Range <- do.call(rbind,lapply(1:Nind,function(x) Range * (runif(Nvar) < MutR)))
  ## compute, if + or - sign
  Range <- Range * (1 - 2 * (matrix(runif(Nind*Nvar),Nind,Nvar) < 0.5))
  ## used for later computing, here only ones computed
  ACCUR <- 20
  Vect <- 2 ^ (-1 * (0:(ACCUR-1)))
  Delta <- (matrix(runif(Nind*ACCUR),Nind,ACCUR) < 1/ACCUR) %*% Vect
  Delta <- matrix(rep(Delta,Nvar),Nind,Nvar)
  
  NewChrom <- OldChrom + Range * Delta
  
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