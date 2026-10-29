#' @title Line recombination with mutation features
#'
#' @description
#' This function performs line recombination with mutation features between
#' pairs of individuals.
#' @inheritParams mutate
#' @param MutOpt  an optional vector containing recombination rate and shrink value:
#' \enumerate{
#'   \item MutOpt[1] a number containing the recombination rate in the range \eqn{[0, 1]}.
#' By default this value is set to 1.
#'   \item MutOpt[2] a number for shrinking the
#' mutation range in the range \eqn{[0, 1]}.
#' By default this value is assumed to be 1 (no shrinking).
#' }
#' @details
#' \code{recmut} is a function only applicable to populations of real-value variables
#' (and not binary or integer). The pairs are mated in order, odd row with the next even row.
#' If the number of rows in the matrix \code{OldChrom} is odd then the last row is not
#' mated and added at the end of \code{NewChrom}. The population should therefore be organised
#' into contiguous pairs that require mating. This can be achieved by using the function
#' \code{\link{ranking}} to assign a fitness level to each individual and a selection function
#' (e.g. \code{\link{select}}) to select individuals with a probability related to their
#' fitness in the current population.
#' \cr\cr
#' The offspring of a pair of two parents are computed as follows:\cr
#' \eqn{offspring1 = parent1 + RecMx * range * MutOpt[2] * delta * Diff}\cr
#' \eqn{offspring2 = parent2 + RecMx * range * MutOpt[2] * delta * -Diff}\cr
#' where,\cr \eqn{RecMx} = 1 or -1 with probability \eqn{MutOpt[1]} (- with
#' probability 0.9), else 0. \cr
#' \eqn{range = 0.5 *} domain of variable (search interval defined
#' by \code{FieldDR}). \cr
#' \eqn{delta = sum(0 to m-1)pi*2^-i, pi = 1} with probability
#' \eqn{1/m}, else 0, m = 20.\cr
#' \eqn{Diff = (parent2 - parent1)/|norm(parent2) - norm(parent1)|}.
#' Pairs of parents with equal norms are not moved.
#' \cr\cr
#' The recombination operator \code{recmut} generates offspring in a direction defined by
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
#' The original 'MATLAB' implementation of recmut was written by Hartmut Pohlheim.
#' The R implementation was written by Liang Zhao.
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
    warning("Parameter for recombination rate must be in [0, 1], this parameter has been temporarily set to 1")
  }

  if (MutShrink < 0 || MutShrink > 1) {
    MutShrink <- 1
    warning("Parameter for shrinking mutation range must be in [0, 1], this parameter has been temporarily set to 1, which means no shrinking")
  }

  ## Identify the number of matings
  Xops <- floor(Nind/2)

  ## perform mutation
  ## get range matrix for every variable
  Range <- as.vector(0.5 * MutShrink * diff(FieldDR))
  upper <- matrix(rep(FieldDR[2,],Nind),Nind,Nvar,byrow = TRUE)
  lower <- matrix(rep(FieldDR[1,],Nind),Nind,Nvar,byrow = TRUE)

  ## zeros and ones for mutate or not this variable, together with Range
  Range <- do.call(rbind,lapply(1:Xops,function(x) Range * (runif(Nvar) < MutR)))

  ## compute, if + or - sign
  Range <- Range * (1 - 2 * (matrix(runif(Xops*Nvar),Xops,Nvar) < 0.9))

  odd <- seq(1, by = 2, length.out = Xops)
  even <- odd + 1

  ## compute distance between mating pairs
  NormO <- abs(sqrt(rowSums(OldChrom[even,,drop=FALSE]^2)) -
                 sqrt(rowSums(OldChrom[odd,,drop=FALSE]^2)))

  ## compute difference between variables divided by distance
  ## (pairs with equal norms are not moved)
  ChromDiff <- (OldChrom[even,,drop=FALSE] - OldChrom[odd,,drop=FALSE]) / NormO
  ChromDiff[!is.finite(ChromDiff)] <- 0

  ## compute delta value for all individuals
  ACCUR <- 20
  Vect <- 2 ^ (-1 * (0:(ACCUR-1)))
  Delta <- (matrix(runif(Xops*ACCUR),Xops,ACCUR) < 1/ACCUR) %*% Vect
  Delta <- matrix(rep(Delta,Nvar),Xops,Nvar)

  ## perform recombination
  NewChrom <- matrix(NA,Nind,Nvar)
  NewChrom[odd,] <- OldChrom[odd,,drop=FALSE] + Range * Delta * ChromDiff
  NewChrom[even,] <- OldChrom[even,,drop=FALSE] + Range * Delta * (-ChromDiff)

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