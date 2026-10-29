#' @title Linear fitness scaling
#' @description
#' This function implements a linear fitness scaling algorithm. For multiple
#' subpopulations \code{scaling} is performed separately for each
#' subpopulation.
#'
#' @param ObjV a vector containing the non-negative objective values of the
#' individuals.
#' @param Smul an optional value used to determine the upper bound,
#' default is set to 2.
#' @param SUBPOP an optional value indicating subpopulations.
#' Default is set to 1 subpopulation.
#' @details
#' \code{scaling} converts the objective values, \code{ObjV}, into a fitness
#' measure with a known upper bound, determined by the value of \code{Smul}, such that,
#' \deqn{F(x_i) = \alpha f(x_i) + b}{%
#' F(xi) = a*f(xi) + b}
#' where \eqn{f(x_i)}{f(xi)} is the objective value of individual \eqn{x_i}{xi},
#' \eqn{\alpha}{a} is a scaling coefficient, \eqn{b} is an offset and \eqn{F(x_i)}{F(xi)} is
#' the resulting fitness value of individual \eqn{x_i}{xi}. If \eqn{f_{ave}}{fave} is the average
#' objective value in the current generation, then the maximum fitness of the scaled population
#' is upper bounded at \eqn{f_{ave} * Smul}{fave * Smul}. If \code{Smul} is omitted then the default
#' value is set to 2. The average fitness of the scaled population is also set to
#' \eqn{f_{ave}}{fave}.
#' \cr \cr
#' If the upper bound cannot be reached without producing negative fitness
#' values, \code{scaling} instead chooses the offset \eqn{b} such that the
#' smallest scaled fitness value is zero. If all objective values of a
#' subpopulation are equal, every individual is assigned the average value.
#' @references
#' Goldberg, D. E. (1989) \emph{Genetic Algorithms in Search, Optimization and Machine
#' Learning}. Addison Wesley.
#' @return
#' a vector containing the individual fitnesses
#' for the current population.
#' @export
#' @author
#' The original 'MATLAB' implementation of scaling was written by Andrew Chipperfield.
#' The R implementation was written by Liang Zhao.
#' @seealso
#' \code{\link{ranking}}, \code{\link{reins}}, \code{\link{rws}}, \code{\link{select}},
#' \code{\link{sus}}
#' @note
#' Negative objective values are not supported, as linear scaling becomes
#' unreliable in that case. \code{\link{ranking}} is generally recommended
#' instead; \code{scaling} is included only for the sake of completeness.
#' @examples
#' Chrom = crtbp(40,10)$Chrom
#'
#' ## Objective function is the sum of individuals
#' ObjV = apply(Chrom,1,sum)
#' FitnV = scaling(ObjV)

scaling <- function(ObjV, Smul=2, SUBPOP=1){
  ## Check parameter consistency
  if(!is.vector(ObjV)) stop("ObjV must be a vector")
  if(min(ObjV) < 0) stop("all values that objective function returns must be non-negative")
  Nind <- NROW(ObjV)
  if(length(SUBPOP) != 1) stop("SUBPOP must be a scalar")
  if(Nind %% SUBPOP) stop("ObjV and SUBPOP disagree")
  Nind <- Nind/SUBPOP

  FitnV <- numeric(0)
  for(ix in 1:SUBPOP){
    ObjVSub <- ObjV[((ix-1)*Nind+1):(ix*Nind)]
    Oave <- sum(ObjVSub)/Nind
    Omin <- min(ObjVSub)
    Omax <- max(ObjVSub)
    if (Omax == Omin){
      FitnVSub <- rep(Oave, Nind)
    } else {
      if (Omin > (Smul * Oave - Omax)/(Smul - 1.0)){
        delta <- Omax - Oave
        a <- (Smul - 1.0) * Oave/delta
        b <- Oave * (Omax - Smul * Oave)/delta
      } else{
        delta <- Oave - Omin
        a <- Oave/delta
        b <- -Omin * Oave/delta
      }
      FitnVSub <- ObjVSub * a + b
    }
    FitnV <- c(FitnV,FitnVSub)
  }
  return(FitnV)
}
