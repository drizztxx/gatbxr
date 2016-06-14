#' @title Linear fitness scaling
#' @description
#' This function implements a linear fitness scaling algorithm. For multiple
#' subpopulations \code{scaling} is performed separately for each
#' subpopulation.
#' @usage
#' scaling(ObjV, Smul = 2, SUBPOP = 1)
#' 
#' @param Objv a vector containing the values of individuals fitness.
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
#' the resulting fitness value of individual \eqn{x_i}{xi}. If \eqn{f_ave}{fave} is the average
#' objective value in the current generation, then the maximum fitness of the scaled population
#' is upper bounded at \eqn{f_ave * Smul}{fave * Smul}. If \code{Smul} is ommited the the default
#' value is set to 2. The average fitness of the scaled population is also set to \eqn{f_ave}{fave}.
#' \cr \cr
#' In the case of some of the objective values being negative, scailing attempts to 
#' provide an offset, \eqn{b}, such that the scaled fitness values are greater than zero.
#' @references
#' Goldberg, D. E. (1989) \emph{Genetic Algorithms in Search, Optimization and Machine
#' Learning}. Addison Wesley.
#' @return
#' a vector containing the individual fitnesses
#' for the current population.
#' @export
#' @author 
#' The original matlab implementation of scaling was written by Andrew Chipperfield.
#' The R implementation was written by David Zhao. 
#' @seealso
#' \code{\link{ranking}}, \code{\link{reins}}, \code{\link{rws}}, \code{\link{select}}, 
#' \code{\link{sus}}
#' @note
#' \code{scaling} is not recommended when fitness
#' functions produce negative results as it will become unreliable.
#' It is included in this version of package only for the sake of
#' completeness.
#' @examples
#' 
#' Chrom = crtbp(40,10)$Chrom
#' 
#' #Ojbective function is the sum of individuals
#' ObjV = apply(Chrom,1,sum)
#' FitnV = scaling(ObjV)

scaling <- function(ObjV, Smul=2, SUBPOP=1){
  ## Check parameter consistency
  if(!is.vector(ObjV)) stop("Objv must be a vector")
  if(min(ObjV) < 0) stop("all vaules that objective function returns must be greater than 0")
  Nind <- NROW(ObjV)
  if(length(SUBPOP) != 1) stop("SUBPOP must be a scalar")
  if(Nind %% SUBPOP) stop("ObjV and SUBPOP disagree")
  Nind <- Nind/SUBPOP  

  FitnV <- vector()
  for(ix in 1:SUBPOP){
    ObjVSub <- ObjV[((ix-1)*Nind+1):(ix*Nind)]
    Oave <- sum(ObjVSub)/Nind
    Omin <- min(ObjVSub)
    Omax <- max(ObjVSub)
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
    FitnV <- c(FitnV,FitnVSub)
  }
  return(FitnV)
}