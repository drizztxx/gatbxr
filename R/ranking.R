#' @title Rank-based fitness assignment
#'
#' @description
#' This function ranks individuals represented by their associated
#' cost, to be \emph{minimized}, and returns a vector \code{FitnV}
#' containing the corresponding individual fitnesses. For multiple
#' subpopulations the ranking is performed separately for each
#' subpopulation.
#'
#' @param ObjV a vector containing the objective values of the
#' individuals in the current population (cost values).
#' @param RFun if \code{RFun} is a number in \eqn{[1, 2]}, linear ranking is
#' assumed and the value indicates the selective pressure.
#' If \code{RFun} is a vector of 2 elements, then
#' \code{RFun[1]} indicates the selective pressure and
#' \code{RFun[2]} indicates the ranking method: when \code{RFun[2] = 0},
#' linear ranking is used; when \code{RFun[2] = 1},
#' non-linear ranking is used.
#' If \code{RFun} is a vector with \code{length(RFun) > 2}, it contains
#' the fitness to be assigned to each rank. It should have
#' the same length as the number of individuals per subpopulation.
#' Usually \code{RFun} is monotonically increasing.
#' If \code{RFun} is \code{NULL}, linear ranking
#' and a selective pressure of 2 are assumed.
#' @param SUBPOP an optional number indicating subpopulations.
#' Default is set to 1 subpopulation.
#'
#' @details
#' Individuals with equal objective values receive the same fitness,
#' namely the average of the fitness values assigned to their ranks.
#'
#' @return
#' a vector containing the fitness values of the
#' individuals in the current population.
#' @export
#' @author
#' The original 'MATLAB' implementation of ranking was written by Hartmut Pohlheim and
#' Carlos Fonseca.
#' The R implementation was written by Liang Zhao.
#' @seealso \code{\link{scaling}}, \code{\link{select}}
#' @examples
#' Chrom = crtbp(40,10)$Chrom
#'
#' ## Objective function is the sum of individuals
#' ObjV = apply(Chrom,1,sum)
#'
#' ## Linear ranking with selective pressure 2
#' FitnV = ranking(ObjV)
#'
#' ## Non-linear ranking with selective pressure 3
#' FitnV = ranking(ObjV, RFun = c(3,1))

ranking <- function(ObjV,
                    RFun=c(2,0),
                    SUBPOP=1){

  if(!is.vector(ObjV)) stop("ObjV must be a vector")
  Nind <- NROW(ObjV)
  if(length(SUBPOP) != 1) stop("SUBPOP must be a scalar")
  if(Nind %% SUBPOP != 0) stop("ObjV and SUBPOP disagree")
  Nind <- Nind/SUBPOP

  if(is.null(RFun)) RFun <- c(2,0)
  if(!is.vector(RFun)) stop("RFun must be defined as NULL or a vector")

  if(length(RFun) == 2) {
    if(RFun[2] %in% c(0,1)) NonLin <- RFun[2] else stop("Parameter for ranking method must be 0 or 1")
  } else if(length(RFun) == 1) {
    NonLin <- 0
  } else {
    if (length(RFun) != Nind) stop("ObjV and RFun disagree")
    Prss <- RFun
    NonLin <- NA
  }

  if(isTRUE(NonLin == 0)){ ##Linear method
    if(RFun[1] < 1 || RFun[1] > 2) stop("Selective pressure for linear ranking must be between 1 and 2")
    Prss <- (2-RFun[1]) + 2*(RFun[1]-1)*(0:(Nind-1))/(Nind-1)
  } else if(isTRUE(NonLin == 1)){ ##Non-Linear method
    if(RFun[1] < 1) stop("Selective pressure must be greater than 1")
    if(RFun[1] > (Nind - 2)) stop("Selective pressure too big")
    Root1 <- polyroot(c(RFun[1]*rep(1,Nind-1),RFun[1]-Nind))
    Prss <- Re(Root1)[which.min(abs(Im(Root1)))] ^ (0:(Nind-1))
    Prss <- Prss / sum(Prss) * Nind
  }

  FitnV <- numeric(0)
  for(ix in 1:SUBPOP){
    ObjVSub <- ObjV[((ix-1)*Nind+1):(ix*Nind)]
    ## The worst individual gets Prss[1], the best gets Prss[Nind];
    ## individuals with equal objective values share the mean fitness
    FitnVSub <- Prss[rank(-ObjVSub, ties.method = "first")]
    FitnVSub <- ave(FitnVSub, ObjVSub)
    FitnV <- c(FitnV,FitnVSub)
  }
  return(FitnV)
}
