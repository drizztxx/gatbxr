#' @title RANK-based fitness assignment
#'
#' @description
#' This function ranks individuals represented by their associated
#' cost, to be *minimized*, and returns a vector FitnV
#' containing the corresponding individual fitnesses. For multiple
#' subpopulations the ranking is performed separately for each
#' subpopulation.
#'
#' @usage
#' ranking(ObjV,RFun=c(2,0),SUBPOP=1)
#'
#' @param ObjV a vector containing the objective values of the
#' individuals in the current population (cost values).
#' @param RFun if set RFun to a number range in [1, 2], linear ranking is
#' assumed and the value indicates the selective pressure.
#' If RFun is set to a 2 elements vector, then
#' RFun[1] indicates the selective pressure and
#' RFun[2] indicates the ranking method, when RFun[2]
#' = 0, linear ranking is used, when RFun[2] = 1, 
#' non-linear ranking is used.
#' If RFun is a vector with length(RFun) > 2, it will contain
#' the fitness to be assigned to each rank. It should have
#' the same length as ObjV. Usually RFun is monotonously
#' increasing.
#' If RFun is omitted, linear ranking
#' and a selective pressure of 2 are assumed.
#' @param SUBPOP an optional number indicating subpopulations.
#' Default is set to 1 subpopulation.
#'
#' @return
#' a vector containing the fitness values of the
#' individuals in the current population.
#' @export
#' @author 
#' The original matlab implementation of ranking was written by Hartmut Pohlheim and
#' Carlos Fonseca.
#' The R implementation was written by David Zhao. 
#' @examples
#' 
#' Chrom = crtbp(40,10)$Chrom
#' 
#' #Ojbective function is the sum of individuals
#' ObjV = apply(Chrom,1,sum)
#' FitnV = ranking(ObjV)

ranking <- function(ObjV,
                    RFun=c(2,0),
                    SUBPOP=1){

  if(!is.vector(ObjV)) stop("Objv must be a vector")
  Nind <- NROW(ObjV)
  if(length(SUBPOP) != 1) stop("SUBPOP must be a scalar")
  if(Nind %% SUBPOP != 0) stop("ObjV and SUBPOP disagree")
  Nind = Nind/SUBPOP
  
  if(!is.null(RFun)){
    if(!is.vector(RFun)) stop("RFun must be defined as NULL or a vector")
    
    if(length(RFun) == 2) {
      if(RFun[2] %in% c(0,1)) NonLin = RFun[2] else stop("Parameter for ranking method must be 0 or 1")
    } else if(length(RFun) == 1) NonLin = 0 else if(length(RFun) > 2) {
      if (length(RFun) != NROW(ObjV)) stop("ObjV and RFun disagree")
      Prss <- RFun
      NonLin <- 9999
    }
    
    if(NonLin == 0){ ##Linear method
      if(RFun[1] < 1 || RFun[1] > 2) stop("Selective pressure for linear ranking must be between 1 and 2")
      Prss <- (2-RFun[1]) + 2*(RFun[1]-1)*(0:(Nind-1))/(Nind-1)
    } else if(NonLin == 1){ ##Non-Linear method
      if(RFun[1] < 1) stop("Selective pressure must be greater than 1")
      if(RFun[1] > (Nind - 2)) stop("Selective pressure too big")
      Root1 <- polyroot(c(RFun[1]*rep(1,Nind-1),RFun[1]-Nind))
      Prss <- Re(Root1)[which.min(abs(Im(Root1)))] ^ (0:(Nind-1))
      Prss= Prss / sum(Prss) * Nind;
    }
  }
  
  FitnV <- vector()
  for(ix in 1:SUBPOP){
    ObjVSub <- ObjV[((ix-1)*Nind+1):(ix*Nind)]
    Prss <- Prss[1:Nind]
    FitnVSub <- Prss[rank(-rank(ObjVSub))]
    FitnV <- c(FitnV,FitnVSub)
  }
  return(FitnV)
}