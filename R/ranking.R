## RANKING.R      (RANK-based fitness assignment)
##
## This function performs ranking of individuals.
##
## Syntax:  FitnV = ranking(ObjV, RFun, SUBPOP)
##
## This function ranks individuals represented by their associated
## cost, to be *minimized*, and returns a column vector FitnV
## containing the corresponding individual fitnesses. For multiple
## subpopulations the ranking is performed separately for each
## subpopulation.
##
## Input parameters:
##    ObjV      - Vector containing the objective values of the
##                individuals in the current population (cost values).
##    RFun      - (optional) If RFun is a scalar in [1, 2] linear ranking is
##                assumed and the scalar indicates the selective pressure.
##                If RFun is a 2 element vector:
##                RFun(1): SP - scalar indicating the selective pressure
##                RFun(2): RM - ranking method
##                         RM = 0: linear ranking
##                         RM = 1: non-linear ranking
##                If RFun is a vector with length(RFun) > 2 it contains
##                the fitness to be assigned to each rank. It should have
##                the same length as ObjV. Usually RFun is monotonously
##                increasing.
##                If RFun is omitted or NaN, linear ranking
##                and a selective pressure of 2 are assumed.
##    SUBPOP    - (optional) Number of subpopulations
##                if omitted or NaN, 1 subpopulation is assumed
##
## Output parameters:
##    FitnV     - Column vector containing the fitness values of the
##                individuals in the current population.
##                
##
## Author:     Hartmut Pohlheim (Carlos Fonseca)
##             David Zhao (Modified for R)
##
## Date: 13May2016

ranking <- function(ObjV,
                    RFun=c(2,0),
                    SUBPOP=1){
  # ObjV = c(NaN,1:10,NaN,20:16);RFun=c(2,0);SUBPOP=1 ##Test
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
    FitnVSub <- Prss[rank(-rank(ObjV))]
    FitnV <- c(FitnV,FitnVSub)
  }
  return(FitnV)
}