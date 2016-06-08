#' @title RE-INSertion of offspring in population replacing parents
#'
#' @description
#' This function reinserts offspring in the population.
#'
#' @usage
#' reins(Chrom,SelCh,SUBPOP=1,InsOpt=c(0,1),ObjVCh=NULL,ObjVSel=NULL)
#'
#' @param Chrom matrix containing the individuals (parents) of the current
#' population. Each row corresponds to one individual.
#' @param SelCh matrix containing the offspring of the current
#' population. Each row corresponds to one individual.
#' @param SUBPOP an optional number indicating subpopulations
#' 1 subpopulation is as default.
#' @param InsOpt  an optional vector containing the insertion method parameters
#' InsOpt[1] number indicating kind of insertion. 0 - uniform insertion and
#' 1 - fitness-based insertion. if omitted, 0 is assumed.
#' InsOpt[2] rate of offspring to be inserted per subpopulation (% of subpopulation)
#' if omitted, 1.0 (100%) is assumed
#' @param ObjVCh  an optional vector containing the objective values
#' of the individuals (parents - Chrom) in the current 
#' population, needed for fitness-based insertion
#' saves recalculation of objective values for population
#' @param ObjVSel an optional vector containing the objective values
#' of the offspring (SelCh) in the current population, needed for
#' partial insertion of offspring,
#' saves recalculation of objective values for population
#'
#' @return a list containing following components:
#' \item{Chrom}{matrix containing the individuals of the current
#'  population after reinsertion.}
#' \item{ObjVCh}{if ObjVCh and ObjVSel are input parameters, then column 
#'  vector containing the objective values of the individuals
#'  of the current generation after reinsertion.}
#' @author 
#' The original matlab implementation of bs2rv was written by Hartmut Pohlheim and 
#' tested by Alex Shenfield. 
#' The R implementation was written by David Zhao.       

reins <- function(Chrom,
                  SelCh,
                  SUBPOP=1,
                  InsOpt=c(0,1),
                  ObjVCh=NULL,
                  ObjVSel=NULL){
  ## Check parameter consistency
  NindP <- NROW(Chrom); NvarP <- NCOL(Chrom)
  NindO <- NROW(SelCh); NvarO <- NCOL(SelCh)
  
  if (length(SUBPOP) != 1) stop("SUBPOP must be a scalar")
  if (NindP %% SUBPOP) stop("Chrom and SUBPOP disagree")
  if (NindO %% SUBPOP) stop("SelCh and SUBPOP disagree")
  NIND <- NindP/SUBPOP  ## Compute number of individuals per subpopulation
  NSEL <- NindO/SUBPOP  ## Compute number of offspring per subpopulation
  
  if (is.null(ObjVCh)){
    IsObjVCh <- 0
  } else {
    if (!is.vector(ObjVCh)) stop("ObjVCh must be a vector")
    if (NindP != length(ObjVCh)) stop("Chrom and ObjVCh disagree")
    IsObjVCh <- 1
  }
  if (is.null(ObjVSel)){
    IsObjVSel <- 0
  } else {
    if (!is.vector(ObjVSel)) stop("ObjVSel must be a vector")
    if (NindO != length(ObjVSel)) stop("Chrom and ObjVSel disagree")
    IsObjVSel <- 1
  }
  
  if (length(InsOpt) > 2) stop("Parameter InsOpt too long")
  Select <- InsOpt[1]
  INSR <- InsOpt[2]
  if (INSR < 0 || INSR > 1) stop("Parameter for insertion rate must be a scalar in [0, 1]")
  if (!INSR %in% c(0,1) && IsObjVSel != 1) stop("For selection of offspring ObjVSel is needed")
  if (!Select %in% c(0,1)) stop("Parameter for selection method must be 0 or 1")
  if (Select == 1 && IsObjVCh == 0) stop("ObjVCh for fitness-based exchange needed")
  
  if (INSR == 0) {
    if (!IsObjVCh) ObjVCh <- NULL
    return(list(ObjVCh=ObjVCh,Chrom=Chrom))
  } 
  NIns = min(max(floor(INSR*NSEL+0.5),1),NIND)  ## Number of offspring to insert
  
  ## perform insertion for each subpopulation
  for(ix in 1:SUBPOP){
    ## Calculate positions in old subpopulation, where offspring are inserted
    if (Select == 1){
      ChIx <- order(-ObjVCh[((ix-1)*NIND+1):(ix*NIND)])
    } else {
      ChIx <- order(runif(NIND))
    }
    PopIx <- ChIx[1:NIns] + (ix-1)*NIND
    ## Calculate position of Nins-% best offspring
    if (NIns < NSEL){  ## select best offspring
      OffIx <- order(ObjVSel[((ix-1)*NSEL+1):(ix*NSEL)])
    } else {
      OffIx <- 1:NIns
    }
    SelIx <- OffIx[1:NIns] + (ix-1)*NSEL
    ## Insert offspring in subpopulation -> new subpopulation
    Chrom[PopIx,] <- SelCh[SelIx,]
    if (IsObjVCh == 1 && IsObjVSel == 1) ObjVCh[PopIx] <- ObjVSel[SelIx]
  }
  if (!IsObjVCh) ObjVCh <- NULL
  return(list(ObjVCh=ObjVCh,Chrom=Chrom))
}