#' @title Reinsertion of offspring in population replacing parents
#'
#' @description
#' This function reinserts offspring in the population.
#'
#' @param Chrom matrix containing the individuals (parents) of the current
#' population. Each row corresponds to one individual.
#' @param SelCh matrix containing the offspring of the current
#' population. Each row corresponds to one individual.
#' @param SUBPOP an optional number indicating subpopulations.
#' Default is set to 1 subpopulation.
#' @param InsOpt  an optional vector containing the insertion method parameters.
#' \code{InsOpt[1]} is a number indicating the kind of insertion: 0 stands for
#' uniform insertion and 1 for fitness-based insertion. If omitted, 0 is assumed.
#' \code{InsOpt[2]} is the rate of offspring to be inserted per subpopulation.
#' If omitted, 1.0 is assumed.
#' @param ObjVCh  an optional vector containing the objective values
#' of the individuals (parents, \code{Chrom}) in the current
#' population, needed for fitness-based insertion. It
#' saves recalculation of objective values for the population.
#' @param ObjVSel an optional vector containing the objective values
#' of the offspring (\code{SelCh}) in the current population, needed for
#' partial insertion of offspring. It
#' saves recalculation of objective values for the population.
#'
#' @return a list containing following components:
#' \item{ObjVCh}{if \code{ObjVCh} is given, the
#'  vector containing the objective values of the individuals
#'  of the current generation after reinsertion (updated only if
#'  \code{ObjVSel} is also given); otherwise \code{NULL}.}
#' \item{Chrom}{matrix containing the individuals of the current
#'  population after reinsertion.}
#' @export
#' @author
#' The original 'MATLAB' implementation of reins was written by Hartmut Pohlheim and
#' tested by Alex Shenfield.
#' The R implementation was written by Liang Zhao.
#' @seealso \code{\link{select}}, \code{\link{recombin}}, \code{\link{mutate}}
#' @examples
#' ## Parents and their objective values (to be minimized)
#' Chrom = crtbp(10,8)$Chrom
#' ObjVCh = apply(Chrom,1,sum)
#'
#' ## Offspring and their objective values
#' SelCh = crtbp(6,8)$Chrom
#' ObjVSel = apply(SelCh,1,sum)
#'
#' ## Replace the worst parents with all offspring (fitness-based insertion)
#' res = reins(Chrom, SelCh, InsOpt = c(1,1), ObjVCh = ObjVCh, ObjVSel = ObjVSel)
#'
#' ## Insert the best 50% of the offspring at uniformly random positions
#' res = reins(Chrom, SelCh, InsOpt = c(0,0.5), ObjVSel = ObjVSel)

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
  Select <- if (length(InsOpt) >= 1) InsOpt[1] else 0
  INSR <- if (length(InsOpt) >= 2) InsOpt[2] else 1
  if (INSR < 0 || INSR > 1) stop("Parameter for insertion rate must be a scalar in [0, 1]")
  if (!INSR %in% c(0,1) && IsObjVSel != 1) stop("For selection of offspring ObjVSel is needed")
  if (!Select %in% c(0,1)) stop("Parameter for selection method must be 0 or 1")
  if (Select == 1 && IsObjVCh == 0) stop("ObjVCh for fitness-based exchange needed")

  if (INSR == 0) {
    if (!IsObjVCh) ObjVCh <- NULL
    return(list(ObjVCh=ObjVCh,Chrom=Chrom))
  }
  NIns = min(max(floor(INSR*NSEL+0.5),1),NIND)  ## Number of offspring to insert
  if (NIns < NSEL && IsObjVSel != 1) stop("For selection of offspring ObjVSel is needed")

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
    Chrom[PopIx,] <- SelCh[SelIx,,drop=FALSE]
    if (IsObjVCh == 1 && IsObjVSel == 1) ObjVCh[PopIx] <- ObjVSel[SelIx]
  }
  if (!IsObjVCh) ObjVCh <- NULL
  return(list(ObjVCh=ObjVCh,Chrom=Chrom))
}