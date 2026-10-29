#' @title Mutation (high-level function)
#'
#' @description
#' This function takes a matrix \code{OldChrom} containing the
#' representation of the individuals in the current population,
#' mutates the individuals and returns the resulting population.
#'
#' The function handles multiple populations and calls the low-level
#' mutation function for the actual mutation process.
#'
#' @param MUT_F character string containing the name of the low-level mutation
#' function, e.g. \code{"mut"} (default) or \code{"mutbga"}.
#' @param OldChrom  matrix containing the chromosomes of the old
#' population. Each line corresponds to one individual.
#' @param FieldDR matrix describing the boundaries of each variable
#' (real-values) or defining the base of the variables of
#' each individual (discrete values). Optional for (binary) discrete values
#' @param MutOpt  an optional vector containing mutation rate and shrink value:
#' \enumerate{
#'   \item MutOpt[1] a number containing the mutation rate in the range \eqn{[0, 1]}.
#' By default this value is assumed to \eqn{1/Nvar}, where \eqn{Nvar} is the number
#' of variables per individual defined by \code{NCOL(FieldDR)}.
#'   \item MutOpt[2] a number for shrinking the
#' mutation range in the range \eqn{[0, 1]}.
#' By default this value is assumed to be 1 (no shrinking).
#' }
#' @param SUBPOP an optional number of subpopulations.
#' Default is set to 1.
#' @param ... other arguments passed on to the mutation function.
#'
#' @return
#' a matrix containing the chromosomes of the population
#' after mutation in the same format as \code{OldChrom}.
#' @note
#' \code{mutate} has the same name as functions in other packages (for example
#' \code{dplyr::mutate}); use \code{gatbxr::mutate} if another package masks it.
#' @export
#' @author
#' The original 'MATLAB' implementation of mutate was written by Hartmut Pohlheim.
#' The R implementation was written by Liang Zhao.
#' @seealso \code{\link{mut}}, \code{\link{mutbga}}
#' @examples
#' ## Binary population
#' Chrom = crtbp(40,10)$Chrom
#' NewChrom = mutate(MUT_F="mut",OldChrom=Chrom)
#'
#' ## Real-valued population
#' FieldDR = matrix(c(-100,-50,-30,-20,100,50,30,20),2,4,byrow=TRUE)
#' Chrom = crtrp(6,FieldDR)
#' NewChrom = mutate("mutbga",Chrom,FieldDR)

mutate <- function(MUT_F="mut",
                   OldChrom,
                   FieldDR=NULL,
                   MutOpt=NULL,
                   SUBPOP=1,
                   ...){
  ## Identify the population size (Nind) and the number of variables (Nvar)
  Nind <- NROW(OldChrom); Nvar <- NCOL(OldChrom)

  if (is.null(FieldDR)) {
    IsDiscret = 1
  } else {
    mF <- NROW(FieldDR); nF <- NCOL(FieldDR)
    if (nF != Nvar) stop("FieldDR and OldChrom disagree")
    if (mF == 2) {
      IsDiscret = 0
    } else if (mF == 1){
      IsDiscret = 1
    } else stop("FieldDR must be a matrix with 1 or 2 rows")
  }

  if (length(SUBPOP) != 1) stop("SUBPOP must be a scalar")

  if (Nind %% SUBPOP) stop("OldChrom and SUBPOP disagree")

  Nind <- Nind/SUBPOP ## Compute number of individuals per subpopulation

  ## Select individuals of one subpopulation and call low level function
  NewChrom <- NULL
  for(ix in 1:SUBPOP){
    ChromSub <- OldChrom[((ix-1)*Nind+1):(ix*Nind),,drop=FALSE]
    if (IsDiscret == 1) {
      NewChromSub <- do.call(MUT_F,list(ChromSub,MutOpt,FieldDR,...))
    } else {
      NewChromSub <- do.call(MUT_F,list(ChromSub,FieldDR,MutOpt,...))
    }
    NewChrom <- rbind(NewChrom,NewChromSub)
  }
  return(NewChrom)
}