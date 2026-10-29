#' @title Recombination (high-level function)
#'
#' @description
#' This function performs recombination between pairs of individuals
#' and returns the new individuals after mating. The function handles
#' multiple populations and calls the low-level recombination function
#' for the actual recombination process.
#'
#' @param REC_F character string containing the name of the recombination or
#' crossover function, e.g. \code{"xovsp"} or \code{"recdis"}.
#' @param Chrom a matrix containing the chromosomes of the old
#' population. Each row corresponds to one individual.
#' @param RecOpt  an optional value containing the probability of
#' recombination/crossover occurring between pairs of individuals.
#' It is only used by the crossover functions (\code{xov*}).
#' Default is set to 0.7.
#' @param SUBPOP  an optional number indicating subpopulations.
#' Default is set to 1.
#' @param ... other arguments passed on to the crossover function.
#'
#' @return
#' a matrix containing the chromosomes of the population
#' after recombination in the same format as \code{Chrom}.
#'
#' @note
#' This function doesn't work with low level recombination function \code{\link{recmut}}.
#' @export
#' @author
#' The original 'MATLAB' implementation of recombin was written by Hartmut Pohlheim.
#' The R implementation was written by Liang Zhao.
#' @seealso \code{\link{xovmp}}, \code{\link{recdis}}, \code{\link{recint}},
#' \code{\link{reclin}}
#' @examples
#' Selch = crtbp(40,10)$Chrom
#'
#' NewChrom = recombin(REC_F="xovsp",Chrom=Selch)

recombin <- function(REC_F,
                     Chrom,
                     RecOpt=0.7,
                     SUBPOP=1,
                     ...){
  ## Check parameter consistency
  if (!is.character(REC_F) || length(REC_F) != 1) stop("REC_F must be a string of function name")
  match.fun(REC_F)
  if (REC_F == "recmut") stop("recmut cannot be called by this high level function, please use it directly")
  NindCh <- NROW(Chrom)
  Nvar <- NCOL(Chrom)

  if(length(SUBPOP) != 1) stop("SUBPOP must be a scalar")
  if(NindCh %% SUBPOP) stop("Chrom and SUBPOP disagree")
  Nind = NindCh/SUBPOP ## Compute number of individuals per subpopulation

  if(length(RecOpt) != 1) stop("RecOpt must be a scalar")
  if(RecOpt < 0 || RecOpt > 1) stop("RecOpt must be a scalar in [0, 1]")

  ## Select individuals of one subpopulation and call low level function
  NewChrom <- NULL
  for(ix in 1:SUBPOP){
    ChromSub <- Chrom[((ix-1)*Nind+1):(ix*Nind),,drop=FALSE]
    if (grepl("xov",REC_F)) {
      NewChromSub <- do.call(REC_F,list(ChromSub,RecOpt,...))
      } else NewChromSub <- do.call(REC_F,list(ChromSub,...))
    NewChrom <- rbind(NewChrom,NewChromSub)
  }

  return(NewChrom)
}