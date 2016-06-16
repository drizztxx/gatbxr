#' @title RECOMBINation high-level function
#'
#' @description
#' This function performs recombination between pairs of individuals
#' and returns the new individuals after mating. The function handles
#' multiple populations and calls the low-level recombination function
#' for the actual recombination process.
#'
#' @usage
#' recombin(REC_F, Chrom, RecOpt = 0.7, SUBPOP = 1, ...)
#'
#' @param REC_F character string containing the name of the recombination or
#' crossover function.
#' @param Chrom a matrix containing the chromosomes of the old
#' population. Each line corresponds to one individual
#' @param RecOpt  an optional value containing the probability of 
#' recombination/crossover occurring between pairs of individuals.
#' Default is set to 0.7.
#' @param SUBPOP  an optional number indicating subpopulations.
#' Default is set to 1.
#' @param ... ohter aurguments passed on to crossover function.
#'
#' @return 
#' a matrix containing the chromosomes of the population
#' after recombination in the same format as OldChrom.
#' 
#' @note
#' This function doesn't work with low level recombination function \code{\link{recmut}}
#' @export
#' @author 
#' The original matlab implementation of recombin was written by Hartmut Pohlheim.
#' The R implementation was written by David Zhao. 
#' 
#' @examples
#' 
#' Selch = crtbp(40,10)$Chrom
#' 
#' NewChrom = recombin(REC_F="xovsp",Chrom=Selch)

recombin <- function(REC_F,
                     Chrom,
                     RecOpt=0.7,
                     SUBPOP=1,
                     ...){
  ## Check parameter consistency
  if (!is.function(get(REC_F))) stop("REC_F must be a string of function name")
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
    ChromSub <- Chrom[((ix-1)*Nind+1):(ix*Nind),]
    NewChromSub <- do.call(REC_F,list(ChromSub,RecOpt,...))
    NewChrom <- rbind(NewChrom,NewChromSub)
  }
  
  return(NewChrom)
}