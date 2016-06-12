#' @title CROSSOVer Single-Point
#'
#' @description
#' This function performs single-point crossover between pairs of 
#' individuals and returns the current generation after mating.
#'
#' @usage
#' xovsp(OldChrom,XOVR=0.7)
#'
#' @param OldChrom a matrix containing the chromosomes of the old
#' population. Each line corresponds to one individual.
#' In any form, not necessarily real values.
#' @param XOVR  a number indicating the probability of recombination 
#' occurring between pairs of individuals.
#'
#' @return
#' a matrix containing the chromosomes of the population
#' after mating, ready to be mutated and/or evaluated,
#' in the same format as OldChrom.
#' @export
#' @author 
#' The original matlab implementation of mutate was written by Hartmut Pohlheim.
#' The R implementation was written by David Zhao. 
#' @examples
#' 
#' Chrom = crtbp(40,10)$Chrom
#' 
#' Selch = xovsp(OldChrom=Chrom)

xovsp <- function(OldChrom,
                  XOVR=0.7){
  ## call low level function with appropriate parameters
  NewChrom <- xovmp(OldChrom,XOVR,1,TRUE)
  return(NewChrom)
}