#' @title Discrete mutation (low-level function)
#'
#' @description
#' This function takes the representation of the current population,
#' mutates each element with given probability and returns the resulting
#' population.
#'
#' @param OldChrom a matrix containing the chromosomes of the
#' current population. Each row corresponds to
#' an individual's string representation.
#' @param Pm a value indicating mutation probability. Default value
#' is \code{Pm = 0.7/Lind}, where \code{Lind} is the chromosome length.
#' @param BaseV	an optional vector of the same length as the
#' chromosome structure defining the base of the
#' individual elements of the chromosome. Default is set to binary
#' representation.
#'
#' @return
#' a matrix containing a mutated version of OldChrom.
#' @export
#' @author
#' The original 'MATLAB' implementation of mut was written by Andrew Chipperfield.
#' The R implementation was written by Liang Zhao.
#' @seealso \code{\link{mutate}}, \code{\link{mutbga}}
#' @examples
#' Chrom = crtbp(40,10)$Chrom
#'
#' NewChrom = mut(OldChrom=Chrom)

mut <- function(OldChrom,
                Pm=NULL,
                BaseV=NULL){
  ## get population size (Nind) and chromosome length (Lind)
  Nind <- NROW(OldChrom);Lind <- NCOL(OldChrom)

  ## check input parameters
  if(is.null(Pm)) Pm <- 0.7/Lind
  if(is.null(BaseV)) BaseV <- crtbase(Lind)

  if(Lind != length(BaseV)) stop("OldChrom and BaseV are incompatible")

  ## create mutation mask matrix
  BaseM <- matrix(rep(1,Nind)) %*% BaseV

  ## perform mutation on chromosome structure
  Mask <- (matrix(runif(Nind*Lind),Nind,Lind)<
             Pm)*ceiling(matrix(runif(Nind*Lind),Nind,Lind)*(BaseM-1))
  NewChrom <- (OldChrom+Mask) %% BaseM

  return(NewChrom)
}