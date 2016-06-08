## MUT.R
##
## This function takes the representation of the current population,
## mutates each element with given probability and returns the resulting
## population.
##
## Syntax:	NewChrom = mut(OldChrom,Pm,BaseV)
##
## Input parameters:
##
##		OldChrom - A matrix containing the chromosomes of the
##			         current population. Each row corresponds to
##			         an individuals string representation.
##
##		Pm	     - Mutation probability (scalar). Default value
##			         of Pm = 0.7/Lind, where Lind is the chromosome
##			         length is assumed if omitted.
##
##		BaseV	   - Optional row vector of the same length as the
##			         chromosome structure defining the base of the 
##			         individual elements of the chromosome. Binary
##			         representation is assumed if omitted.
##
## Output:
##
##		           A Matrix containing a mutated version of
##			         OldChrom.
##
## Author: Andrew Chipperfield
##         David Zhao (Modified for R)
##
## Date: 13May2016

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