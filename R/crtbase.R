#' @title Create base vector 
#'
#' @description
#' This function creates a vector containing the base of the loci
#' in a chromosome.
#'
#' @usage
#' crtbase(Lind, Base=NULL)
#'
#' @param Lind	a number or vector containing the lengths
#' of the alleles. Sum(Lind) is the length of
#' the corresponding chromosome.
#' @param Base a number or vector containing the base of
#' the loci contained in the Alleles. Default is set to binary
#' representation.
#'
#' @return
#' a vector whose elements correspond to the base
#' of the loci of the associated chromosome structure.
#' @export
#' @author 
#' The original matlab implementation of mutate was written by Andrew Chipperfield.
#' The R implementation was written by David Zhao. 
#' @examples
#' 
#' Base = crtbase(10)

crtbase <- function(Lind,
                    Base=NULL) {
  ##Check parameter consistency  
  if(is.vector(Lind)) LenL <- NROW(Lind) else stop("Lind is not a vector") 
  if(is.null(Base)) {
    Base = 2 * rep(1,LenL) #Default to Base 2
  } else if(!is.vector(Base)) stop("Base is not a vector")
  LenB <- NROW(Base)
  if((LenL > 1 && LenB > 1 && LenL != LenB) || (LenL == 1 && LenB > 1)) stop("Vector dimensions must agree")
  if(LenB == 1 & LenL > 1){
    Base <- Base * rep(1,LenL)
  }
  
  ##main body
  BaseVec <- vector()
  for (i in 1:LenL){
    BaseVec <- c(BaseVec,Base[i]*rep(1,Lind[i]))
  }
  return(BaseVec)
}