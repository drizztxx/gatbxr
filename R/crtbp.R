#' @title Create an initial population
#'
#' @description
#' This function creates a binary population of given size and structure.
#'
#' @param Nind	- either a number containing the number of individuals
#' in the new population or a row vector of length two
#' containing the number of individuals and their length.
#' @param Lind	an optional number containing the length of the individual
#' chromosomes.
#' @param Base	an optional number containing the base of the chromosome 
#' elements or a row vector containing the base(s) 
#' of the loci of the chromosomes.
#'
#' @return
#' a list containing following components:
#' \item{Chrom}{a matrix containing the random valued chromosomes 
#' row wise.}
#' \item{Lind}{a scalar containing the length of the chromosome.}
#' \item{BaseV}{a row vector containing the base of the 
#' chromosome loci.}
#' @export
#' @author 
#' The original matlab implementation of mutate was written by Andrew Chipperfield.
#' The R implementation was written by David Zhao. 
#' @examples
#' 
#' Chrom = crtbp(40,10)

crtbp <- function(Nind,
                  Lind=NULL,
                  Base=NULL){
  ## Check parameter consistency 
  nN <- NROW(Nind)
  if(!is.null(Lind)) nL <- NROW(Nind)
  if(!is.null(Base)) nB <- NROW(Base)
  
  if(nN == 2){
    if(is.null(Lind) && is.null(Base)) {
      Lind <- Nind[2]
      Nind <- Nind[1]
      BaseV <- crtbase(Lind)
    } else if(is.null(Base) && !is.null(Lind) && nL == 1){
      BaseV <- crtbase(Nind[2],Lind)
      Lind <- Nind[2]
      Nind <- Nind[1]
    } else if(is.null(Base) && !is.null(Lind) && nL > 1){
      if(!length(Lind) %in% Lind) stop("Lind and Base disagree")
      BaseV <- Lind
      Lind <- Nind[2]
      Nind <- Nind[1]
    }
  } else if(nN == 1){
    if(is.null(Base) && !is.null(Lind)){
      if(nL == 1) BaseV <- crtbase(Lind) else {
        BaseV <- Lind
        Lind <- nL
      }
    } else if(!is.null(Base) && !is.null(Lind)){
      if(nB == 1) {
        BaseV <- crtbase(Lind,Base)
      } else if(nB != Lind){
        stop("Lind and Base disagree")
      } else BaseV <- Base
    }
  } else stop("Input parameters inconsistent")
  if(length(Nind) != 1 || length(Lind) != 1) stop("Input parameters inconsistent")
  if(Lind != length(BaseV)) stop("Lind and Base disagree")
  
  ## Create a structure of random chromosomes in row wise order, dimensions
  ## Nind by Lind. The base of each chromosomes loci is given by the value
  ## of the corresponding element of the row vector base.
  rand.matrix <- matrix(runif(Nind * Lind),Nind,Lind)
  Base.matrix <- matrix(rep(BaseV,Nind),Nind,Lind,byrow = T)
  Chrom <- floor(unlist(rand.matrix) * unlist(Base.matrix))
  result <- list(Chrom=Chrom,Lind=Lind,BaseV=BaseV)
  return(result)
}