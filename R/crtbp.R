#' @title Create an initial population
#'
#' @description
#' This function creates a binary population of given size and structure.
#'
#' @param Nind	either a number containing the number of individuals
#' in the new population or a vector of length two
#' containing the number of individuals and their length.
#' @param Lind	an optional number containing the length of the individual
#' chromosomes, or a vector containing the base of each locus of the
#' chromosomes.
#' @param Base	an optional number containing the base of the chromosome
#' elements or a vector containing the base(s)
#' of the loci of the chromosomes.
#'
#' @return
#' a list containing following components:
#' \item{Chrom}{a matrix containing the random valued chromosomes
#' row wise.}
#' \item{Lind}{a scalar containing the length of the chromosome.}
#' \item{BaseV}{a vector containing the base of the
#' chromosome loci.}
#' @export
#' @author
#' The original 'MATLAB' implementation of crtbp was written by Andrew Chipperfield.
#' The R implementation was written by Liang Zhao.
#' @seealso \code{\link{crtbase}}, \code{\link{crtrp}}
#' @examples
#' ## 40 individuals with 10 binary loci
#' Chrom = crtbp(40,10)$Chrom
#'
#' ## 6 individuals with 5 loci of base 8
#' Chrom = crtbp(6,5,8)$Chrom
#'
#' ## 6 individuals whose loci have bases 2, 4 and 8
#' Chrom = crtbp(6,c(2,4,8))$Chrom

crtbp <- function(Nind,
                  Lind=NULL,
                  Base=NULL){
  ## Check parameter consistency
  nN <- length(Nind)

  if(nN == 2){
    if(is.null(Lind) && is.null(Base)) {
      Lind <- Nind[2]
      Nind <- Nind[1]
      BaseV <- crtbase(Lind)
    } else if(is.null(Base) && length(Lind) == 1){
      BaseV <- crtbase(Nind[2],Lind)
      Lind <- Nind[2]
      Nind <- Nind[1]
    } else if(is.null(Base)){
      if(Nind[2] != length(Lind)) stop("Lind and Base disagree")
      BaseV <- Lind
      Lind <- Nind[2]
      Nind <- Nind[1]
    } else stop("Base must not be given when Nind has length 2")
  } else if(nN == 1){
    if(is.null(Lind)) stop("Lind is required when Nind is a scalar")
    if(is.null(Base)){
      if(length(Lind) == 1) BaseV <- crtbase(Lind) else {
        BaseV <- Lind
        Lind <- length(Lind)
      }
    } else {
      if(length(Lind) != 1) stop("Lind must be a scalar when Base is given")
      if(length(Base) == 1) {
        BaseV <- crtbase(Lind,Base)
      } else if(length(Base) != Lind){
        stop("Lind and Base disagree")
      } else BaseV <- Base
    }
  } else stop("Input parameters inconsistent")

  ## Create a structure of random chromosomes in row wise order, dimensions
  ## Nind by Lind. The base of each chromosomes loci is given by the value
  ## of the corresponding element of the row vector base.
  rand.matrix <- matrix(runif(Nind * Lind),Nind,Lind)
  Base.matrix <- matrix(BaseV,Nind,Lind,byrow = TRUE)
  Chrom <- floor(rand.matrix * Base.matrix)
  result <- list(Chrom=Chrom,Lind=Lind,BaseV=BaseV)
  return(result)
}
