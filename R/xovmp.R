#' @title Crossover operators for binary representation
#'
#' @description
#' These functions take a matrix \code{OldChrom} containing the binary
#' representation of the individuals in the current population,
#' apply crossover to consecutive pairs of individuals with
#' defined probability and return the resulting population.
#'
#' @param OldChrom  a matrix containing the chromosomes of the old
#' population. Each row corresponds to one individual.
#' @param Px a number indicating the probability of recombination
#' occurring between pairs of individuals. If omitted, 0.7 is used.
#' @param Npt an optional number indicating the number of crossover points to use.
#' Specifically, 1 indicates single point recombination; 2 indicates
#' double point recombination; and 0 indicates shuffle point recombination.
#' Default is set to 0.
#' @param Rs  an optional logical value indicating whether or not to
#' force the production of offspring different from their parents.
#' Default is set to FALSE.
#' @return
#' a matrix containing the offspring after mating,
#' ready to be mutated and/or evaluated, in the same
#' format as OldChrom.
#' @export
#' @author
#' The original 'MATLAB' implementation of xovmp was written by Carlos Fonseca and
#' updated by Andrew Chipperfield; the other crossover functions were written
#' by Hartmut Pohlheim.
#' The R implementation was written by Liang Zhao.
#' @details
#' \code{xovmp} is the general crossover function for binary representation and
#' is called by all other crossover functions. When \code{xovmp} is called by
#' the high-level function \code{\link{recombin}}, its additional arguments can
#' be passed through \code{...}. With different combinations of arguments, it can
#' be used identically to \code{xovsp}, \code{xovdp}, \code{xovsh}, \code{xovsprs},
#' \code{xovdprs} and \code{xovshrs}. Although these functions are redundant, they
#' are kept for completeness.
#' \cr\cr
#' \code{xovsp} performs single-point crossover between pairs of individuals. When called by
#' \code{\link{recombin}}, \code{recombin("xovsp",Chrom)} is identical to
#' \code{recombin("xovmp",Chrom,Npt=1)}.
#' \cr\cr
#' \code{xovdp} performs double-point crossover between pairs of individuals. When called by
#' \code{\link{recombin}}, \code{recombin("xovdp",Chrom)} is identical to
#' \code{recombin("xovmp",Chrom,Npt=2)}.
#' \cr\cr
#' \code{xovsh} performs shuffle-point crossover between pairs of individuals. When called by
#' \code{\link{recombin}}, \code{recombin("xovsh",Chrom)} is identical to
#' \code{recombin("xovmp",Chrom)}.
#' \cr\cr
#' \code{xovsprs} performs single-point reduced surrogate crossover between
#' pairs of individuals. When called by \code{\link{recombin}},
#' \code{recombin("xovsprs",Chrom)} is identical to
#' \code{recombin("xovmp",Chrom,Npt=1,Rs=TRUE)}.
#' \cr\cr
#' \code{xovdprs} performs double-point reduced surrogate crossover between
#' pairs of individuals. When called by \code{\link{recombin}},
#' \code{recombin("xovdprs",Chrom)} is identical to
#' \code{recombin("xovmp",Chrom,Npt=2,Rs=TRUE)}.
#' \cr\cr
#' \code{xovshrs} performs shuffle-point reduced surrogate crossover between
#' pairs of individuals. When called by \code{\link{recombin}},
#' \code{recombin("xovshrs",Chrom)} is identical to
#' \code{recombin("xovmp",Chrom,Rs=TRUE)}.
#' @seealso
#' \code{\link{recombin}}
#' @examples
#' Chrom = crtbp(40,10)$Chrom
#'
#' ## Shuffle crossover
#' Selch = xovmp(OldChrom=Chrom)
#'
#' ## Single-point crossover
#' Selch = xovsp(Chrom, Px = 0.8)


xovmp <- function(OldChrom,
                  Px=0.7,
                  Npt=0,
                  Rs=FALSE){
  ## Identify the population size (Nind) and the chromosome length (Lind)
  Nind <- NROW(OldChrom); Lind <- NCOL(OldChrom)

  if(Lind < 2 || Nind < 2) return(OldChrom)

  Xops <- floor(Nind/2)
  DoCross <- runif(Xops) < Px
  odd <- seq(1, by = 2, length.out = Xops)
  even <- odd + 1

  ## Compute the effective length of each chromosome pair
  Mask <- (OldChrom[odd,,drop=FALSE] != OldChrom[even,,drop=FALSE]) | !Rs
  Mask <- matrix(t(apply(Mask,1,cumsum)),Xops,Lind)

  ## Compute cross sites for each pair of individuals, according to their
  ## effective length and Px (two equal cross sites mean no crossover)
  xsites <- Mask[,Lind]
  if(Npt >= 2) xsites <- ceiling(xsites * runif(Xops))
  xsites2 <- (xsites + ceiling((Mask[,Lind]-1) * runif(Xops)) * DoCross - 1) %% Mask[,Lind] + 1

  ## Express cross sites in terms of a Boolean mask
  Mask <- (matrix(rep(xsites,Lind),Xops,Lind) < Mask) ==
          (matrix(rep(xsites2,Lind),Xops,Lind) < Mask)
  Mask[is.na(Mask)] <- TRUE  #in case that all elements of individuals in one pair exactly equal

  ##shuffle individuals
  if(!Npt){
    shuff <- matrix(runif(Lind*Xops),Xops,Lind)
    shuff <- matrix(t(apply(shuff,1,order)),Xops,Lind)
    unshuff <- matrix(t(apply(shuff,1,order)),Xops,Lind)
    OldChrom[odd,] <- do.call(rbind,lapply(1:Xops,function(i) OldChrom[odd[i],shuff[i,]]))
    OldChrom[even,] <- do.call(rbind,lapply(1:Xops,function(i) OldChrom[even[i],shuff[i,]]))
  }

  ## Perform crossover
  NewChrom <- matrix(NA,Nind,Lind)
  NewChrom[odd,] <- (OldChrom[odd,,drop=FALSE] * Mask) + (OldChrom[even,,drop=FALSE] * !Mask)
  NewChrom[even,] <- (OldChrom[odd,,drop=FALSE] * !Mask) + (OldChrom[even,,drop=FALSE] * Mask)

  ## If the number of individuals is odd, the last individual cannot be mated
  ## but must be included in the new population
  if(Nind %% 2){
    NewChrom[Nind,] <- OldChrom[Nind,]
  }

  ##unshuffle
  if(!Npt){
    NewChrom[odd,] <- do.call(rbind,lapply(1:Xops,function(i) NewChrom[odd[i],unshuff[i,]]))
    NewChrom[even,] <- do.call(rbind,lapply(1:Xops,function(i) NewChrom[even[i],unshuff[i,]]))
  }
  return(NewChrom)
}

#' @rdname xovmp
#' @export
xovsp <- function(OldChrom,
                  Px=0.7){
  ## call low level function with appropriate parameters
  NewChrom <- xovmp(OldChrom,Px,1,FALSE)
  return(NewChrom)
}

#' @rdname xovmp
#' @export
xovdp <- function(OldChrom,
                  Px=0.7){
  ## call low level function with appropriate parameters
  NewChrom <- xovmp(OldChrom,Px,2,FALSE)
  return(NewChrom)
}

#' @rdname xovmp
#' @export
xovsh <- function(OldChrom,
                  Px=0.7){
  ## call low level function with appropriate parameters
  NewChrom <- xovmp(OldChrom,Px,0,FALSE)
  return(NewChrom)
}

#' @rdname xovmp
#' @export
xovsprs <- function(OldChrom,
                  Px=0.7){
  ## call low level function with appropriate parameters
  NewChrom <- xovmp(OldChrom,Px,1,TRUE)
  return(NewChrom)
}

#' @rdname xovmp
#' @export
xovdprs <- function(OldChrom,
                    Px=0.7){
  ## call low level function with appropriate parameters
  NewChrom <- xovmp(OldChrom,Px,2,TRUE)
  return(NewChrom)
}

#' @rdname xovmp
#' @export
xovshrs <- function(OldChrom,
                    Px=0.7){
  ## call low level function with appropriate parameters
  NewChrom <- xovmp(OldChrom,Px,0,TRUE)
  return(NewChrom)
}