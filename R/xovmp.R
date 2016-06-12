#' @title Multi-point crossover
#'
#' @description
#' This function takes a matrix OldChrom containing the binary
#' representation of the individuals in the current population,
#' applies crossover to consecutive pairs of individuals with
#' defined probability and returns the resulting population.
#'
#' @usage
#' xovmp(OldChrom,Px=0.7,Npt=0,Rs=TRUE)
#' 
#' @param OldChrom  a matrix containing the chromosomes of the old
#' population. Each row corresponds to one individual.
#' @param Px a number indicating the probability of recombination 
#' ocurring between pairs of individuals. If omitted, 0.7 is set.
#' @param Npt an optinal number indicating the number of crossover points to use.
#' Specifically, 1 indicates single point recombination; 2 indicates 
#' double point recombination; and 0 indicates shuffle point recombination.
#' Default is set to 0.
#' @param Rs  an optional logical value indicating whether or not to 
#' force the production of offspring different from their parents.
#' Default is set to TRUE.
#' @return
#' a matrix containing the offspring after mating, 
#' ready to be mutated and/or evaluated, in the same 
#' format as OldChrom.
#' @export
#' @author 
#' The original matlab implementation of mutate was written by Carlos Fonseca and
#' updated by Andrew Chipperfield.
#' The R implementation was written by David Zhao. 
#' @examples
#' 
#' Chrom = crtbp(40,10)$Chrom
#' 
#' Selch = xovmp(OldChrom=Chrom)


xovmp <- function(OldChrom,
                  Px=0.7,
                  Npt=0,
                  Rs=TRUE){
  ## Identify the population size (Nind) and the chromosome length (Lind)
  Nind <- NROW(OldChrom); Lind <- NCOL(OldChrom)
  
  if(Lind < 2) return(OldChrom)
  
  Xops <- floor(Nind/2)
  DoCross <- runif(Xops) < Px
  odd <- (1:(Nind-1))[as.logical(1:(Nind-1) %% 2)]
  even <- (1:Nind)[!as.logical(1:Nind %% 2)]
  
  ## Compute the effective length of each chromosome pair
  Mask <- outer(OldChrom[odd,] != OldChrom[even,],Rs,"|")
  if(Xops==1) Mask <- t(Mask)
  Mask <- t(apply(Mask,1,cumsum))
  
  ## Compute cross sites for each pair of individuals, according to their
  ## effective length and Px (two equal cross sites mean no crossover)
  xsites <- Mask[,Lind]
  if(Npt >= 2) xsites <- ceiling(xsites * runif(Xops))
  xsites2 <- (xsites + ceiling((Mask[,Lind]-1) * runif(Xops)) * (DoCross-1)) %% Mask[,Lind] + 1
  
  ## Express cross sites in terms of a bull mask
  Mask <- (matrix(rep(xsites,Lind),Xops,Lind) < Mask) == 
          (matrix(rep(xsites2,Lind),Xops,Lind) < Mask)
  Mask[is.na(Mask)] <- TRUE  #in case that all elements of individuals in one pair exactly equal
  
  ##shuffle individuals
  if(!Npt){
    shuff <- matrix(runif(Lind*Xops),Xops,Lind)
    shuff <- t(apply(shuff,1,order))
    unshuff <- t(apply(shuff,1,order))
    OldChrom[odd,] <- do.call(rbind,lapply(1:Xops,function(i) OldChrom[odd[i],shuff[i,]]))
    OldChrom[even,] <- do.call(rbind,lapply(1:Xops,function(i) OldChrom[even[i],shuff[i,]]))
  }
  
  ## Perform crossover
  NewChrom <- matrix(NA,Nind,Lind)
  NewChrom[odd,] <- (OldChrom[odd,] * Mask) + (OldChrom[even,] * !Mask)
  NewChrom[even,] <- (OldChrom[odd,] * !Mask) + (OldChrom[even,] * Mask)
  
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