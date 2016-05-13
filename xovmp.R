## XOVMP.R                Multi-point crossover
##
## This function takes a matrix OldChrom containing the binary
## representation of the individuals in the current population,
## applies crossover to consecutive pairs of individuals with
## probability Px and returns the resulting population.
##
## Npt indicates how many crossover points to use (1 or 2, zero
## indicates shuffle crossover).
## Rs indicates whether or not to force the production of
## offspring different from their parents.
##
## Syntax: NewChrom =  xovmp(OldChrom, Px, Npt, Rs)
## 
## Input parameters:
##   OldChrom  - Matrix containing the chromosomes of the old
##               population. Each row corresponds to one individual
##               (in binary form).
##   Px        - Probability of recombination ocurring between pairs
##               of individuals.
##   Npt       - Scalar indicating the number of crossover points
##               1: single point recombination
##               2: double point recombination
##               0: shuffle point recombination
##   Rs        - reduced surrogate
##               False: no reduced surrogate
##               True : reduced surrogate
## 
## Output:
##   NewChrom  - Matrix containing the offspring after mating, 
##               ready to be mutated and/or evaluated, in the same 
##               format as OldChrom.
##
## Author: Carlos Fonseca, 	Updated: Andrew Chipperfield
##         David Zhao (Modified for R)
##
## Date: 13May2016

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