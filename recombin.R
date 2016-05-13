## RECOMBIN.R       (RECOMBINation high-level function)
##
## This function performs recombination between pairs of individuals
## and returns the new individuals after mating. The function handles
## multiple populations and calls the low-level recombination function
## for the actual recombination process.
##
## Syntax:  NewChrom = recombin(REC_F, OldChrom, RecOpt, SUBPOP)
##
## Input parameters:
##    REC_F     - String containing the name of the recombination or
##                crossover function
##                (NOTE : doesn't work with low level recmut)
##    Chrom     - Matrix containing the chromosomes of the old
##                population. Each line corresponds to one individual
##    RecOpt    - (optional) Scalar containing the probability of 
##                recombination/crossover occurring between pairs
##                of individuals.
##                if omitted or NaN, 0.7 is assumed
##    SUBPOP    - (optional) Number of subpopulations
##                if omitted or NaN, 1 subpopulation is assumed
##    ...       - Ohter aurguments passed on to crossover function
##
## Output:
##    NewChrom  - Matrix containing the chromosomes of the population
##                after recombination in the same format as OldChrom.
##
##  Author:    Hartmut Pohlheim
##             David Zhao (Modified for R)
##
## Date: 13May2016

recombin <- function(REC_F=c("xovmp","xovsp"),
                     Chrom,
                     RecOpt=0.7,
                     SUBPOP=1,
                     ...){
  ## Check parameter consistency
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