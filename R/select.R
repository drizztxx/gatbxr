## SELECT.R          (universal SELECTion)
##
## This function performs universal selection. The function handles
## multiple populations and calls the low level selection function
## for the actual selection process.
##
## Syntax:  SelCh = select(SEL_F, Chrom, FitnV, GGAP, SUBPOP)
##
## Input parameters:
##    SEL_F     - Name of the selection function
##    Chrom     - Matrix containing the individuals (parents) of the current
##                population. Each row corresponds to one individual.
##    FitnV     - Column vector containing the fitness values of the
##                individuals in the population.
##    GGAP      - (optional) Rate of individuals to be selected
##                if omitted 1.0 is assumed
##    SUBPOP    - (optional) Number of subpopulations
##                if omitted 1 subpopulation is assumed
##    ...       - Ohter aurguments passed on to selection function
##
## Output:
##              - Matrix containing the selected individuals.
##
## Author:     Hartmut Pohlheim
##             David Zhao (Modified for R)
##
## Date: 12May2016

select <- function(SEL_F=c("sus","rws"),
                   Chrom,
                   FitnV,
                   GGAP=1,
                   SUBPOP=1,
                   ...){
  ## Check parameter consistency
  NindCh <- NROW(Chrom)
  Nvar <- NCOL(Chrom)
  if(!is.vector(FitnV)) stop("FitnV must be a vector")
  NindF <- NROW(FitnV)
  if(NindCh != NindF) stop("Chrom and FitnV disagree")
  
  if(length(SUBPOP) != 1) stop("SUBPOP must be a scalar")
  if(NindCh %% SUBPOP != 0) stop("Chrom and SUBPOP disagree")
  Nind = NindCh/SUBPOP ## Compute number of individuals per subpopulation
  
  if(length(GGAP) != 1) stop("GGAP must be a scalar")
  if(GGAP < 0) stop("GGAP must be a scalar bigger than 0")
  NSel <- max(floor(Nind*GGAP+0.5),2) ## Compute number of new individuals (to select)
  
  ## Select individuals from population
  SelCh <- NULL
  for(ix in 1:SUBPOP){
    FitnVSub <- FitnV[((ix-1)*Nind+1):(ix*Nind)]
    ChrIx <- do.call(SEL_F,list(FitnVSub,NSel,...)) + (ix-1)*Nind
    SelChSub <- Chrom[ChrIx,]
    SelCh <- rbind(SelCh,SelChSub)
  }
  
  return(SelCh)
}