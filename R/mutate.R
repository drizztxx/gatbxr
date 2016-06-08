## MUTATE.R       (MUTATion high-level function)
##
## This function takes a matrix OldChrom containing the 
## representation of the individuals in the current population,
## mutates the individuals and returns the resulting population.
##
## The function handles multiple populations and calls the low-level
## mutation function for the actual mutation process.
##
## Syntax:  NewChrom = mutate(MUT_F, OldChrom, FieldDR, MutOpt, SUBPOP)
##
## Input parameter:
##    MUT_F     - String containing the name of the mutation function
##    OldChrom  - Matrix containing the chromosomes of the old
##                population. Each line corresponds to one individual.
##    FieldDR   - Matrix describing the boundaries of each variable 
##                (real-values) or defining the base of the variables of 
##                each individual (discrete values).
##                optional for (binary) discrete values
##    MutOpt    - (optional) Vector containing mutation rate and shrink value
##                if omitted or NaN, MutOpt = NaN is assumed
##                MutOpt(1): MutR - number containing the mutation rate -
##                           probability for mutation of a variable
##                MutOpt(2): MutShrink - (optional) number for shrinking the
##                           mutation range in the range [0, 1], possibility to
##                           shrink the range of the mutation depending on,
##                           for instance actual generation (only for
##                           real-values).
##    SUBPOP    - (optional) Number of subpopulations
##                if omitted or NaN, 1 subpopulation is assumed
##    ...       - Ohter aurguments passed on to mutation function
##
## Output parameter:
##    NewChrom  - Matrix containing the chromosomes of the population
##                after mutation in the same format as OldChrom.
##
## Author:     Hartmut Pohlheim
##             David Zhao (Modified for R)
##
## Date: 6Jun2016

mutate <- function(MUT_F=c("mut","mutbga"),
                   OldChrom,
                   FieldDR=NULL,
                   MutOpt=NULL,
                   SUBPOP=1,
                   ...){
  ## Identify the population size (Nind) and the number of variables (Nvar)
  Nind <- NROW(OldChrom); Nvar <- NCOL(OldChrom)
  
  if (is.null(FieldDR)) {
    IsDiscret = 1
  } else {
    mF <- NROW(FieldDR); nF <- NCOL(FieldDR)
    if (nf != Nvar) stop("FieldDR and OldChrom disagree")
    if (mF == 2) {
      IsDiscret = 0
    } else if (mF == 1){
      IsDiscret = 1
    } else stop("FieldDR must be a matrix with 1 or 2 rows")  
  }
  
  if (length(SUBPOP) != 1) stop("SUBPOP must be a scalar")
  
  if (Nind %% SUBPOP) stop("OldChrom and SUBPOP disagree")
  
  Nind <- Nind/SUBPOP ## Compute number of individuals per subpopulation
  
  ## Select individuals of one subpopulation and call low level function
  NewChrom <- NULL
  for(ix in 1:SUBPOP){
    ChromSub <- OldChrom[((ix-1)*Nind+1):(ix*Nind),]
    if (IsDiscret == 1) {
      NewChromSub <- do.call(MUT_F,list(ChromSub,MutOpt,FieldDR,...))
    } else {
      NewChromSub <- do.call(MUT_F,list(ChromSub,FieldDR,MutOpt,...))
    }    
    NewChrom <- rbind(NewChrom,NewChromSub)
  }
  return(NewChrom)
}