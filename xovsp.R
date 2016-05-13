## XOVSP.R        (CROSSOVer Single-Point)
##
## This function performs single-point crossover between pairs of 
## individuals and returns the current generation after mating.
##
## Syntax:  NewChrom = xovsp(OldChrom, XOVR)
##
## Input parameters:
##    OldChrom  - Matrix containing the chromosomes of the old
##                population. Each line corresponds to one individual
##                (in any form, not necessarily real values).
##    XOVR      - Probability of recombination occurring between pairs
##                of individuals.
##
## Output:
##    NewChrom  - Matrix containing the chromosomes of the population
##                after mating, ready to be mutated and/or evaluated,
##                in the same format as OldChrom.
##
##  Author:    Hartmut Pohlheim
##             David Zhao (Modified for R)
##
## Date: 13May2016

xovsp <- function(OldChrom,
                  XOVR=0.7){
  ## call low level function with appropriate parameters
  NewChrom <- xovmp(OldChrom,XOVR,1,TRUE)
  return(NewChrom)
}