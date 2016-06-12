#' @title MUTATion high-level function
#'
#' @description
#' This function takes a matrix OldChrom containing the 
#' representation of the individuals in the current population,
#' mutates the individuals and returns the resulting population.
#'
#' The function handles multiple populations and calls the low-level
#' mutation function for the actual mutation process.
#'
#' @usage 
#' mutate(MUT_F=c("mut","mutbga"),OldChrom,FieldDR=NULL,MutOpt=NULL,SUBPOP=1,...)
#'
#' @param MUT_F character string containing the name of the mutation function.
#' @param OldChrom  matrix containing the chromosomes of the old
#' population. Each line corresponds to one individual.
#' @param FieldDR matrix describing the boundaries of each variable 
#' (real-values) or defining the base of the variables of 
#' each individual (discrete values). Optional for (binary) discrete values
#' @param MutOpt  vector containing mutation rate and shrink value.
#' MutOpt[1] is the number containing the mutation rate -
#' probability for mutation of a variable.
#' MutOpt[2] is the number for shrinking the
#' mutation range in the range [0, 1], possibility to
#' shrink the range of the mutation depending on,
#' for instance actual generation.
#' @param SUBPOP an optional number of subpopulations.
#' Default is set to 1.
#' @param ... ohter aurguments passed on to mutation function.
#'
#' @return 
#' a matrix containing the chromosomes of the population
#' after mutation in the same format as OldChrom.
#' @export
#' @author 
#' The original matlab implementation of mutate was written by Hartmut Pohlheim.
#' The R implementation was written by David Zhao. 
#' 
#' @examples
#' 
#' Chrom = crtbp(40,10)$Chrom
#' 
#' NewChrom = mutate(MUT_F="mut",OldChrom=Chrom)

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