#' @title Universal selection high-level function
#'
#' @description
#' This function performs universal selection. The function handles
#' multiple populations and calls the low level selection function
#' for the actual selection process.
#'
#' @usage
#' select(SEL_F=c("sus","rws"),Chrom,FitnV,GGAP=1,SUBPOP=1,...)
#'
#' @param SEL_F character string indicating the selection function.
#' @param Chrom a matrix containing the individuals (parents) of the current
#' population. Each row corresponds to one individual.
#' @param FitnV a vector containing the fitness values of the
#' individuals in the population.
#' @param GGAP  an optional number indicating rate of individuals to be selected.
#' Default is set to 1.
#' @param SUBPOP  an optional number indicating subpopulations.
#' Default is set to 1 subpopulation.
#' @param ... ohter aurguments passed on to selection function.
#'
#' @return 
#' a matrix containing the selected individuals.
#' @export
#' @author 
#' The original matlab implementation of mutate was written by Hartmut Pohlheim.
#' The R implementation was written by David Zhao. 
#' @examples
#' 
#' Chrom = crtbp(40,10)$Chrom
#' 
#' #Ojbective function is the sum of individuals
#' FitnV = ranking(apply(Chrom,1,sum))
#' 
#' Selch = select(SEL_F="sus",Chrom=Chrom,FitnV=FitnV,GGAP=0.9)

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