#' @title Migration of individuals between subpopulations
#' @description
#' This function performs migration of individuals.
#' @usage
#' migrate(Chrom, SUBPOP = 1, ObjV, MIGR = 0.2, Select = FALSE, 
#' Structure = c("net","neighbourhood","ring"))
#' @param Chrom a matrix containing the individuals of the current
#' population. Each row corresponds to one individual.
#' @inheritParams ranking
#' @param MIGR a number indicating rate of individuals to be 
#' migrated per subpopulation. Default is set to 0.2.
#' @param Select an logical value indicating whether fitness-based selection
#' will be used. Default is set to FALSE. If set to TRUE, then \code{ObjV} must also 
#' be defined.
#' @param Structure an character string indicating the name of structure of the subpopulations 
#' for migration. The default structure "net" indicates net structure (unconstrained 
#' migration).
#' @details
#' \pkg{gatbxr} provides support for multiple subpopulations through the use of 
#' high-level genetic operator functions and a routine for exchanging individuals
#' between subpopulations. In the literature, the use of multiple populations has
#' been shown, in most cases, to improve the quality of the results obtained using GAs
#' compared to the single population GA.
#' \cr\cr
#' \pkg{gatbxr} supports the use of a single population divided into a number of 
#' subpopulations or \eqn{demes} by modifying the use of data structures such that 
#' subpopulations are stored in contiguous blocks wihtin a single matrix.
#' \cr\cr
#' This model enables each subpopulation evolved over generations by a traditional GA and 
#' from time to time individuals migrate from one subpopulation to another. The amount of 
#' migration of individuals and the pattern of that migration determines how much genetic
#' diversity can occur.
#' \cr\cr
#' \code{migrate} implements the transfer of individuals between subpopulations, in which
#' \code{Structure} specifies the population topology over which migration will tha place.
#' The default structure unrestricted migration topology (\code{net}) is the most general 
#' migration strategy. Here, individuals may migrate from any subpopulation to another. For each 
#' subpopulation, a pool of potential immigrants is constructed from the other subpopulations.
#' The individual migrants are then determined according to the appropriate selection
#' strategy.
#' \cr\cr
#' Structure \code{ring} transfers individuals between directionally adjacent subpopulations.
#' And \code{neighbourhood} is a similar strategy. Like the \code{ring} topology, migration
#' is made only between nearest neighbours, however, migration may occur in either direction
#' between subpopulations. For each subpopulation, the possible immigrants are determined,
#' according to the desired selection method, from the adjacent subpopulations and a final
#' selection made from this pool of individuals. This ensures that individuals will not 
#' migrate from a subpopulation to the same population.
#' @return a list containing following components:
#' \item{Chrom}{a matrix containing the individuals of the current
#'  population after reinsertion.}
#' \item{ObjV}{if ObjV is input parameters, then return 
#'  vector containing the objective values of the individuals
#'  of the current generation after reinsertion. Else return NULL.}
#' @seealso
#' \code{\link{select}}, \code{\link{recombin}}, \code{\link{mutate}},
#' \code{\link{reins}}     
#' @export
#' @author 
#' The original matlab implementation of reins was written by Hartmut Pohlheim.
#' The R implementation was written by David Zhao.  
#' @examples
#' ## create initial population
#' Chrom = crtbp(50,10)$Chrom
#' 
#' ## calculate objective value with sum function
#' objv = apply(Chrom,1,sum)
#' 
#' ## assume there are 5 subpopulations, then perform ranking
#' ObjV = ranking(objv, SUBPOP = 5)
#' 
#' ## chosse 20% of the individuals of one subpopulation and replaces
#' ## these individuals with the fittest individuals from an adjacent
#' ## subpopulation in a unidirectional ring structure
#' res = migrate(Chrom,SUBPOP = 5, ObjV, Select = TRUE, Structure = "ring")

migrate <- function(Chrom, SUBPOP = 1, ObjV, MIGR = 0.2, Select = FALSE, 
                    Structure = c("net","neighbourhood","ring")){
  ## Check parameter consistency
  Nind <- NROW(Chrom); Nvar <- NCOL(Chrom)
  if (length(SUBPOP) != 1) {
    stop("SUBPOP must be a scalar")
  } 
  if(SUBPOP == 1 || MIGR == 0){
    if (missing(ObjV)) ObjV <- NULL
    return(list(Chrom = Chrom, ObjV = ObjV))
  }
  if (Nind %% SUBPOP) stop("Chrom and SUBPOP disagree")
  NIND <- Nind/SUBPOP
  
  if (!missing(ObjV)){
    if (!is.vector(ObjV)) stop("ObjV must be a vector")
    mO <- NROW(ObjV)
    if (Nind != mO) stop("Chrom and ObjV disagree")
  }
  
  if (MIGR < 0 || MIGR > 1) {
    MIGR <- 0.2
    warning("Parameter for migration rate must be a scalar in [0, 1] and is temporarily
            set to default 0.2")
  }
  
  if (!Structure %in% c("net","neighbourhood","ring")){
    Structure <- "net"
    warning("Strucutre name illeagal, and this value is temporarily set to
            the default net")
  }
  
  if (Select && missing(ObjV)) {
    Select <- FALSE
    warning("ObjV for fitness-based migration needed, Select has temporarily set to FALSE")
  }
  
  MigTeil <- max(floor(NIND * MIGR), 1) ## number of individuals to migrate
  
  ## Perform migration between subpopulations --> create a matrix for migration
  ## in every subpopulation from best individuals of the other subpopulations
  ChromMigAll <- NULL
  if (!missing(ObjV)) ObjVAll <- NULL
  for(x in 1:SUBPOP){
    ## sort ObjV of actual subpopulation
    if (Select) {
      IndMigSo <- order(ObjV[((x-1)*NIND+1):(x*NIND)])
    } else IndMigSo <- order(runif(NIND))
    ## take MigTeil (best) individuals, copy individuals and objective values
    IndMigTeil <- IndMigSo[1:MigTeil] + (x-1)*NIND
    ChromMigAll <- rbind(ChromMigAll, Chrom[IndMigTeil,])
    if (!missing(ObjV)) ObjVAll <- c(ObjVAll,ObjV[IndMigTeil])
  }
  
  ## perform migration
  for (i in 1:SUBPOP){
    ChromMig <- ChromMigAll
    if (!missing(ObjV)) ObjVMig <- ObjVAll
    if (Structure == "neighbourhood"){
      ## select individuals of neighbourhood subpopulations for ChromMig and ObjVMig
      popnum <- c(SUBPOP,1:SUBPOP,1)
      ins1 <- popnum[i]; ins2 <- popnum[i+2]
      InsRows <- c(((ins1-1)*MigTeil+1):(ins1*MigTeil),
                   ((ins2-1)*MigTeil+1):(ins2*MigTeil))
      ChromMig <- ChromMig[InsRows,]
      if (!missing(ObjV)) ObjVMig <- ObjVMig[InsRows]
    } else if (Structure == "ring"){
      ## select individuals of actual-1 subpopulation for ChromMig and ObjVMig
      popnum <- c(SUBPOP,1:SUBPOP,1)
      ins1 <- popnum[i]
      InsRows <- ((ins1-1)*MigTeil+1):(ins1*MigTeil)
      ChromMig <- ChromMig[InsRows,]
      if (!missing(ObjV)) ObjVMig <- ObjVMig[InsRows]    
    } else {
      ## delete individuals of actual subpopulation from ChromMig and ObjVMig
      DelRows <- ((i-1)*MigTeil+1):(i*MigTeil)
      ChromMig <- ChromMig[-DelRows,]
      if (!missing(ObjV)) ObjVMig <- ObjVMig[-DelRows]
    }
    ## Create an index from a sorted vector with random numbers
    IndMigRa <- order(runif(NROW(ChromMig)))
    ## Take MigTeil numbers from the random vector
    IndMigN <- IndMigRa[1:MigTeil]
    ## copy MigTeil individuals into Chrom and ObjV
    Chrom[(1:MigTeil)+(i-1)*NIND,] <- ChromMig[IndMigN,]
    if (!missing(ObjV)) ObjV[(1:MigTeil)+(i-1)*NIND] <- ObjVMig[IndMigN]
  }
  if (missing(ObjV)) ObjV <- NULL
  return(list(Chrom = Chrom, ObjV = ObjV))
}