## SGA.M          (Simple Genetic Algorithm)
##
## This script implements the Simple Genetic Algorithm.
## Test objective function is De Jong's first test funciton.
## Binary representation for the individuals is used.
##
## Author:     Hartmut Pohlheim
## History:    23.03.94     file created
##             15.01.03     tested under MATLAB v6 by Alex Shenfield
##             06.06.16     Re-modified for R by David Zhao

require("gatbxr")

NIND = 40;           ## Number of individuals per subpopulations
MAXGEN = 300;        ## max Number of generations
GGAP = 0.9;          ## Generation gap, how many new individuals are created
SEL_F = 'sus';       ## Name of selection function
XOV_F = 'xovsp';     ## Name of recombination function for individuals
MUT_F = 'mut';       ## Name of mutation function for individuals
OBJ_F = 'objfun1';   ## Name of function for objective values

objfun1 <- function(Chrom){
  res <- apply(Chrom * Chrom,1,sum)
  return(res)
}

## set boundaries of objective function
   lb=-512;
   ub=512;

## Number of variables of objective function, in OBJ_F defined
   NVAR = 20;   

## Build fielddescription matrix
   PRECI = 20;    ## Precisicion of binary representation
   FieldDD = list(prec=PRECI,lb=lb,ub=ub,code="gray",scale="arith",lbin=TRUE,ubin=TRUE);

## Create population
   Chrom = crtbp(NIND, NVAR*PRECI)$Chrom;
   ObjV = objfun1(bs2rv(Chrom,FieldDD));

## reset count variables
   gen = 0;
   Best = rep(NA,MAXGEN)

## Iterate population
   while (gen < MAXGEN){
     ## Calculate objective function for population
     Best[gen+1] = min(ObjV);
     
     ## Fitness assignement to whole population
     FitnV = ranking(ObjV);
     
     ## Select individuals from population
     SelCh = select(SEL_F, Chrom, FitnV, GGAP);
     
     ## Recombine selected individuals (crossover)
     SelCh=recombin(XOV_F, SelCh);
     
     ## Mutate offspring
     SelCh=mutate(MUT_F, SelCh);
     
     ##Evaluate offspring, call objective function
     ObjVSel = objfun1(bs2rv(SelCh,FieldDD));
     
     ## Insert offspring in population replacing parents
     rs = reins(Chrom, SelCh,1,c(1,1),ObjV,ObjVSel);
     Chrom = rs$Chrom;
     ObjV  = rs$ObjVCh

     gen = gen + 1;
   }
## End of script
