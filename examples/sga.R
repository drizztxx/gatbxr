##
## Test objective function is De Jong's first test function.
## Binary representation for the individuals is used.
##
## Author:     Hartmut Pohlheim
## History:    23.03.94     file created
##             15.01.03     tested under MATLAB v6 by Alex Shenfield
##             06.06.16     Re-modified for R by Liang Zhao

NIND = 40;           ## Number of individuals per subpopulations
MAXGEN = 100;        ## max Number of generations
GGAP = 0.9;          ## Generation gap, how many new individuals are created
SEL_F = 'sus';       ## Name of selection function
XOV_F = 'xovsp';     ## Name of recombination function for individuals
MUT_F = 'mut';       ## Name of mutation function for individuals

## Objective function: De Jong's first test function (sphere)
objfun1 <- function(Phen){
  res <- apply(Phen * Phen,1,sum)
  return(res)
}

## set boundaries of objective function
   lb=-512;
   ub=512;

## Number of variables of objective function
   NVAR = 20;

## Build field description list
   PRECI = 20;    ## Precision of binary representation
   FieldDD = list(prec=PRECI,lb=lb,ub=ub,code="gray",scale="arith",lbin=TRUE,ubin=TRUE);

## Create population
   Chrom = crtbp(NIND, NVAR*PRECI)$Chrom;
   ObjV = objfun1(bs2rv(Chrom,FieldDD));

## reset count variables
   gen = 0;
   Best = rep(NA,MAXGEN)

## Iterate population
   while (gen < MAXGEN){
     ## Record the best objective value of the population
     Best[gen+1] = min(ObjV);

     ## Fitness assignment to whole population
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

## Best objective value found in each generation
plot(Best, type = "l", xlab = "Generation", ylab = "Best objective value")
## End of script
