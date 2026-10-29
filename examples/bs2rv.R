## Create a population of 40 chromosomes of length 400
Chrom = crtbp(Nind = 40, Lind = 400)$Chrom

## Convert to reals in the range [-512, 512] with a precision of 20
## (20 bits represent one real number) on an arithmetic scale by binary decoding.
Phen = bs2rv(Chrom = Chrom,
             FieldD = list(prec=20,lb=-512,ub=512,code="binary",scale="arith",lbin=TRUE,ubin=TRUE))

## Two variables with different bounds, decoded with Gray code
Chrom = crtbp(5, 16)$Chrom
Phen = bs2rv(Chrom,
             FieldD = list(prec=8,lb=c(-1,10),ub=c(1,100),code="gray",
                           scale=c("arith","log"),lbin=TRUE,ubin=TRUE))
