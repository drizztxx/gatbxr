## Create a chromsome 
Chrom = crtbp(Nind = 40, Lind = 400)$Chrom

## Converting to reals range from [-512, 512] with a precision of 20 (20 bits mean a real) in 
## arithmetic scale by binary decode.
Phen = bs2rv(Chrom = Chrom,
             FieldD = list(prec=20,lb=-512,ub=512,code="binary",scale="arith",lbin=T,ubin=T))