## BS2RV.R - Binary string to real vector
##
## This function decodes binary chromosomes into vectors of reals. The
## chromosomes are seen as the concatenation of binary strings of given
## length, and decoded into real numbers in a specified interval using
## either standard binary or Gray decoding.
##
## Syntax:       Phen = bs2rv(Chrom,FieldD)
##
## Input parameters:
##
##               Chrom    - Matrix containing the chromosomes of the current
##                          population. Each line corresponds to one
##                          individual's concatenated binary string
##           			   representation. Leftmost bits are MSb and
##			               rightmost are LSb.
##
##               FieldD   - List describing how to decode
##           			   each substring in the chromosome. It has the
##           			   following structure:
##
##				[prec;	(num)
##				 lb;		(num)
##				 ub;		(num)
##				 code;		(binary        | gray)
##				 scale;		(arith         | log)
##				 lbin;		(False=excluded| True=included)
##				 ubin];		(False=excluded| True=included)
##
##			   where
##				prec - a scalar describing the precision when 
##               converting chromosome to reals.
##				lb,
##				ub    - Lower and upper bounds for each
##			    		variable. 
##				code  - binary row vector indicating how each
##				    	substring is to be decoded.
##				scale - binary row vector indicating where to
##				    	use arithmetic and/or logarithmic
##				    	scaling.
##				lbin,
##				ubin  - binary row vectors indicating whether
##				    	or not to include each bound in the
##				    	representation range
##
## Output parameter:
##
##               Phen     - Real matrix containing the population phenotypes.
##
## Author: Carlos Fonseca, 	Updated: Andrew Chipperfield,   
##         David Zhao (Modified for R)
##
## Date: 7Jun2016

bs2rv <- function(Chrom,
                  FieldD=list(prec=numeric(),
                              lb=numeric(),
                              ub=numeric(),
                              code=c("binary","gray"),
                              scale=c("arith","log"),
                              lbin=T,
                              ubin=T)){
  ## Identify the population size (Nind) and the chromosome length (Lind) 
  Nind <- NROW(Chrom); Lind <- NCOL(Chrom)
  
  ## Identify the number of decision variables (Nvar)
  if (length(FieldD) != 7) stop("FieldD must have 7 elements")
  
  ## Get substring properties
  prec  <- FieldD$prec
  lb    <- FieldD$lb
  ub    <- FieldD$ub
  code  <- FieldD$code
  scale <- FieldD$scale
  lin   <- FieldD$lbin
  uin   <- FieldD$ubin

  ## Check substring properties for consistency
  if (Lind %% prec) stop("prec disagree with chromosome length")
  Nvar <- Lind/prec
  if (scale == "log" && (lb * ub <= 0)){
    stop("Log-scaled variables must not include 0 in their range")
  }
  
  ## Decode chromosomes
  Phen <- matrix(0,Nind,Nvar)
  
  Prec <- 0.5^prec
  if (scale == "log"){
    lb <- log(abs(lb))
    ub <- log(abs(ub))
  }
  delta <- ub - lb
  
  num <- !lin * Prec
  den <- (lin + uin - 1) * Prec
  
  for(i in 1:Nvar){
    idx <- ((i-1)*prec+1):(i*prec)
    if (code == "gray"){ ## Gray decoding
      Chrom[,idx] <- t(apply(Chrom[,idx],1,cumsum)) %% 2
    }
    Phen[,i] <- Chrom[,idx] %*% (0.5 ^ (1:prec))
    Phen[,i] <- lb + delta * (Phen[,i] + num) / (1 - den)
  }
  if (scale == "log"){
    if (FieldD$lb < 0 ) {
      Phen = -1 * exp(Phen)
    } else {
      Phen = exp(Phen)
    }
  }
  return(Phen)
}