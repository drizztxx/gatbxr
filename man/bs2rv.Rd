% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bs2rv.R
\name{bs2rv}
\alias{bs2rv}
\title{Binary string to real vector}
\usage{
bs2rv(Chrom,FieldD=list(prec=numeric(),lb=numeric(),ub=numeric(),code=c("binary","gray"),
  scale=c("arith","log"),lbin=T,ubin=T))
}
\arguments{
\item{Chrom}{matrix containing the chromosomes of the current
population. Each line corresponds to one
individual's concatenated binary string
representation. Leftmost bits are MSb and
rightmost are LSb.}

\item{FieldD}{a list of additional options describing how to decode
each substring in the chromosome. It has 7 elements,
prec: a scalar describing the precision when converting chromosome to reals;
lb, ub: Lower and upper bounds for each variable;
code: binary row vector indicating how each substring is to be decoded;
scale: binary row vector indicating where to use arithmetic and/or logarithmic scaling;
lbin, ubin: binary row vectors indicating whether or not to include each bound in the
representation range}
}
\value{
real matrix containing the population phenotypes.
}
\description{
This function decodes binary chromosomes into vectors of reals. The
chromosomes are seen as the concatenation of binary strings of given
length, and decoded into real numbers in a specified interval using
either standard binary or Gray decoding.
}
\examples{

## Create a chromsome 
Chrom = crtbp(Nind = 40, Lind = 400)$Chrom

## Converting to reals range from [-512, 512] with a precision of 20 
##(20 bits mean a real number) in  arithmetic scale by binary decoding.
Phen = bs2rv(Chrom = Chrom,
             FieldD = list(prec=20,lb=-512,ub=512,code="binary",scale="arith",lbin=TRUE,ubin=TRUE))
}
\author{
The original matlab implementation of bs2rv was written by Carlos Fonseca and 
updated by Andrew Chipperfield. 
The R implementation was written by David Zhao.
}
\seealso{
\code{\link{crtbp}}
}

