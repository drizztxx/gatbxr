% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crtrp.R
\name{crtrp}
\alias{crtrp}
\title{Create a real-valued initial population}
\usage{
crtrp(Nind,FieldDR)
}
\arguments{
\item{Nind}{a value containing the number of individuals in the new 
population.}

\item{FieldDR}{a matrix of 2 rows by number of variables describing the
boundaries of each variable.}
}
\value{
a matrix containing the random valued individuals of the
new population of size \code{Nind} by number of variables.
}
\description{
This function creates a population of given size of random real-values.
}
\details{
The first step in a genetic algorithm is to create an initial population consisting of 
random individuals. \code{crtrp} produces a matrix, \code{Chrom}, containing uniformly
distributed random values in its elements.

\code{Chrom = crtrp(Nind, FieldDR)} creates a random real-valued matrix of
size \code{Nind * Nvar}, where \code{Nind} specifies the number of individuals in the
population and \code{Nvar} the number of variables of each individual. \code{Nvar} is derived
from \code{FieldDR} with \code{Nvar = NCOL(FieldDR)}.

\code{FieldDR} is a matrix of size 2 * \code{Nvar} and contains
the boundaries of each variable of an individual. The first row contains the lower
bounds, the second row the upper bounds.

\code{FieldDR} is used in other functions \code{\link{mutate}}
}
\examples{
## To create a random populatin of 6 individuals with 4 variables each:
## Define boundaries on the variables,
FieldDR = matrix(c(-100,-50,-30,-20,
                    100,50,30,20),2,4,byrow=TRUE)
## Create initial population
Chrom = crtrp(6,FieldDR)                     
}
\author{
The original matlab implementation of crtrp was written by Hartmut Pohlheim and 
tested by Alex Shenfield. 
The R implementation was written by David Zhao.
}
\seealso{
\code{\link{mutbga}}, \code{\link{recdis}}, \code{\link{recint}}, \code{\link{reclin}}
}

