# gatbxr

[![R-CMD-check](https://github.com/drizztxx/gatbxr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/drizztxx/gatbxr/actions/workflows/R-CMD-check.yaml)

Genetic Algorithm Toolbox for R, ported from the Sheffield Genetic Algorithm
Toolbox for MATLAB.

gatbxr provides operators for creating populations, fitness assignment,
selection, crossover, recombination, mutation, reinsertion and migration
between subpopulations for binary, integer and real-valued representations.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("drizztxx/gatbxr")
```

## Example

```r
library(gatbxr)

## Minimise the sphere function in 10 dimensions
objfun <- function(Phen) rowSums(Phen^2)
FieldD <- list(prec = 20, lb = -512, ub = 512, code = "gray",
               scale = "arith", lbin = TRUE, ubin = TRUE)

Chrom <- crtbp(40, 10 * 20)$Chrom
ObjV <- objfun(bs2rv(Chrom, FieldD))
for (gen in 1:100) {
  FitnV <- ranking(ObjV)
  SelCh <- select("sus", Chrom, FitnV, GGAP = 0.9)
  SelCh <- recombin("xovsp", SelCh)
  SelCh <- mutate("mut", SelCh)
  ObjVSel <- objfun(bs2rv(SelCh, FieldD))
  res <- reins(Chrom, SelCh, 1, c(1, 1), ObjV, ObjVSel)
  Chrom <- res$Chrom
  ObjV <- res$ObjVCh
}
min(ObjV)
```

## Origin

The [Genetic Algorithm Toolbox for MATLAB](https://github.com/UoS-CODeM/GA-Toolbox)
was developed at the Department of Automatic Control and Systems Engineering of
The University of Sheffield, UK, in order to make GAs accessible to the control
engineer within the framework of an existing computer-aided control system design
package. The toolbox was written with the support of a UK SERC grant, and the
final version (v1.2) was completed in 1994. It was written by Andrew Chipperfield,
Peter Fleming, Hartmut Pohlheim and Carlos Fonseca.

The original toolbox is copyright the original authors and The University of
Sheffield, and is published under the GNU General Public License. As a derived
work, gatbxr is also distributed under the GNU General Public License
(version 2 or later).
