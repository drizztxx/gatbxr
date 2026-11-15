# gatbxr 0.2.1

## Bug fixes

* `recombin()` could not find its own operators (e.g. `"xovsp"`) when gatbxr
  was not attached, for example when called as `gatbxr::recombin()` from
  another package. The function names given to `recombin()`, `select()` and
  `mutate()` are now looked up in gatbxr first and then from the calling
  environment, so user-defined operators no longer need to be defined in the
  global environment.

## Package metadata

* The reference in the Description is now written as
  `authors (year) <https:...>`, as requested by CRAN.

# gatbxr 0.2.0

## Package metadata

* Credited the authors and copyright holders of the original Sheffield
  'MATLAB' toolbox in `Authors@R`.
* Changed the license to GPL (>= 2), as required for a derived work of the
  GPL-licensed original toolbox.
* Revised Title and Description following CRAN review.

## Bug fixes

* `xovmp()` and all `xov*()` crossover functions: the crossover probability
  `Px` was inverted (pairs were crossed when they should not be), and the
  single-point cut was always placed after the first locus.
* `recmut()`: the second offspring of each pair now moves in the opposite
  direction (`-Diff`) as in the original toolbox, and identical parent norms
  no longer cause division by zero.
* `bs2rv()`: fixed decoding when `lbin = FALSE`; `lb`, `ub`, `code`, `scale`,
  `lbin` and `ubin` may now be given per variable.
* `crtbp()`: fixed argument handling when `Lind` contains a base vector, and
  added informative errors for unsupported argument combinations.
* `ranking()`: individuals with equal objective values now receive the same
  (average) fitness; `RFun = NULL` now works as documented.
* `migrate()`, `mutate()` and `select()` failed with their default arguments
  on R >= 4.2.0.
* `reins()`: `InsOpt` of length 1 now works as documented.
* `scaling()`: no longer returns `NaN` when all objective values are equal.

## Documentation

* Added examples to all exported functions and fixed spelling throughout.
* Corrected author attributions in several help pages.

# gatbxr 0.1

* First version.
