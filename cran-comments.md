## Resubmission

This is a resubmission of gatbxr, which was first submitted as version 0.1 in
August 2016. The comments from that review (Uwe Ligges, 2016-08-20/21) have
been addressed as follows:

* "Please omit 'Implemented in R' which is redundant" / "write at least
  'Implemented in R' (rather than 'by')": the Title is now
  "Port of the Sheffield Genetic Algorithm Toolbox for 'MATLAB'".

* "Please omit 'R version of'": the Description no longer starts with
  "R version of" and now describes the functionality of the package.

* "Please also indicate who is author / copyright holder. Since it is derived
  work, the other people may still be authors / copyright holders":
  the authors of the original 'MATLAB' toolbox (Andrew Chipperfield, Peter
  Fleming, Hartmut Pohlheim, Carlos Fonseca) are now listed in Authors@R with
  roles "aut" and "cph", the University of Sheffield as "cph", Alex Shenfield
  (tester of the original toolbox) as "ctb", and the maintainer as
  "aut", "cre" and "trl".

* The original toolbox is distributed under the GNU General Public License,
  so the license of this derived work has been changed from MIT to
  GPL (>= 2).

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* local Windows 11, R 4.5.1
* win-builder (devel and release)
* GitHub Actions: macOS, Windows and Ubuntu (R devel, release and oldrel-1)
