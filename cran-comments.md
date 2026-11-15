## Update

This is a minor update (0.2.1) of gatbxr 0.2.0:

* Bug fix: `recombin()` could not find its own operators when gatbxr was not
  attached (e.g. when called as `gatbxr::recombin()` from another package).
  Operator names in `recombin()`, `select()` and `mutate()` are now looked up
  in gatbxr first and then from the calling environment.

* As requested by CRAN on acceptance of 0.2.0, the reference in the
  Description is now written as `authors (year) <https:...>`.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are currently no reverse dependencies on CRAN.

## Test environments

* local Windows 11, R 4.5.1
* win-builder (devel and release)
* GitHub Actions: macOS, Windows and Ubuntu (R devel, release and oldrel-1)
