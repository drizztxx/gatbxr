test_that("built-in operators are found when gatbxr is not visible from the caller", {
  set.seed(1)
  Chrom <- crtbp(10, 8)$Chrom
  ## The caller sees nothing but these objects, as inside another package that
  ## calls gatbxr::recombin() without importing the operators. (An environment
  ## on the search path would see gatbxr, which the test runner attaches.)
  env <- list2env(list(recombin = recombin, select = select, mutate = mutate,
                       Chrom = Chrom, FitnV = ranking(rowSums(Chrom))),
                  parent = emptyenv())
  expect_equal(dim(evalq(recombin("xovsp", Chrom), env)), dim(Chrom))
  expect_equal(dim(evalq(recombin("recdis", Chrom), env)), dim(Chrom))
  expect_equal(dim(evalq(select("sus", Chrom, FitnV), env)), dim(Chrom))
  expect_equal(dim(evalq(mutate("mut", Chrom), env)), dim(Chrom))
})

test_that("user-defined operators are found from the calling environment", {
  set.seed(1)
  Chrom <- crtbp(6, 4)$Chrom
  FitnV <- ranking(rowSums(Chrom))
  swap <- function(OldChrom) OldChrom[nrow(OldChrom):1, , drop = FALSE]
  first <- function(FitnV, Nsel) rep(1L, Nsel)
  flip <- function(OldChrom, Pm, BaseV) 1 - OldChrom
  expect_equal(recombin("swap", Chrom), Chrom[6:1, ])
  expect_equal(select("first", Chrom, FitnV), Chrom[rep(1, 6), ])
  expect_equal(select(first, Chrom, FitnV), Chrom[rep(1, 6), ])
  expect_equal(mutate("flip", Chrom), 1 - Chrom)
  expect_error(recombin("no_such_operator", Chrom), "no_such_operator")
  expect_error(select("no_such_operator", Chrom, FitnV), "no_such_operator")
  expect_error(mutate(1, Chrom), "MUT_F must be a string")
})
