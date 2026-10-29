## Population with alternating all-zero / all-one individuals, so that any
## crossover is visible as a mix of zeros and ones
pairs_population <- function(Npairs, Lind) {
  Chrom <- matrix(0, 2 * Npairs, Lind)
  Chrom[seq(2, 2 * Npairs, 2), ] <- 1
  Chrom
}

test_that("crossover respects Px", {
  set.seed(1)
  Chrom <- pairs_population(200, 10)
  for (f in c("xovsp", "xovdp", "xovsh", "xovsprs", "xovdprs", "xovshrs")) {
    expect_equal(match.fun(f)(Chrom, Px = 0), Chrom, info = f)
  }
  New <- xovsp(Chrom, Px = 1)
  crossed <- rowSums(New[seq(1, 400, 2), ]) %in% 1:9
  expect_true(all(crossed))
})

test_that("single-point crossover cuts at random positions", {
  set.seed(1)
  Chrom <- pairs_population(500, 10)
  New <- xovsp(Chrom, Px = 1)
  ## leading zeros kept by the first offspring = position of the cut
  cuts <- rowSums(New[seq(1, 1000, 2), ] == 0)
  expect_setequal(unique(cuts), 1:9)
})

test_that("crossover preserves genes and handles odd population sizes", {
  set.seed(1)
  Chrom <- crtbp(7, 12)$Chrom
  New <- xovmp(Chrom, Px = 1, Npt = 2)
  expect_equal(dim(New), dim(Chrom))
  expect_equal(colSums(New[1:6, ]), colSums(Chrom[1:6, ]))
  expect_equal(New[7, ], Chrom[7, ])
})

test_that("real-valued recombination keeps the shape", {
  set.seed(1)
  FieldDR <- matrix(c(-100, -50, 100, 50), 2, 2, byrow = TRUE)
  Chrom <- crtrp(5, FieldDR)
  for (f in c("recdis", "recint", "reclin")) {
    New <- match.fun(f)(Chrom)
    expect_equal(dim(New), dim(Chrom), info = f)
    expect_equal(New[5, ], Chrom[5, ], info = f)
  }
  New <- recmut(Chrom, FieldDR)
  expect_false(anyNA(New))
  expect_true(all(New[, 1] >= -100 & New[, 1] <= 100))
})

test_that("recmut moves offspring in opposite directions", {
  set.seed(1)
  FieldDR <- matrix(c(-100, 100), 2, 1)
  Chrom <- matrix(c(0, 10), 2, 1)
  for (i in 1:20) {
    ## shrink the range so that no offspring is clipped at the bounds
    New <- recmut(Chrom, FieldDR, MutOpt = c(1, 0.01))
    expect_equal(New[1, 1] - Chrom[1, 1], -(New[2, 1] - Chrom[2, 1]))
  }
  ## parents with equal norm must not produce NaN
  expect_false(anyNA(recmut(matrix(c(5, -5), 2, 1), FieldDR)))
})

test_that("recombin dispatches to low-level functions", {
  set.seed(1)
  Chrom <- crtbp(20, 8)$Chrom
  expect_equal(dim(recombin("xovsp", Chrom)), dim(Chrom))
  expect_equal(dim(recombin("xovmp", Chrom, Npt = 2, SUBPOP = 2)), dim(Chrom))
  expect_equal(recombin("xovsp", Chrom, RecOpt = 0), Chrom)
  expect_error(recombin("recmut", Chrom))
})

test_that("mutation works for discrete and real values", {
  set.seed(1)
  Chrom <- crtbp(20, 10)$Chrom
  New <- mut(Chrom, Pm = 1)
  expect_equal(New, 1 - Chrom)
  expect_equal(mutate(OldChrom = Chrom, MutOpt = 0), Chrom)

  FieldDR <- matrix(c(-1, -1, 1, 1), 2, 2, byrow = TRUE)
  Chrom <- crtrp(10, FieldDR)
  New <- mutate("mutbga", Chrom, FieldDR, MutOpt = c(1, 1))
  expect_equal(dim(New), dim(Chrom))
  expect_true(all(abs(New) <= 1))
})
