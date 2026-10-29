test_that("crtbase builds base vectors", {
  expect_equal(crtbase(4), rep(2, 4))
  expect_equal(crtbase(c(3, 2), c(8, 4)), c(8, 8, 8, 4, 4))
  expect_equal(crtbase(c(2, 2), 5), rep(5, 4))
  expect_error(crtbase(2, c(2, 3)))
})

test_that("crtbp supports all argument forms", {
  set.seed(1)
  res <- crtbp(10, 6)
  expect_equal(dim(res$Chrom), c(10, 6))
  expect_true(all(res$Chrom %in% 0:1))

  res <- crtbp(c(10, 6))
  expect_equal(dim(res$Chrom), c(10, 6))

  res <- crtbp(10, 5, 8)
  expect_equal(res$BaseV, rep(8, 5))
  expect_true(all(res$Chrom >= 0 & res$Chrom < 8))

  res <- crtbp(10, c(2, 4, 8))
  expect_equal(res$Lind, 3)
  expect_equal(res$BaseV, c(2, 4, 8))
  expect_true(all(sweep(res$Chrom, 2, c(2, 4, 8), "<")))

  res <- crtbp(c(10, 3), c(2, 4, 8))
  expect_equal(res$BaseV, c(2, 4, 8))

  expect_error(crtbp(10))
  expect_error(crtbp(10, 3, c(2, 4)))
})

test_that("crtrp respects the bounds", {
  set.seed(1)
  FieldDR <- matrix(c(-1, 10, 1, 20), 2, 2, byrow = TRUE)
  Chrom <- crtrp(50, FieldDR)
  expect_equal(dim(Chrom), c(50, 2))
  expect_true(all(Chrom[, 1] >= -1 & Chrom[, 1] <= 1))
  expect_true(all(Chrom[, 2] >= 10 & Chrom[, 2] <= 20))
})

test_that("bs2rv decodes binary and Gray strings", {
  FieldD <- list(prec = 3, lb = 0, ub = 7, code = "binary", scale = "arith",
                 lbin = TRUE, ubin = TRUE)
  Chrom <- rbind(c(0, 0, 0), c(1, 1, 1), c(1, 0, 1))
  expect_equal(as.vector(bs2rv(Chrom, FieldD)), c(0, 7, 5))

  ## Gray code 111 is binary 101
  FieldD$code <- "gray"
  expect_equal(as.vector(bs2rv(Chrom, FieldD)), c(0, 5, 6))

  ## Excluding the bounds moves the values strictly inside the range
  FieldD <- list(prec = 3, lb = 0, ub = 1, code = "binary", scale = "arith",
                 lbin = FALSE, ubin = FALSE)
  Phen <- bs2rv(Chrom, FieldD)
  expect_true(all(Phen > 0 & Phen < 1))
})

test_that("bs2rv handles per-variable settings", {
  Chrom <- rbind(c(0, 0, 0, 0), c(1, 1, 1, 1))
  FieldD <- list(prec = 2, lb = c(0, 1), ub = c(3, 1000), code = "binary",
                 scale = c("arith", "log"), lbin = TRUE, ubin = TRUE)
  Phen <- bs2rv(Chrom, FieldD)
  expect_equal(Phen[, 1], c(0, 3))
  expect_equal(Phen[, 2], c(1, 1000))
  expect_error(bs2rv(Chrom, list(prec = 3)))
  FieldD$lb <- c(0, -1)
  expect_error(bs2rv(Chrom, FieldD), "Log-scaled")
})
