test_that("linear ranking assigns fitness in [2 - SP, SP]", {
  FitnV <- ranking(c(4, 1, 3, 2, 5))
  expect_equal(FitnV, c(0.5, 2, 1, 1.5, 0))
  FitnV <- ranking(c(4, 1, 3, 2, 5), RFun = 1.5)
  expect_equal(range(FitnV), c(0.5, 1.5))
  expect_equal(ranking(c(4, 1, 3, 2, 5), RFun = NULL), ranking(c(4, 1, 3, 2, 5)))
})

test_that("ranking averages fitness of ties", {
  FitnV <- ranking(c(1, 2, 2, 3))
  expect_equal(FitnV[2], FitnV[3])
  expect_equal(sum(FitnV), 4)
})

test_that("non-linear ranking and subpopulations", {
  FitnV <- ranking(1:10, RFun = c(3, 1))
  expect_equal(sum(FitnV), 10)
  expect_true(all(diff(FitnV) < 0))
  FitnV <- ranking(c(1:5, 1:5), SUBPOP = 2)
  expect_equal(FitnV[1:5], FitnV[6:10])
  expect_error(ranking(1:5, SUBPOP = 2))
})

test_that("scaling keeps the average and handles constant values", {
  ObjV <- c(1, 2, 3, 4, 10)
  FitnV <- scaling(ObjV)
  expect_equal(mean(FitnV), mean(ObjV))
  expect_true(all(FitnV >= 0))
  expect_equal(scaling(rep(3, 4)), rep(3, 4))
  expect_error(scaling(c(-1, 2)))
})

test_that("rws and sus return valid indices", {
  set.seed(1)
  FitnV <- c(1.5, 1.35, 1.21, 1.07, 0.92, 0.78, 0.64, 0.5)
  expect_true(all(rws(FitnV, 6) %in% 1:8))
  idx <- sus(FitnV, 8)
  expect_length(idx, 8)
  expect_true(all(idx %in% 1:8))
  ## zero fitness individuals are never selected
  expect_true(all(sus(c(0, 1, 0, 1), 10) %in% c(2, 4)))
})

test_that("select works with default and explicit selection functions", {
  set.seed(1)
  Chrom <- crtbp(20, 5)$Chrom
  FitnV <- ranking(rowSums(Chrom))
  expect_equal(dim(select(Chrom = Chrom, FitnV = FitnV)), c(20, 5))
  expect_equal(dim(select("rws", Chrom, FitnV, GGAP = 0.5)), c(10, 5))
  expect_equal(dim(select("sus", Chrom, FitnV, SUBPOP = 2)), c(20, 5))
})

test_that("reins replaces the worst parents with the best offspring", {
  Chrom <- matrix(1:5, 5, 1)
  ObjVCh <- c(1, 5, 2, 4, 3)
  SelCh <- matrix(c(100, 200), 2, 1)
  ObjVSel <- c(0.5, 0.1)
  ## all offspring are inserted in place of the two worst parents (2 and 4)
  res <- reins(Chrom, SelCh, InsOpt = c(1, 1), ObjVCh = ObjVCh, ObjVSel = ObjVSel)
  expect_equal(as.vector(res$Chrom), c(1, 100, 3, 200, 5))
  expect_equal(res$ObjVCh, c(1, 0.5, 2, 0.1, 3))

  res <- reins(Chrom, SelCh, InsOpt = c(1, 0.5), ObjVCh = ObjVCh, ObjVSel = ObjVSel)
  expect_equal(as.vector(res$Chrom), c(1, 200, 3, 4, 5))

  expect_silent(reins(Chrom, SelCh, InsOpt = 0))
  expect_error(reins(Chrom, SelCh, InsOpt = c(0, 0.5)), "ObjVSel")
  expect_error(reins(Chrom, SelCh, InsOpt = c(1, 1)), "ObjVCh")
})

test_that("migrate exchanges individuals for every structure", {
  set.seed(1)
  Chrom <- matrix(rep(1:4, each = 5), 20, 3)
  ObjV <- Chrom[, 1]
  for (s in c("net", "neighbourhood", "ring")) {
    res <- migrate(Chrom, SUBPOP = 4, ObjV, Structure = s)
    expect_equal(dim(res$Chrom), dim(Chrom), info = s)
    expect_equal(res$ObjV, res$Chrom[, 1], info = s)
    ## the first individual of every subpopulation came from elsewhere
    expect_true(all(res$Chrom[c(1, 6, 11, 16), 1] != 1:4), info = s)
  }
  res <- migrate(Chrom, SUBPOP = 4, ObjV, Select = TRUE, Structure = "ring")
  expect_equal(res$Chrom[c(1, 6, 11, 16), 1], c(4, 1, 2, 3))
  expect_null(migrate(Chrom, SUBPOP = 4)$ObjV)
  expect_error(migrate(Chrom, SUBPOP = 4, Structure = "star"))
})
