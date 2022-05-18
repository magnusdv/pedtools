
dm = function(...) getMap(distributeMarkers(singleton(1), ...), verbose = FALSE)

mapdf = function(CHROM, MARKER, MB)
  data.frame(CHROM = CHROM, MARKER = MARKER, MB = MB)

test_that("distributeMarkers() trivial cases", {
  expect_identical(dm(n = 1), mapdf("1", "M1", 0))
  expect_identical(dm(n = 1, chromLen = c(a=1)), mapdf("a", "M1", 0))
  expect_identical(dm(n = 1, prefix = "snp:"), mapdf("1", "snp:1", 0))

  expect_identical(dm(dist = 1, chromLen = 0.99), mapdf("1", "M1", 0))
  expect_identical(dm(dist = 1, chromLen = 1), mapdf("1", c("M1", "M2"), c(0,1)))

  expect_equal(nrow(dm(n = 3, chromLen = 1)), 3)
  expect_equal(nrow(dm(dist = .5, chromLen = 1)), 3)


  expect_identical(dm(n = 3, chromLen = 1:2, prefix = ""),
                   mapdf(as.character(c(1,2,2)), as.character(1:3), c(0, 0.5, 2)))
  expect_identical(dm(n = 3, chromLen = c("b"=1, "a"=2)),
                   mapdf(c("b", "a", "a"), paste0("M", 1:3), c(0, 0.5, 2)))
})

test_that("distributeMarkers() cathces bad input", {
  expect_error(dm(n = 0), "`n` must be a positive integer")
  expect_error(dm(n = 1:3), "`n` must be a positive integer")
  expect_error(dm(n = 1.5), "`n` must be a positive integer")
  expect_error(dm(n = "1"), "`n` must be a positive integer")

  expect_error(dm(dist = 0), "`dist` must be a positive number")
  expect_error(dm(dist = 1:3), "`dist` must be a positive number")
  expect_error(dm(dist = "1"), "`dist` must be a positive number")

  expect_error(dm(n=1, dist=1), "Exactly one of `n` and `dist` must be given")
  expect_error(dm(n = 1, chromLen = c(a=1, 2)), "Irregular chromosome names")
  expect_error(dm(n = 1, chromLen = c(a=1, a=2)), "Duplicated chromosome name")
})

test_that("distributeMarkers() assigns correct alleles and freqs", {
  s = singleton(1)
  expect_equal(nAlleles(distributeMarkers(s, n = 2, alleles = 1:3)), c(3,3))
  expect_equal(nAlleles(distributeMarkers(s, n = 1, alleles = 1:4)), c(4))

  afr = c(a=.3,b=.5, c=.2)
  expect_equal(afreq(distributeMarkers(s, n = 2, afreq = afr), 2), afr)

  afr2 = c(b=.3,a=.5, c=.2)
  expect_equal(afreq(distributeMarkers(s, n = 2, afreq = afr2), 2), afr2[sort(names(afr2))])
})

