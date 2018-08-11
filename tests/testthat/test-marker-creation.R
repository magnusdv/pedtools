context("marker creation")

test_that("paramlink syntax raises error", {
  x = nuclearPed(1)
  ERR = "Genotype assignments in `...` must be named."
  expect_error(marker(x, 0), ERR)
  expect_error(marker(x, 1, 'A'), ERR)
  expect_error(marker(x, '1', 1, '2', 1:2, alleles=1:2), ERR)
})


test_that("empty markers are created correctly", {
  x = nuclearPed(2)

  # Adding an empty SNP (all genotypes are missing):
  expect_is(marker(x), "marker")
  expect_is(marker(x, '1'=0), "marker")
})

test_that("alleles and freqs are sorted together", {
  x = nuclearPed(1)
  p = .8; q = 1-p

  # empty
  m1 = marker(x, alleles=2:1, afreq=c(q, p))
  expect_equal(alleles(m1), c('1','2'))
  expect_equal(unname(afreq(m1)), c(p, q))

  # the "2" allele observed
  m2 = marker(x, '1'=2, alleles=1:2, afreq=c(p, q))
  expect_equal(alleles(m2), c('1','2'))
  expect_equal(unname(afreq(m2)), c(p, q))

})
