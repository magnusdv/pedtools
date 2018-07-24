context("marker accessors")


test_that("the simple marker getters work", {
  x = nuclearPed(1)
  m = marker(x, name="m1", chrom=1, posMb=1e7)
  x = setMarkers(x, m)

  expect_equal(name(m), "m1")
  expect_equal(name(x, 1), "m1")

  expect_equal(chrom(m), "1")
  expect_equal(chrom(x, markeridx=1), "1")
  expect_equal(chrom(x, markername="m1"), "1")

  expect_equal(posMb(m), 1e7)
  expect_equal(posMb(x, markeridx=1), 1e7)
  expect_equal(posMb(x, markername="m1"), 1e7)

  expect_equal(posCm(m), NA_real_)
  expect_equal(posCm(x, markeridx=1), NA_real_)
  expect_equal(posCm(x, markername="m1"), NA_real_)
})


test_that("alleles() works", {
  x = nuclearPed(1)
  als = c("p","e","d")
  m1 = marker(x, alleles=1:3, name="m1")
  m2 = marker(x, alleles=als, name="m2")
  x = setMarkers(x, list(m1,m2))

  expect_equal(alleles(m1), as.character(1:3))
  expect_equal(alleles(x, markeridx=1), as.character(1:3))
  expect_equal(alleles(x, markername="m1"), as.character(1:3))

  expect_equal(alleles(m2), sort(als))
  expect_equal(alleles(x, markeridx=2), sort(als))
  expect_equal(alleles(x, markername="m2"), sort(als))
})
