context("marker utils")

x = nuclearPed(1)

test_that("nAlleles of empty marker is 1", {
  m = marker(x)
  expect_equal(nAlleles(m), 1)
})

test_that("is_Xmarker() works", {
  expect_true(is_Xmarker(marker(x, chrom=23L)))
  expect_true(is_Xmarker(marker(x, chrom=23)))
  expect_false(is_Xmarker(marker(x)))
  expect_false(is_Xmarker(marker(x, chrom=1)))
})

test_that("setMarkers() attaches a list of markers", {
  m = marker(x)
  x1 = setMarkers(x, m)
  x2 = setMarkers(x, list(m))
  expect_identical(x1, x2)
  expect_equal(nMarkers(x1), 1)
  expect_is(x1$markerdata, "markerList")
  expect_true(hasMarkers(x1))
})

test_that("getMap() produces a data.frame", {
  m1 = marker(x)
  m2 = marker(x, chrom=1, posCm=1)
  x1 = setMarkers(x, list(m1,m2))
  expect_is(getMap(x1, na=0, verbose=F), "data.frame")
  expect_is(getMap(x1, na=1, verbose=F), "data.frame")
  expect_is(getMap(x1, na=2, verbose=F), "data.frame")
  expect_is(getMap(x1, markers=1, na=0, verbose=F), "data.frame")
  expect_is(getMap(x1, markers=1, na=1, verbose=F), "data.frame")
  expect_is(getMap(x1, markers=1, na=2, verbose=F), "data.frame")
})
