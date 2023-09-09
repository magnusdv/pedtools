
x = nuclearPed(1)

test_that("nAlleles of empty marker is 2", {
  m = marker(x)
  expect_equal(nAlleles(m), 2)
})

test_that("isXmarker() works", {
  expect_true(isXmarker(marker(x, chrom=23L)))
  expect_true(isXmarker(marker(x, chrom=23)))
  expect_false(isXmarker(marker(x)))
  expect_false(isXmarker(marker(x, chrom=1)))
})

test_that("nMarkers() and hasMarkers() works with multiple comps", {
  x = list(singleton(1), singleton(2) |> addMarker())
  expect_error(nMarkers(x), "Pedigree components have different number of markers")
  expect_identical(nMarkers(x, compwise = T), 0:1)

  expect_true(hasMarkers(x), "Pedigree components have different number of markers")
  expect_identical(hasMarkers(x, compwise = T), c(F,T))
})

test_that("setMarkers() attaches a list of markers", {
  m = marker(x)
  x1 = setMarkers(x, m)
  x2 = setMarkers(x, list(m))
  expect_identical(x1, x2)
  expect_equal(nMarkers(x1), 1)
  expect_is(x1$MARKERS, "markerList")
  expect_true(hasMarkers(x1))
})

test_that("getMap() produces a data.frame", {
  m1 = marker(x)
  m2 = marker(x, chrom=1, posMb=1)
  x1 = setMarkers(x, list(m1,m2))
  expect_is(getMap(x1, na=0, verbose=F), "data.frame")
  expect_is(getMap(x1, na=1, verbose=F), "data.frame")
  expect_is(getMap(x1, na=2, verbose=F), "data.frame")
  expect_is(getMap(x1, markers=1, na=0, verbose=F), "data.frame")
  expect_is(getMap(x1, markers=1, na=1, verbose=F), "data.frame")
  expect_is(getMap(x1, markers=1, na=2, verbose=F), "data.frame")
})
