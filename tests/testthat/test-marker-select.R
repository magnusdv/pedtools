context("marker selection")

x = setMarkers(singleton(1), locusAttributes = list(list(name = "m1"),
                                                    list(name = "m2")))
mnames = function(x) name(x, seq_markers(x))

test_that("selectMarkers() works with marker names", {
  expect_equal(mnames(selectMarkers(x, "m1")), "m1")
  expect_equal(mnames(selectMarkers(x, c("m2", "m1"))), c("m2", "m1"))
  expect_equal(mnames(selectMarkers(x, character(0))), character(0))

  expect_equal(mnames(removeMarkers(x, "m1")), "m2")
  expect_equal(mnames(removeMarkers(x, c("m2", "m1"))), character(0))
  expect_equal(mnames(removeMarkers(x, character(0))), c("m1", "m2"))
})

test_that("selectMarkers() works with numeric indices", {
  expect_equal(mnames(selectMarkers(x, 1)), "m1")
  expect_equal(mnames(selectMarkers(x, 2:1)), c("m2", "m1"))
  expect_equal(mnames(selectMarkers(x, numeric(0))), character(0))
})

test_that("selectMarkers() works with negative indices", {
  expect_equal(mnames(selectMarkers(x, -1)), "m2")
  expect_equal(mnames(selectMarkers(x, -2)), "m1")
  expect_equal(nMarkers(selectMarkers(x, -(1:2))), 0)
})
