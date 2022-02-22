
test_that("full subset is identical to starting point", {
  s = singleton(1)
  expect_identical(s, subset(s, 1))

  x = reorderPed(fullSibMating(1), 6:1)
  expect_identical(x, subset(x, labels(x)))
})

test_that("subset() catches errors", {
  x = nuclearPed(1)
  expect_error(subset(x, 4), "Unknown ID label: 4")
  expect_error(subset(x, c(1,1,1,2,2,2)), "Duplicated ID label: 1, 2")
})

test_that("branch() catches errors", {
  x = nuclearPed(1)
  expect_error(branch(x, NULL), "`id` cannot be empty")
  expect_error(branch(x, 1:2), "`id` must contain a single ID label")
  expect_error(branch(x, 4), "Unknown ID label: 4")
})

test_that("subset() preserves marker annotations", {
  x = nuclearPed(1) |>
    addMarker('3' = "2/2", alleles = 1:2, afreq = c(0.1,0.9), name = "m1")

  s = singleton(3) |>
    addMarker('3' = "2/2", alleles = 1:2, afreq = c(0.1,0.9), name = "m1")

  expect_identical(subset(x, 3), s)
})


