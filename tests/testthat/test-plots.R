context("pedigree plots")

test_that("singleton with marker plots without error", {
  pdf(NULL) # avoid Rplots.pdf being created
  expect_silent(plot(s, marker(s)))
})

