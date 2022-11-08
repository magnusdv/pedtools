
test_that("singleton with marker plots without error", {
  s = singleton(1)
  pdf(NULL) # avoid Rplots.pdf being created
  expect_silent(plot(s, marker = marker(s)))
})

