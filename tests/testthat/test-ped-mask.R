test_that("masking/unmasking works for lists", {
  x = list(Ped1 = singleton("Mr. X"), Ped2 = nuclearPed()) |>
    addMarker(name = "bcd", afreq = c(b=0.2, c = 0.3, d = 0.5), `Mr. X` = "c/c", `3` = "d/d") |>
    addMarker(name = "num", alleles = c("6", "6.1", "6.2", "7"), `3` = "6.1/7") |>
    setMutmod(marker = "num", model = "stepwise", rate = 0.1, rate2 = 1e-6, range = 0.1) |>
    setMutmod(marker = "bcd", model = "equal", rate = 0.1)

  y = maskPed(x, seed = 123)
  # plot(y$maskedPed, marker = 1:2)

  z = unmaskPed(y$maskedPed, keys = y$keys)
  expect_identical(x,z)
})
