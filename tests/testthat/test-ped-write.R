
wrped = function(x, ...)
  writePed(x, prefix = tempfile("test"), verbose = F, ...)


test_that("writePed is reversed by readPed - singleton", {
  x = singleton("a")

  # Without famid
  f1 = wrped(x, famid = F, header = T)
  expect_identical(x, readPed(f1))

  # With famid
  famid(x) = "A"
  f2 = wrped(x, famid = T, header = T)
  expect_identical(x, readPed(f2))
})

test_that("writePed is reversed by readPed - ped with marker", {
  x = nuclearPed(fa = "fa", nch = 1) |>
    addMarker(fa = "a/b", name = "m1")

  # Without famid
  f1 = wrped(x, famid = F, header = T)
  expect_identical(x, readPed(f1, sep = "/"))

  # With famid
  famid(x) = "A"
  f2 = wrped(x, famid = T, header = T)
  expect_identical(x, readPed(f2, sep = "/"))
})

test_that("writePed is reversed by readPed - pedlist", {
  x = addMarker(singleton(1), "1" = 1:2, name = "m")
  xx = list(p1 = x, p2 = x)

  f1 = wrped(xx, famid = T, header = T)
  y = readPed(f1, sep = "/")
  y = lapply(y, `famid<-`, value = "") # remove FAMID slots
  expect_identical(xx, y)
})

test_that("writePed is reversed by readPed - ped+freq", {
  x = nuclearPed(fa = "fa", nch = 1) |>
  addMarker(fa = "a/b", afreq = c(a=.1, b=.9), name = "m1")

  # Without famid
  ff = wrped(x, famid = F, header = T, what = c("ped", "freq"))
  expect_identical(x, readPed(ff[1], sep = "/", locus = readFreqDatabase(ff[2])))
})

