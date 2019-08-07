context("allele matrix")


test_that("`getAlleles()` returns a character matrix", {
  df = data.frame(id = c("V1", "V2"), fid = 0, mid = 0, sex = 1, M = "1/2")
  x = as.ped(df, sep = "/")
  a = getAlleles(x, "V1")
  expect_true(is.character(a) && is.matrix(a) && all(dim(a) == 1:2))
})
