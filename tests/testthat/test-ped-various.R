
test_that("ped -> matrix -> ped", {
  x = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1))
  y = restorePed(as.matrix(x))
  expect_identical(x, y)

  x2 = reorderPed(x, 3:1)
  y2 = restorePed(as.matrix(x2))
  expect_identical(x2, y2)
})

test_that("as.ped.data.frame() keeps genotypes", {
  # Without FAMID
  df = data.frame(id=1:2, fid=0, mid=0, sex=1, A1=1, A2=2)
  p = as.ped(df)
  expect_identical(genotype(p[[1]],id = 1, markers = 1), c("1", "2"))
  expect_identical(genotype(p[[2]],id = 2, markers = 1), c("1", "2"))

  # With FAMID, same ID
  df2 = data.frame(famid=1:2, id=1, fid=0, mid=0, sex=1, A1=1, A2=2)
  pp = as.ped(df2)
  expect_identical(genotype(pp[[1]],id = 1, markers = 1), c("1", "2"))
  expect_identical(genotype(pp[[2]],id = 1, markers = 1), c("1", "2"))
})
