context("ped various")

test_that("example ped from identity package", {
  ex = read.table (system.file ("example", "ex.pedigree", package = "identity"))
  x = ped(id=ex[,1], fid=ex[,2], mid=ex[,3], sex=rep(0,nrow(ex)), reorder=F)
  expect_identical(restore_ped(as.matrix(x)), x)

  y = ped(id=ex[,1], fid=ex[,2], mid=ex[,3], sex=rep(0,nrow(ex)), famid="identity-example")
  expect_identical(restore_ped(as.matrix(y)), y)

  z = reorderPed(y, order(as.numeric(labels(y))))
  expect_identical(restore_ped(as.matrix(z)), z)
})

test_that("ped -> matrix -> ped", {
  x = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1))
  y = restore_ped(as.matrix(x))
  expect_identical(x, y)

  x2 = reorderPed(x, 3:1)
  y2 = restore_ped(as.matrix(x2))
  expect_identical(x2, y2)
})

