context("marker transfers")


test_that("single ped is unchaged after self-transfer of markers", {
  x = nuclearPed(1)
  m = marker(x, '1'=1:2)
  x = setMarkers(x, m)

  expect_identical(transferMarkers(x,x), x)
})
test_that("singleton is unchaged after self-transfer of markers", {
  x = singleton(1)
  m = marker(x, '1'=1:2)
  x = setMarkers(x, m)

  expect_identical(transferMarkers(x,x), x)
})

test_that("ped is unchanged after back-and-fourth transfer to singleton", {
  x = nuclearPed(1)
  m = marker(x, '1'=1:2)
  x = setMarkers(x, m)

  y = transferMarkers(x, singleton(id='1'))
  z = transferMarkers(y, x)
  expect_identical(x,z)

  w = transferMarkers(z, y)
  expect_identical(y,w)
})

test_that("pedlist is unchaged after self-transfer of markers", {

  x = nuclearPed(children="boy")
  x = setMarkers(x, marker(x, boy="a", alleles=c("a", "b")))

  y = singleton("otherboy")
  y = setMarkers(y, marker(y, otherboy="b", alleles=c("a", "b")))
  pedlist = list(x, y)

  expect_identical(transferMarkers(pedlist,pedlist), pedlist)
})
