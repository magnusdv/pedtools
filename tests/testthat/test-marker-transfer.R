
test_that("ped is unchanged after self-transfer of markers", {
  x = nuclearPed(1) |> addMarker("1" = "1/2")
  expect_identical(transferMarkers(x,x), x)

  # Singleton
  s = singleton(1) |> addMarker("1" = "1/2")
  expect_identical(transferMarkers(s,s), s)
})

test_that("ped is unchanged after back-and-fourth transfer to singleton", {
  x = nuclearPed(1) |> addMarker("1" = "1/2")

  y = x |> transferMarkers(to = singleton(1)) |> transferMarkers(to = x)
  expect_identical(x,y)
})

test_that("pedlist is unchanged after self-transfer of markers", {

  x = nuclearPed(children="boy") |>
    addMarker(boy="a", alleles=c("a", "b"))

  y = singleton("otherboy") |>
    addMarker(otherboy="b", alleles=c("a", "b"))

  pedlist = list(x, y)

  expect_identical(transferMarkers(pedlist,pedlist), pedlist)
})

test_that("transfer of multiple markers works when erase = F", {
  x = nuclearPed(father = "FA", children = "CH")
  m1 = marker(x, FA = 2, alleles = 1:2, name = 'locus1')
  m2 = marker(x, FA = 'a', alleles = c('a', 'b'), name = 'locus2')
  x = setMarkers(x, list(m1, m2))

  # Extract child as singleton (including locus annotations)
  ch = subset(x, "CH")

  # Add some genotypes and transfer back to trio
  genotype(ch, "locus1", "CH") = 1
  genotype(ch, "locus2", "CH") = 'b'
  y = transferMarkers(ch, x, erase = FALSE)

  # Compare with direct alternative
  xx = x
  genotype(xx, "locus1", "CH") = 1
  genotype(xx, "locus2", "CH") = 'b'
  expect_identical(y, xx)

  # Another test: Marker order in target should not matter
  y = transferMarkers(ch, selectMarkers(x, 2:1), erase = FALSE)
  expect_identical(y, xx)
})

test_that("transferMarkers checks for duplicated IDs", {
  x = singleton("a") |> addMarker(a = "1/2", name = "M")

  expect_error(transferMarkers(list(x,x), x), "Non-unique ID label in source ped: a")
  expect_error(transferMarkers(list(x,x), x, ids = "a"), "Non-unique ID label in source ped: a")

  expect_error(transferMarkers(x, list(x,x)), "Non-unique ID label in target ped: a")
  expect_error(transferMarkers(x, list(x,x), ids = "a"), "Non-unique ID label in target ped: a")

})

test_that("transferMarkers works with changing labels", {
  x = singleton("a") |> addMarker(a = "1/2", name = "M")
  y = singleton("b")
  expect_equal(transferMarkers(x, y, idsFrom = "a", idsTo = "b"),
               addMarker(y, b = "1/2", name = "M"))

  # Roundtrip
  x1 = nuclearPed() |> addMarker("1" = "1/2")

  x2 = transferMarkers(x1, x1, idsFrom=1, idsTo=2)
  x3 = transferMarkers(x2, x2, idsFrom=2, idsTo=3)
  y1 = transferMarkers(x3, x3, idsFrom=3, idsTo=1)
  expect_equal(x1, y1)
})
