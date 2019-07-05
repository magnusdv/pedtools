context("ped from writePed output")

test_that("readPed loads a simple pedigree from the output of writePed", {
  saved = nuclearPed(1, sex = 2)
  famid(saved) = 'TheIncreadibles'
  prefix = paste(tempdir(), 'myped')
  writePed(saved, prefix = prefix)

  loaded = readPed(prefix)
  expect_identical(saved, loaded)
})

test_that('readPed loads a pedigree with marker data', {
  saved = nuclearPed(1, father = 'Father', mother = 'Mother')
  famid(saved) = 1
  m = marker(saved, Father = c(1,3), alleles = 1:4, afreq = c(0.2, 0.3, 0.1, 0.4), name='TH01')
  saved = addMarkers(saved, m)

  m = marker(saved, Mother = c('X', 'X'), alleles = c('X', 'Y'), afreq = c(0.75, 0.25), name = 'AMEL')
  saved = addMarkers(saved, m)

  prefix = paste(tempdir(), '/withmarkers', sep = '')
  writePed(saved, prefix = prefix)

  loaded = readPed(prefix)
  print(loaded$markerdata[[1]])
  print(saved$markerdata[[1]])
  expect_identical(saved, loaded)
})
