
test_that("ped() accepts any of ('', 0, NA) as missing parents", {
  s = singleton(1)
  expect_identical(ped(1, NA, NA, 1), s)
  expect_identical(ped(1, 0, 0, 1), s)
  expect_identical(ped(1, "0", "0", 1), s)
  expect_identical(ped(1, "", "", 1), s)
})

test_that("ped() catches various input errors", {
  id=1:3; fid=c(0,0,1); mid=c(0,0,2); sex=c(1,2,1)

  expect_error(ped(NULL, fid, mid, sex, val = F), "`id` vector has length 0")

  expect_error(ped(id, 0, mid, sex, val = F), "Incompatible input")
  expect_error(ped(id, fid, 0, sex, val = F), "Incompatible input")
  expect_error(ped(id, fid, mid, 0, val = F), "Incompatible input")

  expect_error(ped(id, fid, mid, sex, famid=1:2, val = F), "`famid` must be a character string")
  expect_error(ped(id, fid, mid, sex, famid=NULL, val = F), "`famid` must be a character string")

  expect_error(ped(c(1,2,1), fid, mid, sex, val = F), "Duplicated entry in `id` vector: 1")

  expect_error(ped(letters[1:3], fid, mid, sex, val = F), "`fid` entry does not appear in `id` vector: 1")
  expect_error(ped(id, fid, letters[1:3], sex, val = F), "`mid` entry does not appear in `id` vector: a, b, c")
})


test_that("validatePed() catches malformed pedigrees", {
  id=letters[1:3]; fid=c("b",0,0); mid=c("c",0,0); sex=c(1,1,2)

  expect_error(ped(id, c(0,0,0), mid, sex), "Individual a has exactly 1 parent")
  expect_error(ped(id, fid, c(0,0,0), sex), "Individual a has exactly 1 parent")

  expect_error(ped(id, fid, mid, c(1,3,-1)), "Illegal sex: 3, -1")
  expect_error(ped(id, fid, mid, c(1,2,"?")), "Illegal sex: ?")

  expect_error(ped(id, mid, fid, sex), "Individual c is female, but appear as the father of a")
  expect_error(ped(id, mid, fid, sex), "Individual b is male, but appear as the mother of a")

})

test_that("self-ancestry is detected", {
  expect_error(ped(id=1, fid=1, mid=1, sex=0),
               "Pedigree has no founders")
  expect_error(ped(id=1:3, fid=c(0,1,1), mid=c(0,3,2), sex=c(1,2,2)),
               "Individual 3 is their own ancestor")
})


test_that("simple ped", {
  x1 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), reorder=T)
  x2 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), famid="")
  x3 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), reorder=F)
  expect_is(x1, "ped")
  expect_identical(x1,x2)
  expect_identical(x1,x3)
})


test_that("randomPed works", {
  x = randomPed(6, 3, seed = 3)
  expect_is(x, "ped")
  expect_equal(pedsize(x), 6)

  expect_error(randomPed(2), "The total number of individuals must be at least 3")
  expect_error(randomPed(3,1), "When selfing is disallowed, the number of founders must be at least 2")
  expect_error(randomPed(3,3), "Too many founders")
})

test_that("singleton creation works as expected", {
  x = singleton(1)
  expect_is(x, "ped")
  expect_is(x, "singleton")
})

