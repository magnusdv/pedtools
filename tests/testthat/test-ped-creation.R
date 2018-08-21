context("ped construction")

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

  expect_error(ped(id, fid, mid, sex, famid=1:2, val = F), "`famid` must have length 1")
  expect_error(ped(id, fid, mid, sex, famid=NULL, val = F), "`famid` must have length 1")

  expect_error(ped(c(1,2,1), fid, mid, sex, val = F), "Duplicated entry in `id` vector: 1")

  expect_error(ped(letters[1:3], fid, mid, sex, val = F), "`fid` entry does not appear in `id` vector: 1")
  expect_error(ped(id, fid, letters[1:3], sex, val = F), "`mid` entry does not appear in `id` vector: a, b, c")
})


test_that("validate_ped() catches malformed pedigrees", {
  id=letters[1:3]; fid=c("b",0,0); mid=c("c",0,0); sex=c(1,1,2)

  expect_error(ped(id, c(0,0,0), mid, sex), "Individual a has exactly 1 parent")
  expect_error(ped(id, fid, c(0,0,0), sex), "Individual a has exactly 1 parent")

  expect_error(ped(id, fid, mid, c(1,3,-1)), "Illegal gender code: 3")
  expect_error(ped(id, fid, mid, c(1,3,-1)), "Illegal gender code: -1")

  expect_error(ped(id, mid, fid, sex), "Individual c is female, but appear as the father of a")
  expect_error(ped(id, mid, fid, sex), "Individual b is male, but appear as the mother of a")

})

test_that("self-ancestry is detected", {
  expect_error(ped(id=1, fid=1, mid=1, sex=0),
               "Individual 1 is their own ancestor")
  expect_error(ped(id=1:3, fid=c(0,1,1), mid=c(0,3,2), sex=c(1,2,2)),
               "Individual 3 is their own ancestor")
})

test_that("nuclearPed() works by giving labels", {
  expect_identical(nuclearPed(1), nuclearPed(children='3'))
  expect_identical(nuclearPed(1, sex=2), nuclearPed(children='3', sex=2))
  expect_identical(relabel(nuclearPed(1), letters[3:1]),
                   nuclearPed(fa='c', mo='b', children='a'))
  expect_identical(relabel(nuclearPed(1), old=3, new="foo"),
                   nuclearPed(children="foo"))
})

test_that("nuclearPed() gives sensible error messages", {
  expect_error(nuclearPed(), 'argument "nch" is missing')
  expect_error(nuclearPed(0), '`nch` must be a positive integer: 0')
  expect_error(nuclearPed(nch = 2, child = 3), "`children` must have length `nch`")

  expect_error(nuclearPed(fa = 1:2, nch = 1), "`father` must have length 1")
  expect_error(nuclearPed(mo = 1:2, nch = 1), "`mother` must have length 1")

  expect_error(nuclearPed(nch = 1, sex='a'), "Illegal gender code: a")
  expect_error(nuclearPed(nch = 1, sex=1:2), "`sex` is longer than the number of individuals")
  expect_error(nuclearPed(nch = 1, sex=integer(0)), "`sex` cannot be empty")

  expect_error(nuclearPed(fa = 'a', child = 'a'), "Duplicated ID label: a")
  expect_error(nuclearPed(child = c("b", "b")), "Duplicated ID label: b")

  expect_error(nuclearPed(child = 1), "please specify a different label for the father")
  expect_error(nuclearPed(mo = 1, nch = 1), "please specify a different label for the father")
  expect_error(nuclearPed(child = 2), "please specify a different label for the mother")
  expect_error(nuclearPed(fa = 2, nch = 1), "please specify a different label for the mother")

})

test_that("simple ped", {
  x1 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), reorder=T)
  x2 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), famid="")
  x3 = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1), reorder=F)
  expect_is(x1, "ped")
  expect_identical(x1,x2)
  expect_identical(x1,x3)
})


test_that("random ped", {
  x = randomPed(3, 3)
  expect_is(x, "ped")
  expect_equal(pedsize(x), 6)
})

test_that("singleton creation works as expected", {
  x = singleton(1)
  expect_is(x, "ped")
  expect_is(x, "singleton")

  expect_error(singleton(), 'argument "id" is missing')
})

