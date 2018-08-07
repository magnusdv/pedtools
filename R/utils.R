stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call.=FALSE))
  do.call(stop, a)
}

stopifnotSimpleVector = function(x, argname="x") {
  if(is.null(x))
    return()

  if(!is.vector(x)) {
    errmess = sprintf("argument `%s` must be a vector", argname)

    cl = class(x)[1]
    if(!cl %in% c("numeric", "integer", "character", "logical", "double"))
      errmess = sprintf("%s; received an object of class '%s'", errmess, cl)

    stop2(errmess)
  }
}

# Test that input is a positive (or similar) integer.
is_count = function(x, minimum = 1) {
  isTRUE(length(x) == 1 &&
         is.integer(x) || (is.numeric(x) && x == as.integer(x)) &&
         x >= minimum)
}

# A safer version of base::sample
safe_sample <- function(x, ...) x[sample.int(length(x), ...)]


# Fast setdiff
.mysetdiff = function(x, y) unique.default(x[match(x, y, 0L) == 0L])

# Fast intersection. NB: assumes no duplicates!
.myintersect = function(x, y) y[match(x, y, 0L)]


.comb2 = function(n) {
    if (n < 2)
        return(matrix(nrow = 0, ncol = 2))
    v1 = rep.int(seq_len(n - 1), (n - 1):1)
    v2 = NULL
    for (i in 2:n) v2 = c(v2, i:n)
    cbind(v1, v2, deparse.level = 0)
}

.rand01 = function(n) sample.int(2, size = n, replace = T) - 1  #random 0/1 vector of length n.

# Stripped version of expand.grid
fast.grid = function(argslist, as.list = FALSE) {
  nargs = length(argslist)
  orep = nr = prod(lengths(argslist))
  if (nargs == 0L || nr == 0L)
    return(matrix(ncol = 0, nrow = 0))

  rep.fac = 1L
  res = NULL
  for (x in argslist) {
    nx = length(x)
    orep = orep/nx
    res = c(res, x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)])  #this is res[, i]
    rep.fac = rep.fac * nx
  }
  dim(res) = c(nr, nargs)
  if (as.list)
    res = lapply(seq_len(nr), function(r) res[r, ])
  res
}
