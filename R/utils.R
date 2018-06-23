
# Test that input is a nonegative integer.
#' @importFrom assertthat assert_that
is_count0 <- function(x) {
  assert_that(is.numeric(x), length(x) == 1, x == as.integer(x), x >= 0)
}

# Stripped-down, faster version of expand.grid
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

# Utility function for generating numbered "NN" labels.
# Returns "NN_i" where i increments largest j occuring as NN_j, NN.j or NN-j in input.
nextNN = function(labs) { # labs a character vector
  NNs = grepl("^NN", labs)
  if(!any(NNs))
    return("NN_1")
  NNnum = suppressWarnings(as.numeric(sub("^NN[._-]?", "", labs[NNs])))
  if(all(is.na(NNnum)))
    return("NN_1")
  nextNN = max(NNnum, na.rm=T) + 1
  return(sprintf("NN_%d", nextNN))
}

.mysetdiff = function(x, y) unique.default(x[match(x, y, 0L) == 0L])
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

.prettycat = function(v, andor) switch(min(len <- length(v), 3), toString(v), paste(v, collapse = " and "),
    paste(paste(v[-len], collapse = ", "), andor, v[len]))

