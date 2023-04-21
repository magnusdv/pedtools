#' Printing pedigrees
#'
#' Print a `ped` object using original labels.
#'
#' This first calls [as.data.frame.ped()] and then prints the resulting
#' data.frame. The data.frame is returned invisibly.
#'
#' @param x object of class `ped`.
#' @param ... (optional) arguments passed on to [print.data.frame()].
#' @param markers (optional) vector of marker indices. If missing, and `x` has
#'   less than 10 markers, they are all displayed. If `x` has 10 or more
#'   markers, the first 5 are displayed.
#' @param verbose If TRUE, a message is printed if only the first 5 markers are
#'   printed. (See above).
#' @export
print.ped = function(x, ..., markers, verbose = TRUE) {
  nm = nMarkers(x)
  showmess = FALSE
  if (missing(markers)) {
    if (nm < 10)
      markers = seq_len(nm)
    else {
      markers = 1:5
      showmess = TRUE
    }
  }
  else {
    if (any(markers > nm)) stop2("Markers out of range: ", markers[markers > nm])
  }
  datafr = as.data.frame(x, markers = markers)
  datafr$fid[datafr$fid == "0"] = "*"
  datafr$mid[datafr$mid == "0"] = "*"

  # Add question marks at any genotypes shoing male X heterozygosity
  datafr = questionMaleHetX(x, datafr)

  print(datafr, row.names = FALSE, ...)

  if(showmess && verbose)
    message("Only 5 (out of ", nm, ") markers are shown.")

  invisible(datafr)
}


# Function for adding question mark at genotypes showing male X heterozygosity
questionMaleHetX = function(x, df) {
  if(!hasMarkers(x) || ncol(df) < 5)
    return(df)

  dfnames = names(df)
  xnames = name(x, 1:nMarkers(x))

  # Loop through all marker columns of the data frame
  for(i in 5:ncol(df)) {
    mname = dfnames[i]
    if(startsWith(mname, "<"))
      midx = as.integer(substring(mname, 2, nchar(mname)-1))
    else
      midx = match(mname, xnames)

    if(!is.na(midx) && isXmarker(x, midx)) {
      maleHet = df[["sex"]] == 1 & grepl("/", df[[i]])
      df = commentAndRealign(df, i, maleHet, "?")
    }
  }

  # Returned modified data.frame
  df
}


#### Summary methods ####

#' @export
summary.ped = function(object, ...) {
  x = object
  cat(sprintf("Pedigree with %d members (%d males, %d females, %d unknown).\n",
              pedsize(x), sum(x$SEX == 1), sum(x$SEX == 2), sum(x$SEX == 0)))
  cat(sprintf("%d generations, %d founders, %d leaves.\n",
              generations(x, "max"), length(founders(x)), length(leaves(x))))
  nm = nMarkers(x)
  na = nAlleles.ped(x)
  if(nm == 0)
    cat("0 attached markers.\n")
  else if(nm == 1)
    cat(sprintf("1 attached marker (%d alleles).\n", na))
  else if(min(na) == max(na))
    cat(sprintf("%d attached markers (all with %d alleles).\n", nm, min(na)))
  else
    cat(sprintf("%d attached markers (%d - %d alleles).\n", nm, min(na), max(na)))
  cat(length(typedMembers(x)), "typed members.\n")
}

#' @export
summary.singleton = function(object, ...) {
  x = object
  SX = c("sex unknown", "male", "female")
  cat(sprintf("Singleton (%s) labelled '%s'.\n", SX[x$SEX + 1], x$ID))

  nm = nMarkers(x)
  na = nAlleles.ped(x)
  if(nm == 0)
    cat("0 attached markers.\n")
  else if(nm == 1)
    cat(sprintf("1 attached marker (%d alleles).\n", na))
  else if(min(na) == max(na))
    cat(sprintf("%d attached markers (all with %d alleles).\n", nm, min(na)))
  else
    cat(sprintf("%d attached markers (%d - %d alleles).\n", nm, min(na), max(na)))
}

#' @export
summary.list = function(object, ...) {
  x = object

  if(!is.pedList(x))
    return(summary.default(x))

  nInd = pedsize(x)
  nTot = sum(nInd)
  nMark = sapply(x, nMarkers)

  # Special treatment if all singletons
  if (all(nInd == 1) && all(nMark == nMark[1])) {
    cat(sprintf("List of %d singletons.\n", nTot))
    if(nTot < 10) {
      ids = unlist(labels(x), use.names = FALSE)
      sex = c("sex unknown", "male", "female")[getSex(x) + 1]
      lbs = toString(sprintf("%s (%s)", ids, sex))
      cat(sprintf("Labels: %s.\n", lbs))
    }
    cat(nMark[1], "attached markers.\n")
    return(invisible(NULL))
  }

  cat("List of", length(x), "connected pedigrees")
  if(all(nInd == nInd[1]))
    cat(sprintf(" (each with %d members).\n"), nInd[1])
  else
    cat(sprintf(" (%d - %d members).\n", min(nInd), max(nInd)))

  sex = getSex(x)
  cat(sprintf("In total %d individuals (%d males, %d females, %d unknown).\n",
              nTot, sum(sex == 1), sum(sex == 2), sum(sex == 0)))
  if(length(x) < 5) {
    for(i in seq_along(x)) {
      cat(sprintf("\n--- component %d ---\n", i))
      summary(x[[i]])
    }
  }
}
