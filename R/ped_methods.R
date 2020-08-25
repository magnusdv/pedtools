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
  showmess = F
  if (missing(markers)) {
    if (nm < 10)
      markers = seq_len(nm)
    else {
      markers = 1:5
      showmess = T
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
  cat(sprintf("Pedigree with %d members\n", pedsize(x)))
  cat(nMarkers(x), "attached markers\n")
  cat(length(typedMembers(x)), "typed members\n")
}

#' @export
summary.singleton = function(object, ...) {
  cat("Singleton\n")
  cat(nMarkers(object), "attached markers\n")
}

#' @export
summary.list = function(object, ...) {
  x = object
  if(!is.pedList(x))
    return(summary.default(x))
  nInd = pedsize(x)
  nMark = sapply(x, nMarkers)

  # Special treatment if all singletons
  if (all(nInd == 1) && all(nMark == nMark[1])) {
    cat(sprintf("List of %d singletons\n%d attached markers\n", length(x), nMark[[1]]))
    return(invisible(NULL))
  }

  cat("List of", length(x), "`ped` objects:\n")
  for(i in seq_along(x)) {
    cat(sprintf("\n--- component %d ---\n", i))
    summary(x[[i]])
  }
}
