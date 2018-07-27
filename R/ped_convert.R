#' Convert ped to matrix
#'
#' Converts a `ped` object to a numeric matrix using internal labels, with
#' additional info neccessary to recreate the original `ped` attached as
#' attributes.
#'
#' `restore_ped` is the reverse of `as.matrix.ped`.
#'
#' @param x a `ped` object. In `restore_ped`: A numerical matrix.
#' @param include.attrs a logical indicating if marker annotations and other
#'   info should be attached as attributes. See value.
#' @param attrs a list containing labels and other `ped` info compatible with
#'   `x`, in the format produced by `as.matrix`. If NULL, the attributes of `x`
#'   itself are used.
#' @param check a logical, forwarded to [ped()]. If FALSE, no checks for
#'   pedigree errors are performed.
#' @param \dots not used.
#'
#' @return For `as.matrix`: A numerical matrix with `pedsize(x)` rows.
#'   If `include.attrs = TRUE` the following attributes are added to the matrix,
#'   allowing `x` to be exactly reproduced by `restore_ped`:
#'
#' * `famid` the family identifier (a string)
#' * `labels` the ID labels (a character vector)
#' * `UNBROKEN_LOOPS` a logical indicating whether `x` has unboken loops
#' * `LOOP_BREAKERS` a numerical matrix, or NULL
#' * `markerattr` a list of length `nMarkers(x)`, containing the attributes of each marker
#'
#'  For `restore_ped`: A `ped` object.
#' @author Magnus Dehli Vigeland
#' @seealso [ped()]
#'
#' @examples
#'
#' x = setLabels(nuclearPed(1), letters[1:3])
#'
#' # To examplify the ped -> matrix -> ped trick, we show how to
#' # reverse the internal ordering of the pedigree.
#' m = as.matrix(x, include.attrs=TRUE)
#' m[] = m[3:1, ]
#'
#' # Must reverse the labels also:
#' attrs = attributes(m)
#' attrs$labels = rev(attrs$labels)
#'
#' # Restore ped:
#' y = restore_ped(m, attrs=attrs)
#'
#' # Of course a simpler way is use reorderPed():
#' z = reorderPed(x, 3:1)
#' stopifnot(identical(y, z))
#'
#' @export
as.matrix.ped = function(x, include.attrs = TRUE, ...) {
  m = cbind(x$ID, x$FID, x$MID, x$SEX)
  if(hasMarkers(x)) {
    markermatr = do.call(cbind, x$markerdata)
    m = cbind(m, markermatr)
  }
  if (include.attrs) {
    attr(m, "famid") = x$FAMID
    attr(m, "labels") = x$LABELS
    attr(m, "UNBROKEN_LOOPS") = x$UNBROKEN_LOOPS
    attr(m, "LOOP_BREAKERS") = x$LOOP_BREAKERS
    if(hasMarkers(x)) {
      attr(m, "markerattr") = lapply(x$markerdata, attributes)
    }
  }
  m
}


#' @rdname as.matrix.ped
#' @export
restore_ped = function(x, attrs = NULL, check = TRUE) {
  if (is.null(attrs))
    attrs = attributes(x)
  p = ped(id=x[,1], fid=x[,2], mid=x[,3], sex=x[,4], famid=attr(x, "famid"), check = check, reorder=F)
  p = setLabels(p, attrs$labels)
  p['LOOP_BREAKERS'] = list(attrs$LOOP_BREAKERS) # Trick to keep explicit NULLs

  ### Markers
  if((nc <- ncol(x)) > 4) {
    if(nc %% 2 != 0) stop2("Something is wrong: Odd number of allele columns!")
    markerattr = attrs$markerattr

    mlist = lapply(seq_len((nc-4)/2), function(k) {
      m = x[, c(3 + 2*k, 4 + 2*k), drop = F]
      attr = markerattr[[k]]
      attr$dim = dim(m)
      attributes(m) = attr
      m
    })
    class(mlist) = "markerList"
    p = setMarkers(p, mlist)
  }
  p
}

#' Convert ped to data.frame
#'
#' Convert a `ped` object to a data.frame. The first columns are id, fid, mid
#' and sex, followed by genotype columns for all (or a selection of) markers.
#'
#' Note that the output of [as.data.frame.ped()] is quite different from that of
#' [as.matrix.ped()]. This reflects the fact that these functions have different
#' purposes.
#'
#' Conversion to data.frame is primarily intended for pretty printing. It uses
#' correct labels for pedigree members and marker alleles, and pastes alleles to
#' form nice-looking genotypes.
#'
#' The matrix method, on the other hand, is a handy tool for manipulating the
#' pedigree structure. It produces a numeric matrix, using the internal index
#' labeling both for individuals and alleles, making it very fast. In addition,
#' all neccessary meta information (loop breakers, allele frequencies a.s.o) is
#' kept as attributes, which makes it possible to recreate the original `ped`
#' object.
#'
#' @param x Object of class `ped`.
#' @param ... Further parameters
#' @param markers (Optional) Vector of marker indices. By default, all markers
#'   are included.
#' @return A `data.frame` with `pedsize(x)` rows and `4 + nMarkers(x)` columns.
#' @seealso [as.matrix.ped()]
#'
#' @export
as.data.frame.ped = function(x, ..., markers) {
  lab = x$LABELS
  fid = mid = rep("0", pedsize(x))
  fid[x$FID > 0] = lab[x$FID]
  mid[x$MID > 0] = lab[x$MID]
  df = data.frame(id = lab, fid=fid, mid=mid, sex=x$SEX, stringsAsFactors=FALSE)

  if(hasMarkers(x)) {
    mlist = if(missing(markers)) x$markerdata else getMarkers(x, markers)
    geno = do.call(cbind, lapply(mlist, format))

    # headers of genotype columns: name if present, otherwise <idx>
    nms = vapply(mlist, name, character(1))
    if(any(na_name <- is.na(nms)))
      nms[na_name] = sprintf("<%d>", which(na_name))
    colnames(geno) = nms

    df = cbind(df, geno, stringsAsFactors=FALSE)
  }
  df
}

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
print.ped = function(x, ..., markers, verbose=TRUE) {
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
  datafr = as.data.frame(x, markers=markers, singleCol=TRUE, missing="-")
  datafr$fid[datafr$fid == "0"] = "*"
  datafr$mid[datafr$mid == "0"] = "*"
  print(datafr, row.names=FALSE, ...)

  if(showmess && verbose)
    message("Only 5 (out of ", nm, ") markers are shown. See `?print.ped` for options.")

  invisible(datafr)
}

