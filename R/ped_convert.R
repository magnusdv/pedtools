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
#' * `FAMID` the family identifier (a string)
#' * `LABELS` the ID labels (a character vector)
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
#' attrs$LABELS = rev(attrs$LABELS)
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
    attr(m, "FAMID") = famid(x)
    attr(m, "LABELS") = labels(x)
    attr(m, "UNBROKEN_LOOPS") = has_unbroken_loops(x)
    attr(m, "LOOP_BREAKERS") = x$LOOP_BREAKERS
    attr(m, "FOUNDER_INBREEDING") =
      if(is.null(x$FOUNDER_INBREEDING)) NULL
      else founder_inbreeding(x, named = TRUE)
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
  p = ped(id=x[,1], fid=x[,2], mid=x[,3], sex=x[,4], famid=attrs$FAMID,
          check = check, reorder=F)
  p = relabel(p, new = attrs$LABELS)
  p['LOOP_BREAKERS'] = list(attrs$LOOP_BREAKERS) # Trick to keep explicit NULLs

  # Founder inbreeding
  finb = attrs$FOUNDER_INBREEDING
  if(is.null(finb))
    p['FOUNDER_INBREEDING'] = list(NULL)
  else {
    new_fou = founders(p)
    no_longer_fou = .mysetdiff(names(finb), new_fou)
    if(any(finb[no_longer_fou] > 0))
      stop2("Individual with nonzero founder inbreeding is no longer a founder: ", no_longer_fou)
    finb = finb[intersect(names(finb), new_fou)]
    founder_inbreeding(p) = finb
  }

  ### Markers
  if((nc <- ncol(x)) > 4) {
    if(nc %% 2 != 0) stop2("Something is wrong: Odd number of allele columns!")
    markerattr = attrs$markerattr
    pedlabs = labels(p)
    
    mlist = lapply(seq_len((nc-4)/2), function(k) {
      m = x[, c(3 + 2*k, 4 + 2*k), drop = F]
      attr = markerattr[[k]]
      attr$dim = dim(m)
      attr$pedmembers = pedlabs
      attr$sex = p$SEX
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
  lab = labels(x)
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

#' Conversions to ped objects
#'
#' @param x Any object.
#' @param ... Not used.
#'
#' @return A `ped` object or a list of such.
#'
#' @examples
#' df = data.frame(famid = c("S1", "S2"),
#'                 id = c("A", "B"),
#'                 fid = 0,
#'                 mid = 0,
#'                 sex = 1)
#'
#' # gives a list of two singletons
#' as.ped(df)
#'
#' @export
as.ped = function(x, ...) {
  UseMethod("as.ped")
}

#' @param famid_col Index of family ID column. If NA, the program looks for a
#'   column named "famid" (ignoring case).
#' @param id_col Index of individual ID column. If NA, the program looks for a
#'   column named "id" (ignoring case).
#' @param fid_col Index of father ID column. If NA, the program looks for a
#'   column named "fid" (ignoring case).
#' @param mid_col Index of mother ID column. If NA, the program looks for a
#'   column named "mid" (ignoring case).
#' @param sex_col Index of column with gender codes (0 = unknown; 1 = male; 2 =
#'   female). If NA, the program looks for a column named "sex" (ignoring case).
#'   If this is not found, genders of parents are deduced from the data, leaving
#'   the remaining as unknown.
#' @param marker_col Index vector indicating columns with marker alleles. If NA,
#'   all columns to the right of all pedigree columns are used. If `allele_sep`
#'   (see below) is non-NULL, each column is interpreted as a genotype column
#'   and split into separate alleles with `strsplit(..., split=allele_sep)`.
#' @param locus_annotations Passed on to [setMarkers()] (see explanation there).
#' @param missing Passed on to [setMarkers()] (see explanation there).
#' @param allele_sep Passed on to [setMarkers()] (see explanation there).
#' @param validate A logical indicating if the pedigree structure should be validated.
#'
#' @rdname as.ped
#' @export
as.ped.data.frame = function(x, famid_col=NA, id_col=NA, fid_col=NA, mid_col=NA, sex_col=NA,
                             marker_col=NA, locus_annotations=NULL, missing=0, allele_sep=NULL,
                             validate = TRUE, ...) {

  colnames = tolower(names(x))
  if(is.na(famid_col))
    famid_col = match("famid", colnames)

  if(is.na(famid_col)) {
    multiple_fams = FALSE
  }
  else {
    famid = x[[famid_col]]
    unique_fams = sort.default(unique.default(famid))
    multiple_fams = length(unique_fams) > 1
  }

  # If disconected, treat each component separatly by recursion
  if(multiple_fams) {
    pedlist = lapply(unique_fams, function(fam) {
      comp = x[famid == fam, , drop=FALSE]

      as.ped.data.frame(comp, famid_col = famid_col, id_col = id_col, fid_col = fid_col,
                        mid_col = mid_col, sex_col = sex_col, marker_col = marker_col,
                        locus_annotations = locus_annotations, missing = missing,
                        allele_sep=allele_sep, validate = validate, ...)
    })

    names(pedlist) = unique_fams
    return(pedlist)
  }

  #####################################################
  ### Body of function (for single ped) starts here ###
  #####################################################
  # By this stage, `famid_col` is NA or points to column with identical entries

  if(is.na(id_col)) id_col = match("id", colnames)
  if(is.na(fid_col)) fid_col = match("fid", colnames)
  if(is.na(mid_col)) mid_col = match("mid", colnames)
  if(is.na(sex_col)) sex_col = match("sex", colnames)
  # famid_col has already been dealt with

  ### Various checks
  NC = ncol(x)
  col_idx = c(famid=famid_col, id=id_col, fid=fid_col, mid=mid_col, sex=sex_col)

  # id, fid, mid cannot be missing
  required = col_idx[2:4]
  if(anyNA(required))
    stop2("Cannot find required column: ", names(required)[is.na(required)])

  # Catch duplicated column indices
  dup_idx = anyDuplicated(col_idx, incomparables = NA)
  if(dup_idx > 0)
    stop2("Column ", dup_idx, " has mulitple assignments: ", names(col_idx)[col_idx == dup_idx])

  # Chech that columns exist
  nonexist = !is.na(col_idx) & (col_idx < 1 | col_idx > NC)
  if(any(nonexist))
    stop2("Column index out of range: ", col_idx[nonexist])

  ### Get famid string
  famid = if(is.na(famid_col)) "" else x[[famid_col]][1]

  ### If sex is missing, deduce partially from parental status
  if(is.na(sex_col)) {
    sex = integer(nrow(x))
    sex[match(fid, id)] = 1
    sex[match(mid, id)] = 2
    bisex = intersect(fid, mid)
    sex[match(bisex, id)] = 0
  }
  else
    sex = x[[sex_col]]


  ### Create ped
  id = x[[id_col]]
  fid = x[[fid_col]]
  mid = x[[mid_col]]

  p = ped(id = id, fid = fid, mid = mid, sex = sex, famid = famid, check = validate, reorder = FALSE)

  ### Marker columns
  last_pedcol = max(col_idx, na.rm = TRUE)
  if(is.na(marker_col) && NC > last_pedcol)
    marker_col = seq.int(last_pedcol + 1, NC)
  else
    marker_col = NULL

  if(length(marker_col) > 0) {
    p = setMarkers(p, allele_matrix = x[marker_col],
                   locus_annotations = locus_annotations,
                   missing = missing, allele_sep = allele_sep)
  }

  p
}

