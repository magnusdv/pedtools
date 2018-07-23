#' Plot pedigrees with genotypes
#'
#' This is the main function for pedigree plotting, with many options for
#' controlling the appearance of pedigree symbols and accompanying labels. Most
#' of the work is done by the plotting functionality in the 'kinship2' package.
#'
#' `plot.ped` is in essence a wrapper for `plot.pedigree` in the `kinship2`
#' package.
#'
#' @param x a [ped()] object.
#' @param marker either NULL, a vector of positive integers, a [`marker`]
#'   object, or a list of such. If NULL, no genotypes are plotted.  If a vector
#'   of integers is given, the corresponding marker objects are extracted from
#'   `x$markerdata`. The genotypes are written below each individual in the
#'   pedigree, in the format determined by `sep` and `missing`. See also
#'   `skip.empty.genotypes` below.
#' @param sep a character of length 1 separating alleles for diploid markers.
#' @param missing the symbol (integer or character) for missing alleles.
#' @param skip.empty.genotypes a logical. If TRUE, and `marker` is
#'   non-NULL, empty genotypes (which by default looks like '-/-') are not
#'   printed.
#' @param id.labels a vector with labels for each pedigree member. This defaults
#'   to `x$LABELS` (see [setLabels()]).
#' @param title the plot title. If NULL or '', no title is added to the plot.
#' @param col a vector with color indicators for the pedigree members. Recycled
#'   if necessary. By default everyone is drawn black.
#' @param deceased a numeric containing ID's of deceased pedigree members.
#' @param starred a numeric containing ID's of pedigree members that should be
#'   marked with a star in the pedigree plot.
#' @param margins a numeric of length 4 indicating the plot margins. For
#'   singletons only the first element (the 'bottom' margin) is used.
#' @param \dots arguments passed on to `plot.pedigree` in the `kinship2`
#'   package. In particular `symbolsize` and `cex` can be useful.
#' @author Magnus Dehli Vigeland, Guro Doerum
#' @seealso [plot.pedigree()], [setLabels()]
#'
#' @examples
#'
#' x = cousinsPed(1)
#' plot(x)
#'
#' @export
plot.ped = function(x, marker = NULL, sep = "/", missing = "-", skip.empty.genotypes = FALSE,
                    id.labels = x$LABELS, title = NULL, col = 1, deceased = numeric(0),
                    starred = numeric(0), margins = c(0.6, 1, 4.1, 1), ...) {

  # Labels
  if (is.null(id.labels)) id.labels=rep("", pedsize(x))
  else if(identical(id.labels, "")) id.labels=rep("", pedsize(x))
  else if(identical(id.labels, "num")) id.labels = as.character(x$ID)

  id.labels[is.na(id.labels)] = ""

  text = id.labels

  # Add stars to labels
  starred = internalID(x, starred)
  text[starred] = paste0(text[starred], "*")

  # Marker genotypes
  if (!is.null(marker)) {
    if (is.marker(marker))
      mlist = list(marker)
    else if (is.markerList(marker))
      mlist = marker
    else if (is.numeric(marker))
      mlist = getMarkers(x, markeridx=marker)
    else if (is.character(marker))
      mlist = getMarkers(x, markernames=marker)
    else
      stop("Argument `marker` must be either:\n",
           "  * an integer vector (of marker indices)\n",
           "  * a character vector (of marker names)\n",
           "  * a `marker` or `markerList` object", call.=FALSE)
    checkConsistency(x, mlist)

    gg = .prettyMarkers(mlist, sep = sep, missing = missing, singleCol = TRUE,
                        sex = x$pedigree[, "SEX"])
    geno = apply(gg, 1, paste, collapse = "\n")
    if (skip.empty.genotypes)
      geno[rowSums(do.call(cbind, mlist)) == 0] = ""

    text = if (!any(nzchar(text))) geno else paste(text, geno, sep = "\n")
  }

  # Needed for centered title. Without, par() doesnt equal 'margins'...(why??)
  oldmar = par(mar = margins)

  # Colors
  cols = rep(col, length = pedsize(x))

  # Special treatment for option 'available=shaded'
  #if (identical(available, "shaded")) {
  #    if (any(c("angle", "density") %in% names(list(...))))
  #        stop("Plot parameters 'angle' and 'density' cannot be used in combination with 'available=shaded'")
  #    pedigree = as.kinship2_pedigree(x, deceased = deceased, aff2 = aff2)
  #    pdat = kinship2::plot.pedigree(pedigree, id = text, col = cols, mar = margins, density = 25, angle = 45, ...)

  pedigree = as.kinship2_pedigree(x, deceased = deceased)
  pdat = kinship2::plot.pedigree(pedigree, id = text, col = cols, mar = margins, ...)

  # Add title
  if (!is.null(title)) title(title)

  # par(oldmar)
  invisible(pdat)
}

#' @rdname plot.ped
#' @export
plot.singleton = function(x, marker = NULL, sep = "/", missing = "-", skip.empty.genotypes = FALSE,
                          id.labels = x$LABELS, title = NULL, col = 1, deceased = numeric(0),
                          starred = numeric(0), margins = c(8, 0, 0, 0), ...) {
  assert_that(is.null(id.labels) || is.string(id.labels))

  y = addParents(x, x$LABELS[1], verbose = FALSE) # reorder necessary??

  # If input id.labels is "num" or "" or something else than x$LABELS, pass it directly on.
  if(is.null(id.labels) || id.labels == "num") id = id.labels
  else id = c(id.labels, "", "")

  p = plot.ped(y, marker = marker, sep =sep, missing = missing,
               skip.empty.genotypes = skip.empty.genotypes, id.labels = id,
               title = title, col = col, deceased = numeric(0), starred = starred,
               margins = c(margins[1], 0, 0, 0), ...)

  usr = par("usr")
  rect(usr[1] - 0.1, p$y[1], usr[2] + 0.1, usr[4], border = NA, col = "white")

  if (!is.null(title)) title(title, line = -2.8)
}

#' @rdname plot.ped
#' @export
as.kinship2_pedigree = function(x, deceased = numeric(0)) {
    ped = as.data.frame(x) # to include original labels

    status = ifelse(ped$id %in% deceased, 1, 0)
    kinship2::pedigree(id = ped$id, dadid = ped$fid, momid = ped$mid, sex = ped$sex,
        status = status, missid=0)
}




