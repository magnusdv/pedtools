#' Plot pedigrees with genotypes
#'
#' This is the main function for pedigree plotting, with many options for
#' controlling the appearance of pedigree symbols and accompanying labels.
#' Most of the work is done by the plotting functionality in the 'kinship2'
#' package.
#'
#' \code{plot.ped} is in essence a wrapper for \code{plot.pedigree} in the
#' \code{kinship2} package.
#'
#' @param x a \code{\link{ped}} object.
#' @param id.labels a vector with labels for each pedigree member. This
#' defaults to \code{x$LABELS} (see \code{\link{setLabels}}).
#' @param text a character vector of length \code{x$NIND}.
#' @param title the plot title. If NULL or '', no title is added to the plot.
#' @param col a vector with color indicators for the pedigree members. Recycled
#' if necessary. By default everyone is drawn black.
#' @param deceased a numeric containing ID's of deceased pedigree members.
#' @param starred a numeric containing ID's of pedigree members that should be
#' marked with a star in the pedigree plot.
#' @param margins a numeric of length 4 indicating the plot margins. For
#' singletons only the first element (the 'bottom' margin) is used.
#' @param \dots arguments passed on to \code{plot.pedigree} in the
#' \code{kinship2} package. In particular \code{symbolsize} and \code{cex} can
#' be useful.
#' @author Magnus Dehli Vigeland, Guro Doerum
#' @seealso \code{\link{plot.pedigree}}, \code{\link{setLabels}}
#'
#' @examples
#'
#' x = cousinsPed(1)
#' plot(x)
#'
#' @export
plot.ped = function(x, id.labels = x$LABELS, text=NULL, title = NULL, col = 1, deceased = numeric(0),
                    starred = numeric(0), margins = c(0.6, 1, 4.1, 1), ...) {

  # Labels
  if (is.null(id.labels)) id.labels=rep("", pedSize(x))
  else if(identical(id.labels, "")) id.labels=rep("", pedSize(x))
  else if(identical(id.labels, "num")) id.labels = as.character(x$ID)

  id.labels[is.na(id.labels)] = ""


  #if (!is.null(lb <- x$loop_breakers)) {
  #    origint = lb[, 1])
  #    copyint = lb[, 2])
  #    id.labels[copyint] = paste(id.labels[copyint], id.labels[origint], sep = "=")
  #}

  strid = id.labels

  # Add stars to labels
  starred = internalID(x, starred)
  strid[starred] = paste0(strid[starred], "*")

  if(!is.null(text)) {
    stopifnot(length(text) == pedSize(x))
    strid = paste(strid, text, sep = "\n")
  }

  # Needed for centered title. Without, par() doesnt equal 'margins'...(why??)
  oldmar = par(mar = margins)

  # Colors
  cols = rep(col, length = pedSize(x))

  # Special treatment for option 'available=shaded'
  #if (identical(available, "shaded")) {
  #    if (any(c("angle", "density") %in% names(list(...))))
  #        stop("Plot parameters 'angle' and 'density' cannot be used in combination with 'available=shaded'")
  #    pedigree = as.kinship2_pedigree(x, deceased = deceased, aff2 = aff2)
  #    pdat = kinship2::plot.pedigree(pedigree, id = strid, col = cols, mar = margins, density = 25, angle = 45, ...)

  pedigree = as.kinship2_pedigree(x, deceased = deceased)
  pdat = kinship2::plot.pedigree(pedigree, id = strid, col = cols, mar = margins, ...)

  # Add title
  if (!is.null(title)) title(title)

  # par(oldmar)
  invisible(pdat)
}

#' @rdname plot.ped
#' @importFrom assertthat assert_that is.string
#' @export
plot.singleton = function(x, id.labels = x$LABELS, title = NULL, col = 1, deceased = numeric(0),
                          starred = numeric(0), margins = c(8, 0, 0, 0), ...) {
  assertthat::assert_that(is.null(id.labels) || is.string(id.labels))

  y = addParents(x, x$LABELS[1], verbose = FALSE) # reorder necessary??

  # If input id.labels is "num" or "" or something else than x$LABELS, pass it directly on.
  if(is.null(id.labels) || id.labels == "num") id = id.labels
  else id = c(id.labels, "", "")

  p = plot.ped(y, id.labels = id, title = title, col = col, deceased = numeric(0),
      starred = starred, margins = c(margins[1], 0, 0, 0), ...)
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




