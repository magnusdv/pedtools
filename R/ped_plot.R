#' Plot pedigrees with genotypes
#'
#' This is the main function for pedigree plotting, with many options for
#' controlling the appearance of pedigree symbols and accompanying labels. Most
#' of the work is done by the plotting functionality in the 'kinship2' package.
#'
#' `plot.ped` is in essence an elaborate wrapper for
#' [kinship2::plot.pedigree()].
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
#' @param skip.empty.genotypes a logical. If TRUE, and `marker` is non-NULL,
#'   empty genotypes (which by default looks like '-/-') are not printed.
#' @param id.labels a vector with labels for each pedigree member. This defaults
#'   to `labels(x)`. Several syntaxes are possible:
#'
#'   * If `id.labels` is NULL or the empty character "", then no labels are
#'   drawn.
#'
#'   * If `id.labels` is the word "num", then all individuals are numerically
#'   labelled following the internal ordering.
#'
#'   * If `id.labels` is a subset of `labels(x)`, then only this subset will be
#'   labelled. If the vector is named, then the (non-empty) names are used
#'   instead of the ID label. See Examples.
#'
#'
#' @param title the plot title. If NULL or '', no title is added to the plot.
#' @param col a vector of colors for the pedigree members, recycled if
#'   necessary. Alternatively, `col` can be a list assigning colors to specific
#'   members. For example if `col = list(red = "foo", blue = c("bar", "baz"))`
#'   then "foo" will be red, "bar" and "baz" blue, and everyone else black. By
#'   default everyone is drawn black.
#' @param shaded a vector of ID labels indicating pedigree members whose plot
#'   symbols should appear shaded.
#' @param deceased a vector of ID labels indicating deceased pedigree members.
#' @param starred a vector of ID labels indicating pedigree members that should
#'   be marked with a star in the pedigree plot.
#' @param margins a numeric of length 4 indicating the plot margins. For
#'   singletons only the first element (the 'bottom' margin) is used.
#' @param \dots arguments passed on to `plot.pedigree` in the `kinship2`
#'   package. In particular `symbolsize` and `cex` can be useful.
#' @author Magnus Dehli Vigeland
#' @seealso [kinship2::plot.pedigree()]
#'
#' @examples
#'
#' x = nuclearPed(father = "fa", mother = "mo", child = "boy")
#' m = marker(x, fa = 1, boy = 1:2, name="SNP")
#'
#' plot(x, marker = m)
#'
#' # Alternative syntax if the marker is attached to x
#' x = setMarkers(x, m)
#' plot(x, marker = "SNP")
#'
#' # Other options
#' plot(x, marker = "SNP", shaded = typedMembers(x),
#'      starred = "fa", deceased = "mo")
#'
#' # Labelling only some members
#' plot(x, id.labels = c("fa", "boy"))
#'
#' # Labelling only some members, and renaming the father
#' plot(x, id.labels = c(FATHER = "fa", "boy"))
#'
#' # Colors
#' plot(x, col = list(red = "fa", green = "boy"))
#'
#' @export
plot.ped = function(x, marker = NULL, sep = "/", missing = "-", skip.empty.genotypes = FALSE,
                    id.labels = labels(x), title = NULL, col = 1, shaded = NULL, deceased = NULL,
                    starred = NULL, margins = c(0.6, 1, 4.1, 1), ...) {

  nInd = pedsize(x)

  # Labels
  nms = names(id.labels)
  if (is.null(id.labels) || identical(id.labels, ""))
    id.labels = rep("", nInd)
  else if(identical(id.labels, "num"))
    id.labels = as.character(1:nInd)
  else if(length(id.labels) < nInd && is.null(nms)) {
    labs = rep("", nInd)
    labs[internalID(x, id.labels)] = id.labels
    id.labels = labs
  }
  else if(!is.null(nms)) {
    labs = rep("", nInd)
    int_ids = internalID(x, id.labels)
    labs[int_ids] = id.labels
    # Replace with names where non-trivial
    hasName = (nms != "") & (!is.na(nms))
    labs[int_ids[hasName]] = nms[hasName]
    id.labels = labs
  }
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
    else if (is.numeric(marker) || is.character(marker))
      mlist = getMarkers(x, markers=marker)
    else
      stop2("Argument `marker` must be either:\n",
           "  * a n object of class `marker`\n",
           "  * a list of `marker` objects\n",
           "  * a character vector (names of attached markers)\n",
           "  * an integer vector (indices of attached markers)")
    checkConsistency(x, mlist)

    gg = do.call(cbind, lapply(mlist, format, sep = sep, missing = missing))
    geno = apply(gg, 1, paste, collapse = "\n")
    if (skip.empty.genotypes)
      geno[rowSums(do.call(cbind, mlist)) == 0] = ""

    text = if (!any(nzchar(text))) geno else paste(text, geno, sep = "\n")
  }

  # Needed for centered title. Without, par() doesnt equal 'margins'...(why??)
  oldmar = par(mar = margins)

  # Colors
  if(is.list(col)) {
    colnames = names(col)
    cols = rep(1, nInd)
    for(cc in colnames)
      cols[internalID(x, col[[cc]])] = cc
  }
  else {
     cols = rep(col, length = nInd)
  }

  # Shading
  if (!is.null(shaded)) {
    density = 25
    angle = 45
  } else {
    density = NULL
    angle = NULL
  }

  pedigree = as_kinship2_pedigree(x, deceased = deceased, shaded = shaded)
  pdat = kinship2::plot.pedigree(pedigree, id = text, col = cols, mar = margins,
                                 density = density, angle = angle, ...)

  # Add title
  if (!is.null(title)) title(title)

  # par(oldmar)
  invisible(pdat)
}

#' @rdname plot.ped
#' @export
plot.singleton = function(x, marker = NULL, sep = "/", missing = "-", skip.empty.genotypes = FALSE,
                          id.labels = labels(x), title = NULL, col = 1, shaded = NULL, deceased = NULL,
                          starred = NULL, margins = c(8, 0, 0, 0), ...) {
  if(length(id.labels) > 1)
    stop2("Argument `id.labels` must have length 1 in singleton plot: ", id.labels)

  y = addParents(x, labels(x)[1], verbose = FALSE)

  # Marker genotypes
  if (!is.null(marker)) {
    if (is.marker(marker))
      mlist = list(marker)
    else if (is.markerList(marker))
      mlist = marker
    else if (is.numeric(marker) || is.character(marker))
      mlist = getMarkers(x, markers = marker)
    else
      stop2("Argument `marker` must be either:\n",
           "  * an object of class `marker`\n",
           "  * a list of `marker` objects\n",
           "  * a character vector (names of attached markers)\n",
           "  * an integer vector (indices of attached markers)")
    checkConsistency(x, mlist)

    y = transferMarkers(setMarkers(x, mlist), y)
  }

  # Tweak id.labels if necessary. After addParents, internal index is 3!
  if(length(id.labels) == 0 || id.labels == "")
    id = NULL
  else if(id.labels == "num")
    id = c(`1` = labels(x))
  else if(is.null(names(id.labels))) {
    id = labels(x)
    names(id) = id.labels
  }
  else if(!is.null(names(id.labels))) {
    id = id.labels
  }

  p = plot.ped(y, marker = y$markerdata, sep = sep, missing = missing,
               skip.empty.genotypes = skip.empty.genotypes, id.labels = id,
               title = title, col = col, shaded = shaded, deceased = deceased, starred = starred,
               margins = c(margins[1], 0, 0, 0), ...)

  usr = par("usr")
  rect(usr[1] - 0.1, p$y[3], usr[2] + 0.1, usr[4], border = NA, col = "white")

  if (!is.null(title)) title(title, line = -2.8)
}

#' @rdname plot.ped
#' @export
as_kinship2_pedigree = function(x, deceased = NULL, shaded = NULL) {
    ped = as.data.frame(x)  # not as.matrix()
    ped$sex[ped$sex == 0] = 3 # kinship2 code for "diamond"

    affected = status = NULL

    if(!is.null(shaded))
      affected = ifelse(ped$id %in% shaded, 1, 0)

    if(!is.null(deceased))
      status = ifelse(ped$id %in% deceased, 1, 0)

    suppressWarnings(kinship2::pedigree(id = ped$id, dadid = ped$fid,
                                        momid = ped$mid, sex = ped$sex,
                                        affected = affected, status = status,
                                        missid=0))
}

#' @rdname plot.ped
#' @export
plot.pedList = function(x, frames = F, ...) {
  plotPedList(x, frames = frames, ...)
}

#' Plot a collection of pedigrees.
#'
#' This function creates a row of pedigree plots, each created by [plot.ped()].
#' Any parameter accepted by [plot.ped()] can be applied, either to all plots
#' simultaneously, or to individual plots.  Some effort is made to guess a
#' reasonable window size and margins, but in general the user must be prepared
#' to do manual resizing of the plot window. See various examples in the
#' Examples section below.
#'
#' Note that for tweaking dev.height and dev.width the function [dev.size()] is
#' useful to determine the size of the active device.
#'
#' @param plot.arg.list A list of lists. Each element of `plot.arg.list` is a
#'   list, where the first element is the [ped()] object to be plotted, and the
#'   remaining elements are passed on to `plot.ped`. These elements must be
#'   correctly named. See examples below.
#' @param widths A numeric vector of relative widths of the subplots. Recycled
#'   to `length(plot.arg.list)` if necessary, before passed on to [layout()].
#'   Note that the vector does not need to sum to 1.
#' @param frames Either a single logical (FALSE = no frames; TRUE = automatic
#'   framing) or a list of numeric vectors: Each vector must consist of
#'   consecutive integers, indicating subplots to be framed together. By default
#'   the framing follows the list structure of `plot.arg.list`.
#' @param frametitles A character vector of titles for each frame. If this is
#'   non-NULL, titles for individuals subplots are ignored.
#' @param fmar A single number in the interval \eqn{[0,0.5)} controlling the
#'   position of the frames.
#' @param newdev A logical, indicating if a new plot window should be opened.
#' @param dev.height,dev.width The dimensions of the new device (only relevant
#'   if newdev is TRUE). If these are NA suitable values are guessed from the
#'   pedigree sizes.
#' @param \dots Further arguments passed on to each call to [plot.ped()].
#' @author Magnus Dehli Vigeland
#' @seealso [plot.ped()]
#' @examples
#' # Simplest use: Just give a list of ped objects.
#' # To guess suitable plot window dimensions, use 'newdev=T'
#' peds = list(nuclearPed(3), cousinPed(2), singleton(12), halfSibPed())
#' plotPedList(peds, newdev=TRUE)
#'
#' # Modify the relative widths (which are not guessed)
#' widths = c(2, 3, 1, 2)
#' plotPedList(peds, widths=widths)
#'
#' # In most cases the guessed dimensions are ok but not perfect.
#' # Resize plot window manually, and then plot again with newdev=F (default)
#' # plotPedList(peds, widths=widths)
#'
#' ## Remove frames
#' plotPedList(peds, widths=widths, frames=FALSE)
#'
#' # Non-default frames
#' frames = list(1, 2:3)
#' plotPedList(peds, widths=widths, frames=frames, frametitles=c('First', 'Second'))
#'
#' # To give *the same* parameter to all plots, it can just be added at the end:
#' margins=c(2,4,2,4)
#' title='Same title'
#' id.labels=''
#' symbolsize=1.5
#' plotPedList(peds, widths=widths, frames=frames, margins=margins, title=title,
#'             id.labels=id.labels, symbolsize=symbolsize, newdev=TRUE)
#'
#' \dontrun{
#' # COMPLEX EXAMPLE WITH MARKER DATA AND VARIOUS OPTIONS
#' # For more control of individual plots, each plot and all its parameters
#' # can be specified in its own list:
#' x1 = nuclearPed(nch = 3)
#' m1 = marker(x1, '3' = 1:2)
#' marg1 = c(7, 4, 7, 4)
#' plot1 = list(x1, marker=m1, margins=marg1, title='Plot 1', deceased=1:2, cex=1.3)
#'
#' x2 = cousinPed(2)
#' m2 = marker(x2, alleles='A')
#' genotype(m2, leaves(x2)) = 'A'
#' marg2 = c(3, 4, 2, 4)
#' plot2 = list(x2, marker=m2, margins=marg2, title='Plot 2', symbolsize=1.2,
#'              skip.empty.genotypes=T, id=NULL)
#'
#' x3 = singleton("Mr. X")
#' marg3 = c(10, 0, 0, 0)
#' plot3 = list(x3, margins=marg3, title='Plot 3', symbolsize=1, cex=2)
#'
#' x4 = halfSibPed()
#' shaded = 4:5
#' col = c("black", "black", "black", "blue", "blue")
#' marg4 = marg1
#' plot4 = list(x4, margins=marg4, title='Plot 4', shaded=shaded, col=col)
#'
#' plotPedList(list(plot1, plot2, plot3, plot4), widths=c(2,3,1,2),
#'             frames=list(1,2:3,4), newdev=T)
#'
#' # Different example:
#' plotPedList(list(halfCousinPed(4), cousinPed(7)), title='Many generations',
#'     new=T, dev.height=9, dev.width=9)
#' }
#'
#' @importFrom grDevices dev.new dev.size
#' @importFrom graphics grconvertX grconvertY layout mtext rect par plot
#' @export
plotPedList = function(plot.arg.list, widths = NA, frames = T, frametitles = NULL, fmar = NA,
                       newdev = F, dev.height = NA, dev.width = NA, ...) {

  plot.list.flattened = list()
  if (deduceFrames <- isTRUE(frames)) {
    frames = list()
    k = 0
  }
  for (p in plot.arg.list) {
    if (is.ped(p))
      p = list(p)  # will now be included in next line
    if (is.pedList(p)) {
      plot.list.flattened = c(plot.list.flattened, lapply(p, list))
    }
    else {
        # if list of ped with plot arguments
        if (!is.ped(p[[1]]))
          stop2("First element must be a `ped` object", p[[1]])
        p = list(p)
        plot.list.flattened = append(plot.list.flattened, p)
      }
    if (deduceFrames) {
      group = (k + 1):(k <- k + length(p))
      frames = append(frames, list(group))
    }
  }
  plot.arg.list = plot.list.flattened
  N = length(plot.arg.list)
  if (identical(widths, NA))
    widths = vapply(plot.arg.list, function(p) ifelse(is.singleton(p[[1]]), 1, 2.5), 1)
  else
    widths = rep_len(widths, N)
  maxGen = max(vapply(plot.arg.list, function(arglist) .generations(arglist[[1]]), 1))

  if (hasframetitles <- !is.null(frametitles))
    if(length(frametitles) != length(frames))
      stop2(sprintf("Length of `frametitles` (%d) does not equal number of frames (%d)",
            length(frametitles), length(frames)))

  extra.args = list(...)
  if (!"title" %in% names(extra.args))
    extra.args$title = ""

  defaultmargins = if (N > 2)
    c(0, 4, 0, 4) else c(0, 2, 0, 2)

  plot.arg.list = lapply(plot.arg.list, function(arglist) {
    names(arglist)[1] = "x"
    g = .generations(arglist$x)
    addMargin = 2 * (maxGen - g + 1)
    if (!"margins" %in% c(names(arglist), names(extra.args)))
      arglist$margins = defaultmargins + c(addMargin, 0, addMargin, 0)

    # additional arguments given in (...)
    for (parname in setdiff(names(extra.args), names(arglist)))
      arglist[[parname]] = extra.args[[parname]]
    arglist
  })

  # title: this must be treated specially (in outer margins)
  titles = sapply(plot.arg.list, "[[", "title")
  plot.arg.list = lapply(plot.arg.list, function(arglist) {
    arglist$title = ""
    arglist
  })
  hastitles = hasframetitles || any(titles != "")

  # frame list: check that each vector is consecutive integers, and no duplicates.
  if (is.list(frames)) {
    for (v in frames)
      if (!identical(TRUE, all.equal(v, v[1]:v[length(v)])))
        stop2("Each element of `frames` must consist of consecutive integers: ", v)
    dup = anyDuplicated(unlist(frames))
    if (dup > 0)
      stop2("Plot occurring twice in `frames` list: ", dup)
  }


  # create layout of plot regions and plot!
  if (newdev) {
    if (is.na(dev.height))
      dev.height = max(3, 1 * maxGen) + 0.3 * as.numeric(hastitles)
    if (is.na(dev.width))
      dev.width = 3 * N
    dev.new(height = dev.height, width = dev.width, noRStudioGD = TRUE)
  }

  if (hastitles)
    par(oma = c(0, 0, 3, 0), xpd = NA)
  else
    par(oma = c(0, 0, 0, 0), xpd = NA)

  layout(rbind(1:N), widths = widths)
  for (arglist in plot.arg.list)
    do.call(plot, arglist)

  # leftmost coordinate of each plot region (converted to value in [0,1]).
  ratios = c(0, cumsum(widths)/sum(widths))

  # add frames
  if (is.list(frames)) {
    midpoints = numeric()
    fstart_index = sapply(frames, function(v) v[1])
    fstop_index = sapply(frames, function(v) v[length(v)])
    ratio_start = ratios[fstart_index]
    ratio_stop = ratios[fstop_index + 1]  # fordi 0 foerst
    midpoints = (ratio_start + ratio_stop)/2

    # margin (fmar): if NA, set to 5% of vertical height, but at most 0.25 inches.
    if (is.na(fmar))
      fmar = min(0.05, 0.25/dev.size()[2])
    margPix = grconvertY(0, from = "ndc", to = "device") * fmar
    margXnorm = grconvertX(margPix, from = "device", to = "ndc")
    frame_start = grconvertX(ratio_start + margXnorm, from = "ndc")
    frame_stop = grconvertX(ratio_stop - margXnorm, from = "ndc")
    rect(xleft = frame_start,
         ybottom = grconvertY(1 - fmar, from = "ndc"),
         xright = frame_stop,
         ytop = grconvertY(fmar, from = "ndc"))
  }

  cex.title =
    if ("cex.main" %in% names(extra.args)) extra.args$cex.main
    else NA

  if (hasframetitles) {
    for (i in 1:length(frames))
      mtext(frametitles[i], outer = TRUE, at = midpoints[i], cex = cex.title)
  }
  else if (hastitles) {
    for (i in 1:N)
      mtext(titles[i], outer = TRUE, at = (ratios[i] + ratios[i + 1])/2, cex = cex.title)
  }
}

