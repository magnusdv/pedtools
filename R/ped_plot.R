#' Plot pedigree
#'
#' This is the main function for plotting pedigrees. Many options are available
#' for controlling the appearance of pedigree symbols and accompanying labels.
#' The most important ones are illustrated in the Examples section below; for a
#' complete overview, see the separate page [plotmethods], which also explains
#' the plotting procedure in more detail.
#'
#' The main pedigree layout is calculated with the `kinship2` package, see
#' [kinship2::align.pedigree] for details. Unlike `kinship2`, the implementation
#' here also supports singletons, and plotting pedigrees as DAGs. In addition,
#' some minor adjustments have been made to improve scaling and avoid unneeded
#' duplications.
#'
#' If `x` is a list of `ped` objects, these are plotted next to each other,
#' vertically centred in the plotting window. For finer control, and possibly
#' nested lists of pedigrees, use [plotPedList()].
#'
#' @param x A [ped()] object or a list of such.
#' @param draw A logical, by default TRUE. If FALSE, no plot is produced, only
#'   the plotting parameters are returned.
#' @param keep.par A logical, by default FALSE. If TRUE, the graphical
#'   parameters are not reset after plotting, which may be useful for adding
#'   additional annotation.
#' @param ... Arguments passed on to the internal plot functions. For a complete
#'   list of parameters, see [plotmethods]. The most important ones are
#'   illustrated in the Examples below.
#'
#' @return A list of three lists with various plot details: `alignment`,
#'   `annotation`, `scaling`.
#'
#' @seealso [plotPedList()], [kinship2::plot.pedigree()]. Plot options are
#'   documented in [plotmethods].
#'
#' @examples
#'
#' # Singleton
#' plot(singleton(1))
#'
#' # Trio
#' x = nuclearPed(father = "fa", mother = "mo", child = "boy")
#' plot(x)
#'
#' #' # Modify margins
#' plot(x, margins = 6)
#' plot(x, margins = c(0,0,6,6)) # b,l,t,r
#'
#' # Larger text and symbols
#' plot(x, cex = 1.5)
#'
#' # Enlarge symbols only
#' plot(x, symbolsize = 1.5)
#'
#' # Various annotations
#' plot(x, hatched = "boy", starred = "fa", deceased = "mo", title = "Fam 1")
#'
#' # Swap spouse order
#' plot(x, spouseOrder = c("mo", "fa"))
#'
#' #----- ID labels -----
#'
#' # Label only some members
#' plot(x, labs = c("fa", "mo"))
#'
#' # Label males only
#' plot(x, labs = males)
#'
#' # Rename some individuals
#' plot(x, labs = c(FATHER = "fa", "boy"))
#'
#' # By default, long names are folded to width ~12 characters
#' plot(x, labs = c("Very long father's name" = "fa"), margin = 2)
#'
#' # Folding width may be adjusted ...
#' plot(x, labs = c("Very long father's name" = "fa"), foldLabs = 6)
#'
#' # ... or switched off (requires larger margin!)
#' plot(x, labs = c("Very long father's name" = "fa"), foldLabs = FALSE)
#'
#' # By default, labels are trimmed for initial/trailing line breaks ...
#' plot(x, labs = c("\nFA" = "fa"))
#'
#' # ... but this can be overridden
#' plot(x, labs = c("\nFA" = "fa"), trimLabs = FALSE)
#'
#' #----- Colours -----
#'
#' plot(x, col = c(fa = "red"), fill = c(mo = "green", boy = "blue"))
#'
#' # Non-black hatch colours are specified with the `fill` argument
#' plot(x, hatched = labels, fill = c(boy = "red"))
#'
#' # Use functions to specify colours
#' plot(x, fill = list(red = leaves, blue = ancestors(x, "boy")))
#'
#' #----- Symbol line types and widths -----
#'
#' # Dotted, thick symbols
#' plot(x, lty = 3, lwd = 4, cex = 2)
#'
#' # Detailed specification of line types and width
#' plot(x, lty = list(dashed = founders), lwd = c(boy = 4))
#'
#' #----- Genotypes -----
#'
#' x = nuclearPed(father = "fa", mother = "mo", child = "boy") |>
#'   addMarker(fa = "1/1", boy = "1/2", name = "SNP") |>
#'   addMarker(boy = "a/b")
#'
#' # Show genotypes for first marker
#' plot(x, marker = 1)
#'
#' # Show empty genotypes for untyped individuas
#' plot(x, marker = 1, showEmpty = TRUE)
#'
#' # Markers can also be called by name
#' plot(x, marker = "SNP")
#'
#' # Multiple markers
#' plot(x, marker = 1:2)
#'
#' #----- Further text annotation -----
#'
#' # Founder inbreeding is shown by default
#' xinb = x |> setFounderInbreeding("mo", value = 0.1)
#' plot(xinb)
#'
#' # ... but can be suppressed
#' plot(xinb, fouInb = NULL)
#'
#' # Text can be placed around and inside symbols
#' plot(x, textAnnot = list(topright = 1:3, inside = LETTERS[1:3]))
#'
#' # Use lists to add further options; see `?text()`
#' plot(x, margin = 2, textAnnot = list(
#'   topright = list(1:3, cex = 0.8, col = 2, font = 2, offset = 0.1),
#'   left = list(c(boy = "comment"), cex = 2, col = 4, offset = 2, srt = 20)))
#'
#' # Exhaustive list of annotation positions
#' plot(singleton(1), cex = 3, textAnnot = list(top="top", left="left",
#'   right="right", bottom="bottom", topleft="topleft", topright="topright",
#'   bottomleft="bottomleft", bottomright="bottomright", inside="inside"))
#'
#' #----- Special pedigrees -----
#'
#' # Plot as DAG (directed acyclic graph)
#' plot(x, arrows = TRUE, title = "DAG")
#'
#' # Medical pedigree
#' plot(x, aff = "boy", carrier = "mo")
#'
#' # Miscarriage
#' plot(x, miscarriage = "boy", deceased = "boy", labs = founders)
#'
#' # Twins
#' x = nuclearPed(children = c("tw1", "tw2", "tw3"))
#' plot(x, twins = data.frame(id1 = "tw1", id2 = "tw2", code = 1)) # MZ
#' plot(x, twins = data.frame(id1 = "tw1", id2 = "tw2", code = 2)) # DZ
#'
#' # Triplets
#' plot(x, twins = data.frame(id1 = c("tw1", "tw2"),
#'                            id2 = c("tw2", "tw3"),
#'                            code = 2))
#'
#' # Selfing
#' plot(selfingPed(2))
#'
#' # Complex pedigree: Quadruple half first cousins
#' plot(quadHalfFirstCousins())
#'
#' # Straight legs
#' plot(quadHalfFirstCousins(), align = c(0,0))
#'
#' # Lists of multiple pedigree
#' plot(list(singleton(1), nuclearPed(1), linearPed(2)))
#'
#' # Use of `drawPed()`
#' dat = plot(nuclearPed(), draw = FALSE)
#' drawPed(dat$alignment, dat$annotation, dat$scaling)
#'
#' @export
plot.ped = function(x, draw = TRUE, keep.par = FALSE, ...) {

  if(!(is.logical(draw) && length(draw) == 1))
    stop2("Illegal plot command.\n\nNote: As of version 2.0.0, plot arguments must be named, e.g. `plot(x, marker = ...)` to include genotypes")

  # Alignment parameters
  alignment = .pedAlignment(x, ...)

  # Prepare annotations
  annotation = .pedAnnotation(x, ...)

  # New plot window (do before scaling)
  frame()

  # Scaling parameters
  scaling = .pedScaling(alignment = alignment, annotation = annotation, ...)

  if(!keep.par)
    on.exit(par(scaling$oldpar), add = TRUE)

  if(draw) {
    # Symbols and connectors
    .drawPed(alignment, annotation, scaling)

    # Annotation
    .annotatePed(alignment = alignment, annotation = annotation, scaling = scaling, ...)
  }

  # Prepare output
  output = list(alignment = alignment, annotation = annotation, scaling = scaling)

  # For back compatibility
  output = c(output, alignment[c("plist", "x", "y")], scaling[c("boxw", "boxh")])

  invisible(output)
}



#' @rdname plot.ped
#'
#' @param alignment List of alignment details, as returned by [.pedAlignment()].
#' @param annotation List of annotation details as returned by [.pedAnnotation()].
#' @param scaling List of scaling parameters as returned by [.pedScaling()].

#' @export
drawPed = function(alignment, annotation = NULL, scaling = NULL, keep.par = FALSE, ...) {

  frame()

  if(is.null(scaling))
    scaling = .pedScaling(alignment, annotation, ...)

  if(!keep.par)
    on.exit(par(scaling$oldpar), add = TRUE)

  # Symbols and connectors
  .drawPed(alignment, annotation, scaling)

  # Annotation
  .annotatePed(alignment, annotation, scaling, ...)

  # Prepare output
  output = list(alignment = alignment, annotation = annotation, scaling = scaling)

  # For back compatibility
  output = c(output, alignment[c("plist", "x", "y")], scaling[c("boxw", "boxh")])

  invisible(output)
}



#' Convert pedigree to kinship2 format
#'
#' @inheritParams plotmethods
#'
#' @examples
#' x = nuclearPed()
#' as_kinship2_pedigree(x)
#'
#' @export
as_kinship2_pedigree = function(x, deceased = NULL, aff = NULL, twins = NULL, hints = NULL) {
    ped = as.data.frame(x)  # not as.matrix()

    if(nrow(ped) > 1) # fails (in kinship2::pedigree) for singletons
      ped$sex[ped$sex == 0] = 3 # kinship2 code for "diamond"

    affected = ifelse(ped$id %in% aff, 1, 0) # NULL => affected01 = c(0,0,..)
    status = ifelse(ped$id %in% deceased, 1, 0)

    arglist = list(id = ped$id, dadid = ped$fid, momid = ped$mid,
                   sex = ped$sex, affected = affected,
                   status = status, missid = 0)
    if(!is.null(twins))
      arglist$relation = twins

    # Avoid kinship2 warning about missing genders a.s.o.
    kinped = suppressWarnings(do.call(kinship2::pedigree, arglist))

    # Possible hints for kinship2::align.pedigree
    kinped$hints = hints

    kinped
}


#' @rdname plot.ped
#' @export
plot.pedList = function(x, ...) {
  plotPedList(x, frames = FALSE, ...)
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
#' @param plots A list of lists. Each element of `plots` is a list, where the
#'   first element is a pedigree, and the remaining elements are passed on to
#'   `plot.ped`. These elements must be correctly named. See examples below.
#' @param widths A numeric vector of relative widths of the subplots. Recycled
#'   to `length(plots)` if necessary, before passed on to [layout()]. Note that
#'   the vector does not need to sum to 1.
#' @param groups A list of vectors, each consisting of consecutive integers,
#'   indicating subplots to be grouped. By default the grouping follows the list
#'   structure of `plots`.
#' @param titles A character vector of titles for each group. Overrides titles
#'   given in individuals subplots.
#' @param grouptitlesArgs A list of arguments passed on to [mtext()] for
#'   `titles`.
#' @param frames A logical indicating if groups should be framed.
#' @param fmar A single number in the interval \eqn{[0,0.5)} controlling the
#'   position of the frames.
#' @param source NULL (default), or the name or index of an element of `plots`.
#'   If given, marker data is temporarily transferred from this to all the other
#'   pedigrees. This may save some typing when plotting the same genotypes on
#'   several pedigrees.
#' @param newdev A logical, indicating if a new plot window should be opened.
#' @param dev.height,dev.width The dimensions of the new plot window. If these
#'   are NA suitable values are guessed from the pedigree sizes.
#' @param verbose A logical.
#' @param \dots Further arguments passed on to each call to [plot.ped()].
#'
#' @author Magnus Dehli Vigeland
#'
#' @seealso [plot.ped()]
#'
#' @examples
#' ##################
#' # Basic examples #
#' ##################
#'
#' # Simples use: Just give a list of ped objects.
#' peds = list(nuclearPed(3), cousinPed(2), singleton(12), halfSibPed())
#' plotPedList(peds, newdev = TRUE)
#'
#' # Override automatic determination of relative widths
#' w = c(2, 3, 1, 2)
#' plotPedList(peds, widths = w)
#'
#' # In most cases the guessed dimensions are ok but not perfect.
#' # Resize plot window manually and re-plot with `newdev = FALSE` (default)
#' # plotPedList(peds, widths = w)
#'
#' ## Remove frames
#' plotPedList(peds, widths = w, frames = FALSE)
#'
#' # Non-default grouping
#' plotPedList(peds, widths = w, groups = list(1, 2:3, 4), titles = 1:3)
#'
#' # Parameters added in the main call are used in each sub-plot
#' plotPedList(peds, widths = w, labs = leaves, hatched = leaves,
#'             col = list(blue = males, red = females), symbolsize = 1.3)
#'
#' dev.off()
#'
#' #################################
#' # Example of automatic grouping #
#' #################################
#' H1 = nuclearPed()
#' H2 = singletons(id = c(1,3))
#'
#' plotPedList(list(H1, H2), dev.height = 3, dev.width = 4,
#'             titles = c(expression(H[1]), expression(H[2])),
#'             cex = 1.5, cex.main = 1.3)
#'
#' dev.off()
#'
#' ############################################################
#' # Complex example with individual parameters for each plot #
#' ############################################################
#'
#' # For more control of individual plots, each plot and all
#' # its parameters can be specified in its own list.
#'
#' x1 = nuclearPed(nch = 3) |>
#'   addMarker(`3` = "1/2")
#' plot1 = list(x1, title = "Plot 1", marker = 1, deceased = 1:2, cex = 1.3,
#'              margins = c(7, 4, 7, 4))
#'
#' x2 = cousinPed(2) |>
#'   addMarker(`11` = "A/A", `12` = "A/A")
#' plot2 = list(x2, title = "Family", marker = 1, symbolsize = 1.2, labs = NULL,
#'              margins = c(3, 4, 2, 4))
#'
#' x3 = singleton("NN")
#' plot3 = list(x3, cex = 2, carrier = "NN", lty = c(NN = 2))
#'
#' x4 = halfSibPed()
#' plot4 = list(x4, title = "Half sibs", cex = 1.3, hatched = leaves,
#'              col = list(red = founders), fill = list(blue = leaves),
#'              margins = c(7, 4, 7, 4))
#'
#' plotPedList(list(plot1, plot2, plot3, plot4), widths = c(2,3,1,2),
#'             fmar = 0.03, groups = list(1, 2:3, 4), newdev = TRUE,
#'             cex.main = 1.5)
#'
#' dev.off()
#'
#' ################################
#' # Example with large pedigrees #
#' ################################
#'
#' # Important to set device dimensions here
#'
#' plotPedList(list(halfCousinPed(4), cousinPed(7)),
#'             titles = c("Large", "Very large"), widths = c(1, 1.3),
#'             dev.height = 8, dev.width = 6, margins = 1.5)
#'
#' dev.off()
#'
#' @importFrom grDevices dev.new dev.size
#' @importFrom graphics grconvertX grconvertY layout mtext rect par plot
#' @export
plotPedList = function(plots, widths = NULL, groups = NULL, titles = NULL,
                       grouptitlesArgs = NULL, frames = TRUE, fmar = NULL, source = NULL,
                       dev.height = NULL, dev.width = NULL,
                       newdev = !is.null(dev.height) || !is.null(dev.width),
                       verbose = FALSE, ...) {

  if(!(isTRUE(frames) || isFALSE(frames))) {
    message("`frames` must be either TRUE or FALSE; use `groups` to specify framing groups")
    groups = frames
    frames = TRUE
  }

  # Check each entry for explicit plot args
  hasArgs = function(p) !is.ped(p) && !is.pedList(p)

  # If explicit source given, transfer marker data to all
  if(!is.null(source)) {
    srcPed = plots[[source]]
    if(hasArgs(srcPed))
      srcPed = srcPed[[1]]
    if(is.null(srcPed))
      stop2("Unknown source pedigree: ", source)
    if(nMarkers(srcPed) == 0)
      stop2("The source pedigree has no attached markers")
    plots = lapply(plots, function(p)
      if(hasArgs(p)) {p[[1]] = transferMarkers(from = srcPed, to = p[[1]]); p}
      else transferMarkers(from = srcPed, to = p))
  }

  deduceGroups = is.null(groups)
  if (deduceGroups) {
    groups = list()
    k = 0
  }

  # Flatten plot list
  flatlist = list()

  for (p in plots) {
    if (is.ped(p))
      newpeds = list(list(p))
    else if (is.pedList(p))
      newpeds = lapply(p, list)
    else { # if list of ped with plot arguments
      p1 = p[[1]]
      if(inherits(p1, "pedList"))
        class(p[[1]]) = "list"
      else if (!is.ped(p1))
        stop2("First element must be a `ped` object", p1)
      newpeds = list(p)
    }

    flatlist = c(flatlist, newpeds)
    if (deduceGroups) {
      groups = c(groups, list(k + seq_along(newpeds)))
      k = k + length(newpeds)
    }
  }

  N = length(flatlist)
  NG = length(groups)

  # Groups: check that each vector is consecutive integers, and no duplicates.
  for (v in groups)
    if (!is.numeric(v) || !isTRUE(all.equal.numeric(v, v[1]:v[length(v)])))
      stop2("Each element of `groups` must consist of consecutive integers: ", v)
  dup = anyDuplicated.default(unlist(groups))
  if (dup > 0)
    stop2("Plot occurring twice in `groups`: ", dup)

  # Group titles
  grouptitles = titles %||% names(plots)
  if (!is.null(grouptitles) && length(grouptitles) != NG)
    stop2(sprintf("Length of `titles` (%d) does not equal number of groups (%d)",
                  length(grouptitles), NG))

  # If no group titles, use inner titles if present
  finalTitles = grouptitles %||% sapply(flatlist, function(p) p$title %||% "")
  hasTitles = any(nchar(finalTitles) > 0)

  # Extract `plist$n` for each (contains width and #gen)
  nvec = lapply(flatlist, function(arglist) .pedAlignment(arglist[[1]], ...)$plist$n)

  # Max number of generations
  maxGen = max(lengths(nvec))

  # Max width
  maxWid = vapply(nvec, max, 1)

  extra.args = list(...)

  # Default left/right margin
  marLR = if (N > 2) 3 else 2

  # Prepare plot args for each ped
  plotlist = lapply(flatlist, function(arglist) {
    names(arglist)[1] = "x"

    # Additional arguments given in (...)
    for (parname in setdiff(names(extra.args), names(arglist)))
      arglist[[parname]] = extra.args[[parname]]

    # Remove inner title
    #arglist$title = NULL

    # Margins
    arglist$margins = arglist$margins %||% {
      g = generations(arglist$x)
      if(g == 1) # singleton
        rep(0, 4)
      else {
        marTB = 2*(maxGen - g) + 2
        c(marTB, marLR, marTB, marLR)
      }
    }
    arglist
  })

  # Relative plot widths
  if (is.null(widths))
    widths = sqrt(maxWid - 1) + 1
  else {
    if(!is.numeric(widths) && !length(widths) %in% c(1,N))
      stop2("`widths` must be a numeric of length either 1 or the total number of objects")
    widths = rep_len(widths, N)
  }

  # Layout of plot regions
  if (newdev) {
    dev.height = dev.height %||% {max(3, 1 * maxGen) + 0.3 * as.numeric(hasTitles)}
    dev.width = dev.width %||% {sum(maxWid) + 1}
    dev.new(height = dev.height, width = dev.width, noRStudioGD = TRUE)
  }

  new.oma = if (hasTitles) c(0, 0, 3, 0) else c(0, 0, 0, 0)
  opar = par(oma = new.oma, xpd = NA, mfrow = c(1,1), mar = c(0,0,0,0)) # include mfrow to ensure layout is reverted on exit
  on.exit(par(opar), add = TRUE)

  if(verbose) {
    message("Group structure: ", toString(groups))
    message("Relative widths: ", toString(widths))
    message("Default margins: ", toString(marLR))
    message("Indiv. margins:")
    for(p in plotlist) message("  ", toString(p$margins))
    message("Input width/height: ", toString(c(dev.width %||% NA, dev.height %||% NA)))
    message("Actual dimensions: ", toString(round(dev.size(),3)))
  }

  # Plot!
  layout(rbind(1:N), widths = widths)
  for (arglist in plotlist)
    do.call(plot, arglist)

  # Leftmost coordinate of each plot region (converted to value in [0,1]).
  ratios = c(0, cumsum(widths)/sum(widths))

  # Group coordinates
  grStartIdx = sapply(groups, function(v) v[1])
  grStopIdx = sapply(groups, function(v) v[length(v)])
  grStart = ratios[grStartIdx]
  grStop = ratios[grStopIdx + 1]  # since 0 first

  # Draw frames
  if(frames) {
    # Default margin: 3% of vertical height, but at most 0.25 inches.
    fmar = fmar %||% min(0.03, 0.25/dev.size()[2])

    margPix = grconvertY(0, from = "ndc", to = "device") * fmar
    margXnorm = grconvertX(margPix, from = "device", to = "ndc")
    frame_start = grconvertX(grStart + margXnorm, from = "ndc")
    frame_stop = grconvertX(grStop - margXnorm, from = "ndc")
    rect(xleft = frame_start,
         ybottom = grconvertY(1 - fmar, from = "ndc"),
         xright = frame_stop,
         ytop = grconvertY(fmar, from = "ndc"), xpd = NA)
  }

  # Add titles
  if(hasTitles) {
    midpoints = if(!is.null(grouptitles)) (grStart + grStop)/2 else ratios[1:N] + diff(ratios)/2
    gtArgs = list(text = finalTitles, at = midpoints, outer = TRUE, font = 2,
                cex = extra.args$cex.main %||% NA)
    if(!is.null(grouptitlesArgs))
      gtArgs = modifyList(gtArgs, grouptitlesArgs)
    #cex.title = extra.args$cex.main %||% NA
    #mtext(finalTitles, outer = TRUE, at = midpoints, cex = cex.title,line = 0, font = 2)
    do.call(mtext, gtArgs)
  }
}


#' @rdname plot.ped
#' @export
#' @importFrom graphics plot.default
plot.list = function(x, ...) {
  L = length(x)

  if(L == 0)
    stop2("Empty list, nothing to plot")

  # Variant of is.pedList() for finer control
  for(i in 1:L) {
    if(is.pedList(x[[i]]))
      stop2(sprintf("Nested ped list in component %d. This is not supported by `plot()`; try `plotPedList()` instead.", i))

    if(!is.ped(x[[i]]))
      return(plot.default(x, ...))
  }

  # Individual alignment data
  plists = lapply(x, function(y) .pedAlignment(y, ...)$plist)

  # Disallow arrows
  if(is.null(plists[[1]]$fam)) {
    stop2("`arrows = TRUE` is not supported for lists; use `plotPedList()` instead")
  }

  nlst = lapply(plists, function(pl) pl$n)
  nInd = vapply(nlst, sum, FUN.VALUE = 1)

  # Total number of generations (G)
  mx = max(lengths(nlst))
  G = 2 * mx - 1

  # Each pedigree uses only a subset of rows
  rws = lapply(nlst, function(ni) mx + seq(from = -length(ni)+1, to = length(ni)-1, by = 2))

  # Merge n (= # indivs in each grid row)
  n = integer(G)
  for(i in 1:L) {
    n[rws[[i]]] = n[rws[[i]]] + nlst[[i]]
  }

  # Create empty plist matrices
  nid = pos = fam = sp = matrix(0, nrow = G, ncol = max(n))

  # Storage for rightmost column used in each row (updated in each iter)
  maxcol.tmp = integer(G)

  # Storage for highest id used (updated in each iter)
  maxid.tmp = 0L

  for(i in 1:L) {
    rwi = rws[[i]]
    pl = plists[[i]]
    ni = pl$n
    nidi = pl$nid
    posi = pl$pos
    fami = pl$fam
    spi = pl$spouse
    nipad = `[<-`(integer(G), rwi, ni)

    # Matrix entry indices of ped i (original and in big)
    e = ebig = which(nidi > 0, arr.ind = TRUE)

    # Entry indices in merged matrix
    ebig[, 'row'] = rwi[ebig[, 'row']]
    ebig[, 'col'] = ebig[, 'col'] + maxcol.tmp[ebig[,'row']]

    # Merge nid
    nid[ebig] = nidi[e] + maxid.tmp

    # Merge pos
    pos[ebig] = posi[e] + max(pos) + 1

    # Merge fam: Only consider positive entries
    e.fam = e[fami[e] > 0, , drop = FALSE]
    ebig.fam = ebig[fami[e] > 0, , drop = FALSE]
    fam[ebig.fam] = fami[e.fam] + maxcol.tmp[ebig.fam[, "row"] - 2]

    # Merge spouse
    sp[ebig] = spi[e]

    # Update tmps
    maxcol.tmp = maxcol.tmp + nipad
    maxid.tmp = maxid.tmp + nInd[i]

  }

  # Collect into plist
  plist = list(n = n, nid = nid, pos = pos, fam = fam, spouse = sp)
  align = .extendPlist(x, plist)

  annotList = lapply(x, function(y) .pedAnnotation(y, ...))
  annot1 = annotList[[1]]

  # Merge vectors
  mrg = function(a, def) {
    vecs = lapply(1:L, function(i) annotList[[i]][[a]] %||% rep(def, nInd[i]))
    unlist(vecs)
  }

  # Merge lists
  mrgTxtGen = function() {
    res1 = annotList[[1]]$textAnnot
    for(nm in names(res1))
      res1[[nm]][[1]] = unlist(lapply(1:L, function(i) annotList[[i]]$textAnnot[[nm]][[1]]))

    res1
  }

  annot = list(title = annot1$title,
               textUnder = mrg("textUnder", ""),
               textAbove = mrg("textAbove", ""),
               textInside = mrg("textInside", ""),
               textAnnot = mrgTxtGen(),
               colvec = mrg("colvec", 1),
               densvec = mrg("densvec", 0),
               fillvec =  mrg("fillvec", NA),
               ltyvec =  mrg("ltyvec", 1),
               lwdvec =  mrg("lwdvec", 1),
               density = annot1$density,
               angle = annot1$angle,
               carrierTF = mrg("carrierTF", FALSE),
               deceasedTF = mrg("deceasedTF", FALSE),
               probandTF = mrg("probandTF", FALSE))

  drawPed(align, annot, scaling = NULL, vsep2 = TRUE, ...)
}

