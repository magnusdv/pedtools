#' Plot pedigree
#'
#' This is the main function for pedigree plotting, with many options for
#' controlling the appearance of pedigree symbols and accompanying labels. The
#' main pedigree layout is calculated with the `kinship2` package, see
#' [kinship2::align.pedigree] for details. Unlike `kinship2`, the implementation
#' here also supports singletons, and plotting pedigrees as DAGs. In addition,
#' some minor adjustments have been made to improve scaling and avoid unneeded
#' duplications.
#'
#' For an overview of all plotting options, see the separate page
#' [internalplot], which documents the internal plotting procedure in more
#' detail.
#'
#' @param x A [ped()] object.
#' @param draw A logical, by default TRUE. If FALSE, no plot is produced, only
#'   the plotting parameters are returned.
#' @param keep.par A logical, by default FALSE. If TRUE, the graphical
#'   parameters are not reset after plotting, which may be useful for adding
#'   additional annotation.
#' @param ... Arguments passed on to the internal plot functions. For a complete
#'   list of parameters, see [internalplot]. The most important ones are
#'   illustrated in the Examples below.
#'
#' @return A list of three lists with various plot details: `alignment`,
#'   `annotation`, `scaling`.
#'
#' @seealso [plotPedList()], [kinship2::plot.pedigree()]. Plot options are
#'   documented in [internalplot].
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
#' # Other options
#' plot(x, hatched = "boy", starred = "fa", deceased = "mo", title = "Fam 1")
#'
#' # Medical pedigree
#' plot(x, aff = "boy", carrier = "mo")
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
#' # By default, labels are trimmed for initial/trailing line breaks ...
#' plot(x, labs = c("\nFA" = "fa"))
#'
#' # ... but this can be overridden
#' plot(x, labs = c("\nFA" = "fa"), trimLabs = FALSE)
#'
#' # Colours
#' plot(x, col = c(fa = "red"), fill = c(mo = "green", boy = "blue"))
#'
#' # Non-black hatch colours are specified with the `fill` argument
#' plot(x, hatched = labels, fill = c(boy = "red"))
#'
#' # Use functions to specify colours
#' plot(x, fill = list(red = leaves, blue = ancestors(x, "boy")))
#'
#' # Line type and width
#' plot(x, lty = 2, lwd = 3, cex = 2)
#'
#' # Detailed line type and width
#' plot(x, lty = list(dashed = founders), lwd = c(boy = 4))
#'
#' # Include genotypes
#' x = addMarker(x, fa = "1/1", boy = "1/2", name = "SNP")
#' plot(x, marker = 1)
#'
#' # Markers can also be called by name
#' plot(x, marker = "SNP")
#'
#' # Plot as DAG (directed acyclic graph)
#' plot(x, arrows = TRUE, title = "DAG")
#'
#' # Founder inbreeding is shown by default
#' founderInbreeding(x, "mo") = 0.1
#' plot(x)
#'
#' # ... but can be suppressed
#' plot(x, fouInb = NULL)
#'
#' # Other text above and inside symbols
#' plot(x, textAbove = letters[1:3], textInside = LETTERS[1:3])
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
    on.exit(par(scaling$oldpar))

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
    on.exit(par(scaling$oldpar))

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
#' @inheritParams internalplot
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
#' H2 = list(singleton(1), singleton(3))  # grouped!
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
                       frames = TRUE, fmar = NULL, source = NULL,
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
      if (!is.ped(p[[1]]))
        stop2("First element must be a `ped` object", p[[1]])
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
    arglist$title = NULL

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
  on.exit(par(opar))

  if(verbose) {
    message("Group structure: ", toString(groups))
    message("Relative widths: ", toString(widths))
    message("Default margins: ", toString(marLR))
    message("Indiv. margins:")
    for(p in plotlist) message("  ", toString(p$margins))
    message("Input width/height: ", toString(c(dev.width, dev.height)))
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
    cex.title = extra.args$cex.main %||% NA
    mtext(finalTitles, outer = TRUE, at = midpoints, cex = cex.title, font = 2) # bold
  }
}


#' @rdname plot.ped
#' @export
#' @importFrom graphics plot.default
plot.list = function(x, ...) {

  if(!is.pedList(x))
    return(plot.default(x, ...))

  x1 = x[[1]]
  x2 = x[[2]]

  # Compute and merge plists
  p1 = .pedAlignment(x1, ...)$plist
  p2 = .pedAlignment(x2, ...)$plist

  r1 = length(p1$n)
  c1 = max(p1$n)
  r2 = length(p2$n)
  c2 = max(p2$n)

  nid = pos = fam = sp = matrix(0, nrow = max(r1, r2), ncol = c1 + c2)

  # n
  n1 = p1$n
  n2 = p2$n
  nid1 = p1$nid
  nid2 = p2$nid
  pos1 = p1$pos
  pos2 = p2$pos
  fam1 = p1$fam
  fam2 = p2$fam
  sp1 = p1$spouse
  sp2 = p2$spouse

  nInd1 = max(nid1)
  nInd2 = max(nid2)

  # Merge n
  n1pad = n2pad = integer(max(r1, r2))
  n1pad[seq_along(n1)] = n1
  n2pad[seq_along(n2)] = n2
  n = n1pad + n2pad
  n

  # Merge matrices
  for(i in seq_along(n)) {
    if(i <= r1) {
      seq1 = seq_len(n1[i])
      nid[i, seq1] = nid1[i, seq1]
      pos[i, seq1] = pos1[i, seq1]
      if(i>1) fam[i, seq1] = fam1[i, seq1]
      sp[i, seq1] = sp1[i, seq1]
    }
    if(i <= r2) {
      seq2 = seq_len(n2[i])
      nid[i, n1pad[i] + seq2] = nid2[i, seq2] + nInd1
      pos[i, n1pad[i] + seq2] = pos2[i, seq2] + max(pos1) + 1
      if(i>1) fam[i, n1pad[i] + seq2] = fam2[i, seq2] + n1pad[i-1] * (fam2[i, seq2] > 0)
      sp[i, n1pad[i] + seq2] = sp2[i, seq2]
    }
  }

  # Collect into plist
  plist = list(n = n, nid = nid, pos = pos, fam = fam, spouse = sp)
  align = .extendPlist(x, plist)

  # Adjust y positions
  #idx = match(6:8, plist2$plotord)
  #plist2$yall[idx] = plist2$yall[idx] + 0.5

  # Merge annotations
  annot1 = .pedAnnotation(x1, ...)
  annot2 = .pedAnnotation(x2, ...)

  mrg = function(a, def)
    c(annot1[[a]] %||% rep(def, nInd1), annot2[[a]] %||% rep(def, nInd2))

  annot = list(title = annot1$title,
               textUnder = mrg("textUnder", ""),
               textAbove = mrg("textAbove", ""),
               textInside = mrg("textInside", ""),
               colvec = mrg("colvec", 1),
               densvec = mrg("densvec", 0),
               fillvec =  mrg("fillvec", NA),
               ltyvec =  mrg("ltyvec", 1),
               lwdvec =  mrg("lwdvec", 1),
               density = annot1$density,
               angle = annot1$angle,
               carrierTF = mrg("carrierTF", FALSE),
               deceasedTF = mrg("deceasedTF", FALSE))

  drawPed(align, annot, scaling = NULL, ...)
}

