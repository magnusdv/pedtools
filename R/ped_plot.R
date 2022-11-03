#' Plot pedigrees with genotypes
#'
#' This is the main function for pedigree plotting, with many options for
#' controlling the appearance of pedigree symbols and accompanying labels. The
#' main pedigree layout is calculated with the `kinship2` package, see
#' [kinship2::align.pedigree] for details. Unlike `kinship2`, the implementation
#' here also supports plotting singletons. In addition, some minor adjustments
#' have been made to improve scaling and avoid unneeded duplications.
#'
#' @param x A [ped()] object.
#' @param marker Either a vector of names or indices referring to markers
#'   attached to `x`, a `marker` object, or a list of such. The genotypes for
#'   the chosen markers are written below each individual in the pedigree, in
#'   the format determined by `sep` and `missing`. See also `showEmpty`. If NULL
#'   (the default), no genotypes are plotted.
#' @param sep A character of length 1 separating alleles for diploid markers.
#' @param missing The symbol (integer or character) for missing alleles.
#' @param showEmpty A logical, indicating if empty genotypes should be included.
#' @param labs A vector or function controlling the individual labels included
#'   in the plot. Alternative forms:
#'
#'   * If `labs` is a vector with nonempty intersection with `labels(x)`, these
#'   individuals will be labelled. If the vector is named, then the (non-empty)
#'   names are used instead of the ID label. (See Examples.)
#'
#'   * If `labs` is NULL, or has empty intersection with `labels(x)`, then no
#'   labels are drawn.
#'
#'   * If `labs` is the word "num", then all individuals are numerically
#'   labelled following the internal ordering.
#'
#'   * If `labs` is a function, it is replaced with `labs(x)` and handled as
#'   above. (See Examples.)
#'
#' @param title The plot title. If NULL (default) or '', no title is added to
#'   the plot.
#' @param col A vector of colours for the pedigree members, recycled if
#'   necessary. Alternatively, `col` can be a list assigning colours to specific
#'   members. For example if `col = list(red = "a", blue = c("b", "c"))` then
#'   individual "a" will be red, "b" and "c" blue, and everyone else black. By
#'   default everyone is black.
#' @param aff A vector of labels identifying members whose plot symbols should
#'   be filled. (This is typically used in medical pedigrees to indicate
#'   affected members.)
#' @param carrier A vector of labels identifying members whose plot symbols
#'   should be marked with a dot. (This is typically used in medical pedigrees
#'   to indicate unaffected carriers of the disease allele.)
#' @param hatched A vector of labels identifying members whose plot symbols
#'   should be hatched.
#' @param deceased A vector of labels indicating deceased pedigree members.
#' @param starred A vector of labels indicating pedigree members that should be
#'   marked with a star in the pedigree plot.
#' @param twins A data frame with columns `id1`, `id2` and `code`, passed on to
#'   the `relation` parameter of [kinship2::plot.pedigree()].
#' @param hints A list with alignment hints passed on to
#'   [kinship2::align.pedigree()]. Rarely necessary.
#' @param fouInb Either "autosomal" (default), "x" or NULL. If "autosomal" or
#'   "x", inbreeding coefficients are added to the plot above the inbred
#'   founders. If NULL, or if no founders are inbred, nothing is added.
#' @param textInside,textAbove Character vectors of text to be printed inside or
#'   above pedigree symbols.
#' @param arrows A logical (default = FALSE). If TRUE, the pedigree is plotted
#'   as a DAG, i.e., with an arrow connecting each parent-child pair.
#' @param margins A numeric of length 4 indicating the plot margins. For
#'   singletons only the first element (the 'bottom' margin) is used.
#' @param keep.par A logical (default = FALSE). If TRUE, the graphical
#'   parameters are not reset after plotting, which may be useful for adding
#'   additional annotation.
#' @param \dots Arguments passed to internal methods. In particular `symbolsize`
#'   and `cex` can be useful.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [kinship2::plot.pedigree()]
#'
#' @examples
#'
#' x = nuclearPed(father = "fa", mother = "mo", child = "boy") |>
#'   addMarker(fa = "1/1", boy = "1/2", name = "SNP")
#'
#' plot(x, marker = 1)
#'
#' # or call marker by name
#' plot(x, marker = "SNP")
#'
#' # Modify margins
#' plot(x, margins = 6)
#' plot(x, margins = c(0,0,6,6)) # b,l,t,r
#'
#' # Other options
#' plot(x, marker = 1, hatched = typedMembers(x),
#'      starred = "fa", deceased = "mo")
#'
#' # Medical pedigree
#' plot(x, aff = males(x), carrier = "mo")
#'
#' # Label only some members
#' plot(x, labs = c("fa", "boy"))
#'
#' # Label only some members; rename the father
#' plot(x, labs = c(FATHER = "fa", "boy"))
#'
#' # Label males only
#' plot(x, labs = males)
#'
#' # Colours
#' plot(x, col = list(red = "fa", green = "boy"), hatched = "boy")
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
#' @export
plot.ped = function(x, marker = NULL, sep = "/", missing = "-", showEmpty = FALSE,
                    labs = labels(x), title = NULL, col = 1, aff = NULL, carrier = NULL,
                    hatched = NULL, deceased = NULL, starred = NULL, twins = NULL,
                    textInside = NULL, textAbove = NULL, arrows = FALSE,
                    hints = NULL, fouInb = "autosomal",
                    margins = 1, keep.par = FALSE, ...) {

  if(hasSelfing(x) && !arrows) {
    cat("Pedigree has selfing, switching to DAG mode. Use `arrows = TRUE` to avoid this message\n")
    arrows = TRUE
  }

  # Prepare annotations
  annot = .prepAnnotation(x, marker = marker, sep = sep, missing = missing,
                          showEmpty = showEmpty, labs = labs, starred = starred,
                          textInside = textInside, textAbove = textAbove, fouInb = fouInb,
                          col = col, aff = aff, hatched = hatched, carrier = carrier,
                          deceased = deceased)

  # Twin data: enforce data frame
  if(is.vector(twins))
    twins = data.frame(id1 = twins[1], id2 = twins[2], code = as.integer(twins[3]))

  # Pedigree alignment calculations
  plist = .alignPed(x, dag = arrows, hints = hints, twins = twins, ...)

  # Calculate scaling parameters and prepare plot window
  pdat = plotSetup(plist, textUnder = annot$textUnder, textAbove = annot$textAbove,
                   hasTitle = !is.null(title), mar = margins, ...)

  if(!keep.par)
    on.exit(par(pdat$oldpar))

  # Main plot
  .plotPed(pdat, dag = arrows, aff = annot$aff01, density = annot$density,
           angle = annot$angle, col = annot$colvec, ...)

  # Annotate
  .annotatePed(pdat, title = title, deceased = annot$deceasedTF, carrier = annot$carrierTF,
                textUnder = annot$textUnder, textInside = annot$textInside,
                textAbove = annot$textAbove, col = annot$colvec, ...)

  invisible(pdat)
}


#' @rdname plot.ped
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
#' # Modify the relative widths (which are not guessed)
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
#' plotPedList(peds, widths = w, groups = list(1, 2:3), titles = 1:2)
#'
#' # Parameters added in the main call are used in each sub-plot
#' plotPedList(peds, widths = w, margins = c(6, 3, 6, 3), labs = leaves,
#'             hatched = leaves, symbolsize = 1.3, col = list(red = 1))
#'
#' dev.off()
#'
#' #################################
#' # Example of automatic grouping #
#' #################################
#' H1 = nuclearPed()
#' H2 = list(singleton(1), singleton(3))  # grouped!
#'
#' plotPedList(list(H1, H2), dev.height = 2, dev.width = 4,
#'             titles = c(expression(H[1]), expression(H[2])))
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
#' x1 = nuclearPed(nch = 3)
#' m1 = marker(x1, `3` = "1/2")
#' marg1 = c(7, 4, 7, 4)
#' plot1 = list(x1, marker = m1, margins = marg1, title = "Plot 1",
#'              deceased = 1:2, cex = 1.3)
#'
#' x2 = cousinPed(2)
#' m2 = marker(x2, `11` = "A/A", `12` = "A/A")
#' marg2 = c(3, 4, 2, 4)
#' plot2 = list(x2, marker = m2, margins = marg2, title = "Plot 2",
#'              symbolsize = 1.2, labs = NULL)
#'
#' x3 = singleton("Mr. X")
#' plot3 = list(x3, title = "Plot 3", cex = 2, carrier = "Mr. X")
#'
#' x4 = halfSibPed()
#' hatched = 4:5
#' col = list(red = founders(x4), blue = leaves(x4))
#' marg4 = marg1
#' plot4 = list(x4, margins = marg4, title = "Plot 4", cex = 1.3,
#'              hatched = hatched, col = col)
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
#'             titles = c("Large", "Very large"),
#'             dev.height = 8, dev.width = 5, margins = 1.5)
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

  # If explicit source given, transfer marker data to all
  if(!is.null(source)) {
    srcPed = plots[[source]]
    if(is.null(srcPed))
      stop2("Unknown source pedigree: ", source)
    if(nMarkers(srcPed) == 0)
      stop2("The source pedigree has no attached markers")
    plots = lapply(plots, transferMarkers, from = srcPed)
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
  nvec = lapply(flatlist, function(arglist) .getPlist(arglist[[1]], ...)$n)

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
    dev.width = dev.width %||% sum(vapply(nvec, max, 1)) + 1
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
    mtext(finalTitles, outer = TRUE, at = midpoints, cex = cex.title)
  }
}


#' @rdname plot.ped
#' @export
#' @importFrom graphics plot.default
plot.list = function(x, marker = NULL, sep = "/", missing = "-", showEmpty = FALSE,
                     labs = unlist(labels(x), use.names = FALSE), col = 1, aff = NULL,
                     carrier = NULL, hatched = NULL, deceased = NULL, starred = NULL,
                     textInside = NULL, textAbove = NULL, arrows = FALSE, fouInb = "autosomal",
                     title = NULL, margins = 1, keep.par = FALSE, ...) {

  if(!is.pedList(x))
    return(plot.default(x, ...))

  if(hasSelfing(x) && !arrows) {
    cat("Pedigree has selfing, switching to DAG mode. Use `arrows = TRUE` to avoid this message\n")
    arrows = TRUE
  }

  x1 = x[[1]]
  x2 = x[[2]]

  # Compute and merge plists
  p1 = .getPlist(x1, dag = arrows)
  p2 = .getPlist(x2, dag = arrows)

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
  plist2 = .extendPlist(plist, x)

  # Adjust y positions
  #idx = match(6:8, plist2$plotord)
  #plist2$yall[idx] = plist2$yall[idx] + 0.5

  # Merge annotations
  annot1 = .prepAnnotation(x1, marker = marker, sep = sep, missing = missing,
                           showEmpty = showEmpty, labs = labs, starred = starred,
                           textInside = textInside, textAbove = textAbove, fouInb = fouInb,
                           col = col, aff = aff, hatched = hatched, carrier = carrier,
                           deceased = deceased)
  annot2 = .prepAnnotation(x2, marker = marker, sep = sep, missing = missing,
                           showEmpty = showEmpty, labs = labs, starred = starred,
                           textInside = textInside, textAbove = textAbove, fouInb = fouInb,
                           col = col, aff = aff, hatched = hatched, carrier = carrier,
                           deceased = deceased)

  mrg = function(a, def)
    c(annot1[[a]] %||% rep(def, nInd1), annot2[[a]] %||% rep(def, nInd2))

  annot = list(textUnder = mrg("textUnder", ""),
               textAbove = mrg("textAbove", ""),
               textInside = mrg("textInside", ""),
               colvec = mrg("colvec", 1),
               aff01 = mrg("aff01", 0),
               density = annot1$density,
               angle = annot1$angle,
               carrierTF = mrg("carrierTF", FALSE),
               deceasedTF = mrg("deceasedTF", FALSE))

  # Calculate scaling parameters and prepare plot window
  pdat = plotSetup(plist2, textUnder = annot$textUnder, textAbove = annot$textAbove,
                   hasTitle = !is.null(title), mar = margins, ...)

  if(!keep.par)
    on.exit(par(pdat$oldpar))

  # Main plot
  .plotPed(pdat, dag = arrows, aff = annot$aff01, density = annot$density,
           angle = annot$angle, col = annot$colvec, ...)

  # Annotate
  .annotatePed(pdat, title = title, deceased = annot$deceasedTF, carrier = annot$carrierTF,
               textUnder = annot$textUnder, textInside = annot$textInside,
               textAbove = annot$textAbove, col = annot$colvec, ...)

  invisible(pdat)
}




# Function for preparing various plot annotations
.prepAnnotation = function(x, marker = NULL, sep = "/", missing = "-", showEmpty = FALSE,
                           labs = labels(x), col = 1, aff = NULL, carrier = NULL,
                           hatched = NULL, deceased = NULL, starred = NULL,
                           textInside = NULL, textAbove = NULL, fouInb = "autosomal") {

  res = list()
  nInd = pedsize(x)

  # Labels ------------------------------------------------------------------

  if(is.function(labs))
    labs = labs(x)

  if(identical(labs, "num"))
    labs = setNames(x$ID, 1:nInd)

  text = .prepLabs(x, labs)

  # Add stars to labels
  if(is.function(starred))
    starred = starred(x)
  starred = internalID(x, starred, errorIfUnknown = FALSE)
  starred = starred[!is.na(starred)]
  text[starred] = paste0(text[starred], "*")

  # Marker genotypes
  if (length(marker) > 0) { # excludes NULL and empty vectors/lists
    if (is.marker(marker))
      mlist = list(marker)
    else if (is.markerList(marker))
      mlist = marker
    else if (is.numeric(marker) || is.character(marker) || is.logical(marker))
      mlist = getMarkers(x, markers = marker)
    else
      stop2("Argument `marker` must be either:\n",
            "  * a n object of class `marker`\n",
            "  * a list of `marker` objects\n",
            "  * a character vector (names of attached markers)\n",
            "  * an integer vector (indices of attached markers)",
            "  * a logical vector of length `nMarkers(x)`")
    checkConsistency(x, mlist)

    gg = do.call(cbind, lapply(mlist, format, sep = sep, missing = missing))
    geno = apply(gg, 1, paste, collapse = "\n")
    if (!showEmpty)
      geno[rowSums(do.call(cbind, mlist)) == 0] = ""

    text = if (!any(nzchar(text))) geno else paste(text, geno, sep = "\n")
  }

  res$textUnder = trimws(text, "right", whitespace = "[\t\r\n]")


  # Text above symbols ------------------------------------------------------

  showFouInb = !is.null(fouInb) && hasInbredFounders(x)

  if(is.function(textAbove))
    textAbove = textAbove(x)
  else if(showFouInb) {
    finb = founderInbreeding(x, chromType = fouInb, named = TRUE)
    finb = finb[finb > 0]
    textAbove = sprintf("f = %.4g", finb)
    names(textAbove) = names(finb)
  }
  if(!is.null(textAbove))
    mode(textAbove) = "character"

  nmsAbove = names(textAbove)
  if(!is.null(nmsAbove)) {
    tmp = character(nInd)
    tmp[internalID(x, nmsAbove, errorIfUnknown = FALSE)] = textAbove
    textAbove = tmp
  }

  res$textAbove = textAbove

  # Text inside symbols ------------------------------------------------------

  if(is.function(textInside))
    textInside = textInside(x)
  if(!is.null(textInside))
    mode(textInside) = "character"

  nmsInside = names(textInside)
  if(!is.null(nmsInside)) {
    tmp = character(nInd)
    tmp[internalID(x, nmsInside, errorIfUnknown = FALSE)] = textInside
    textInside = tmp
  }

  res$textInside = textInside

  # Colours -----------------------------------------------------------------

  if(!is.list(col))
    colvec = rep(col, length = nInd)
  else { # E.g. list(red = 1:2, blue = 3:4)
    colvec = rep(1, nInd)
    for(cc in names(col)) {
      thiscol = col[[cc]]
      if(is.function(thiscol))
        idscol = thiscol(x)
      else
        idscol = intersect(x$ID, thiscol)
      colvec[internalID(x, idscol)] = cc
    }
  }

  res$colvec = colvec

  # Affected/hathced --------------------------------------------------------

  if(is.function(aff))
    aff = aff(x)
  if(is.function(hatched))
    hatched = hatched(x)
  if(!is.null(aff) && !is.null(hatched))
    stop2("Both `aff` and `hatched` cannot both be used")

  # filling density and angle
  density = if(!is.null(aff)) -1 else if(!is.null(hatched)) 25 else NULL
  angle = if(!is.null(aff)) 90 else if(!is.null(hatched)) 45 else NULL

  if(!is.null(hatched))
    aff = hatched # for kinship2

  res$aff01 = ifelse(x$ID %in% aff, 1, 0)
  res$density = density
  res$angle = angle

  # Carriers ----------------------------------------------------------------

  if(is.function(carrier))
    carrier = carrier(x)

  # Convert to T/F
  res$carrierTF = x$ID %in% carrier

  # Deceased ----------------------------------------------------------------

  if(is.function(deceased))
    deceased = deceased(x)

  # Convert to T/F
  res$deceasedTF = x$ID %in% deceased

  # Return list -------------------------------------------------------------
  res
}

# Convert `labs` to full vector with "" for unlabels indivs
.prepLabs = function(x, labs) {
  id = rep("", length(x$ID)) # Initialise

  mtch = match(x$ID, labs, nomatch = 0L)
  showIdx = mtch > 0
  showLabs = labs[mtch]

  if(!is.null(nms <- names(labs))) { # use names(labs) if present
    newnames = nms[mtch]
    goodIdx = newnames != "" & !is.na(newnames)
    showLabs[goodIdx] = newnames[goodIdx]
  }

  id[showIdx] = showLabs
  id
}
