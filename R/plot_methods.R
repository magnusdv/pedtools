# Modified from `kinship2::plot.pedigree()`. The main changes are:
#
# * Separate setup (dimensions and scaling), drawing and annotation (see below)
# * Fixed scaling bugs, documented here: https://github.com/mayoverse/kinship2/pull/10
# * Adjust scaling to account for duplication arrows
# * Avoid unwanted duplications, e.g. in 3/4 siblings
# * Don't use bottom labels to calculate inter-generation separation
# * Allow plotting pedigrees as DAGs
#
# * Removed features not used by pedtools, including i) multiple phenotypes and ii) subregions
# * Fixated some parameters at their default values:
#     - pconnect = 0.5
#     - branch = 0.6
#     - packed = TRUE
#
# Internally, the previous single `plot()` function has been refactored into the
# following steps:
# * .pedAlignment(): Builds on `kinship2::align.pedigree()`, but also handles singletons and DAGs
# * .pedAnnotation(): Prepare and collect annotation variables
# * .pedScaling(): Calculate symbol sizes and scaling variables
# * .drawPed(): Draw symbols
# * .annotatePed(): Add labels and other annotation




#' Internal plot methods
#'
#' The main purpose of this page is to document the many options for pedigree
#' plotting. Most of the arguments shown here may be supplied directly in
#' `plot(x, ...)`, where `x` is a pedigree. See [plot.ped()] for many examples.
#'
#' The workflow of `plot.ped(x, ...)` is approximately as follows:
#'
#' ```
#' # Calculate plot parameters
#' align = .pedAlignment(x, ...)
#' annot = .pedAnnotation(x, ...)
#' scale = .pedScaling(align, annot, ...)
#'
#' # Produce plot
#' .drawPed(align, annot, scale)
#' .annotatePed(align, annot, scale)  # if `annot` contains text annotation etc
#' ```
#'
#' The `labs` argument controls the individual ID labels printed below the
#' pedigree symbols. By default the output of `labels(x)` is used, but there are
#' several alternative forms:
#'
#'   * If `labs` is a vector with nonempty intersection with `labels(x)`, only
#' these individuals will be labelled. If the vector is named, then the names
#' are used instead of the ID label. (See Examples.)
#'
#'   * If `labs` is the word "num", then all individuals are numerically
#' labelled following the internal ordering.
#'
#'   * Use `labs = NULL` to remove all labels.
#'
#'   * If `labs` is a function, it is replaced with `labs(x)` and handled as
#' above. (See Examples.)
#'
#' The argument `textAnnot` allows customised annotation around and inside each
#' symbol. This takes a list of lists, whose names may include "topleft",
#' "topright", "left", "right", "bottomleft", "bottom", "bottomright" and
#' "inside". Each inner list should contain a character vector as its first
#' element (with the text to be printed), followed by further arguments passed
#' to [text()]. For example, `textAnnot = list(left = list(c(A = "1"), cex =
#' 2))` prints a large number "1" to the left of individual A (if such an
#' individual exists in the pedigree). See Examples.
#'
#' The arguments `col`, `fill`, `lty` and `lwd` can all be indicated in a number
#' of ways:
#'
#'   * An unnamed vector. This will be recycled and applied to all members. For
#' example, `lty = 2` gives everyone a dashed outline.
#'
#'   * A named vector. Only pedigree members appearing in the names are affected.
#' Example: `fill = c("1" = "red", foo = "blue")` fills individual `1` red and
#' `foo` blue.
#'
#'   * A list of ID vectors, where the list names indicate the parameter values.
#' Example: `col = list(red = 1:2, blue = 3:5)`.
#'
#'   * List entries may also be functions, taking the pedigree `x` as input and
#' producing a vector of ID labels. The many built-in functions in
#' [ped_subgroups] are particularly handy here, e.g.: `fill = list(red =
#' founders, blue = leaves)`.
#'
#' @param x A [ped()] object.
#' @param plist Alignment list with format similar to
#'   [kinship2::align.pedigree()].
#' @param arrows A logical (default = FALSE). If TRUE, the pedigree is plotted
#'   as a DAG, i.e., with arrows connecting parent-child pairs.
#' @param labs A vector or function controlling the individual labels in the
#'   plot. By default, `labels(x)` are used. See Details for valid formats.
#' @param foldLabs A number or function controlling the folding of long labels.
#'   If a number, line breaks are inserted at roughly this width, trying to
#'   break at break-friendly characters. If a function, this is applied to each
#'   label.
#' @param trimLabs A logical, by default TRUE. Removes line breaks and tabs from
#'   both ends of the labels (after adding genotypes, if `marker` is not NULL).
#' @param cex Expansion factor controlling font size. This also affects symbol
#'   sizes, which by default have the width of 2.5 characters. Default: 1.
#' @param symbolsize Expansion factor for pedigree symbols. Default: 1.
#' @param margins A numeric indicating the plot margins. If a single number is
#'   given, it is recycled to length 4.
#' @param addSpace A numeric of length 4, indicating extra padding (in inches)
#'   around the pedigree inside the plot region. Default: 0.
#' @param xlim,ylim Numeric vectors of length 2, used to set `par("usr")`
#'   explicitly. Rarely needed by end users.
#' @param vsep2 A logical; for internal use.
#' @param autoScale A logical. It TRUE, an attempt is made to adjust `cex` so
#'   that the symbol dimensions are at least `minsize` inches. Default: FALSE.
#' @param minsize A positive number, by default 0.15. (See `autoScale`.)
#' @param debug A logical, turning on messages from the autoscale algorithm.
#' @param marker Either a vector of names or indices referring to markers
#'   attached to `x`, a `marker` object, or a list of such. The genotypes for
#'   the chosen markers are written below each individual in the pedigree, in
#'   the format determined by `sep` and `missing`. See also `showEmpty`. If NULL
#'   (the default), no genotypes are plotted.
#' @param sep A character of length 1 separating alleles for diploid markers.
#' @param missing The symbol (integer or character) for missing alleles.
#' @param showEmpty A logical, indicating if empty genotypes should be included.
#' @param title The plot title. If NULL (default) or '', no title is added to
#'   the plot.
#' @param col A vector or list specifying outline colours for the pedigree
#'   members. See Details for valid formats.
#' @param fill A vector or list specifying fill/hatch colours for the pedigree
#'   members. See Details for valid formats. Note that if `fill` is unnamed, and
#'   either `aff` or `hatched` are given, then the fill colour is applied only
#'   to those.
#' @param lty,lwd Vectors or lists specifying linetype and width of pedigree
#'   symbol outlines. See Details for valid formats.
#' @param hatched A vector of labels identifying members whose plot symbols
#'   should be hatched.
#' @param hatchDensity A number specifying the hatch density in lines per inch.
#'   Default: 25.
#' @param aff A vector of labels identifying members whose plot symbols should
#'   be filled. (This is typically used in medical pedigrees to indicate
#'   affected members.)
#' @param carrier A vector of labels identifying members whose plot symbols
#'   should be marked with a dot. (This is typically used in medical pedigrees
#'   to indicate unaffected carriers of the disease allele.)
#' @param deceased A vector of labels indicating deceased pedigree members.
#' @param proband A vector of labels indicating proband individuals, to be
#'   marked with an arrow in the plot.
#' @param starred A vector of labels indicating pedigree members that should be
#'   marked with a star in the pedigree plot.
#' @param twins A data frame with columns `id1`, `id2` and `code`, passed on to
#'   the `relation` parameter of [kinship2::plot.pedigree()].
#' @param miscarriage A vector of labels indicating miscarriages, shown as
#'   triangles in the pedigree plot.
#' @param straight A logical, indicating if the plot should (attempt to) use
#'   straight lines everywhere. Default: FALSE.
#' @param packed,width,align Parameters passed on to
#'   [kinship2::align.pedigree()]. Can usually be left untouched.
#' @param spouseOrder An optional vector (or list of vectors) indicating plot
#'   ordering for spouses. (This is converted into a matrix and forward as
#'   `hints`; see below.)
#' @param hints An optional list of hints passed on to
#'   [kinship2::align.pedigree()].
#' @param fouInb Either "autosomal" (default), "x" or NULL. If "autosomal" or
#'   "x", inbreeding coefficients are added to the plot above the inbred
#'   founders. If NULL, or if no founders are inbred, nothing is added.
#' @param textInside,textAbove Character vectors of text to be printed inside or
#'   above pedigree symbols. \[Soft deprecated; replaced by `textAnnot`.\]
#' @param textAnnot A list specifying further text annotation around or inside
#'   the pedigree symbols. See Details for more information.
#' @param font,fam Arguments passed on to [text()].
#' @param colUnder,colInside,colAbove Colour vectors.
#' @param cex.main,line.main,col.main,font.main Parameters passed on to
#'   [title()].
#' @param alignment List of alignment details, as returned by [.pedAlignment()].
#' @param annotation List of annotation details as returned by
#'   [.pedAnnotation()].
#' @param scaling List of scaling parameters as returned by [.pedScaling()].
#' @param \dots Further parameters passed between methods.
#'
#' @examples
#' x = nuclearPed()
#'
#' align = .pedAlignment(x)
#' annot = .pedAnnotation(x)
#' scale = .pedScaling(align, annot)
#'
#' drawPed(align, annot, scale)
#'
#' @name plotmethods
NULL




# Alignment ---------------------------------------------------------------

#' @rdname plotmethods
#' @importFrom kinship2 align.pedigree
#' @export
.pedAlignment = function(x = NULL, plist = NULL, arrows = FALSE, twins = NULL,
                         miscarriage = NULL, packed = TRUE, width = 10,
                         straight = FALSE, align = NULL,
                         spouseOrder = NULL, hints = NULL, ...) {

  if(hasSelfing(x) && !arrows) {
    message("Pedigree has selfing, switching to DAG mode. Use `arrows = TRUE` to avoid this message.")
    arrows = TRUE
  }

  if(!is.null(plist))
    return(.extendPlist(x, plist, arrows = arrows, miscarriage = miscarriage))

  # Singleton
  if(is.singleton(x)) {
    plist = list(n = 1, nid = cbind(1), pos = cbind(0), fam = cbind(0), spouse = cbind(0))
    return(.extendPlist(x, plist, miscarriage = miscarriage))
  }

  if(arrows)
    return(.alignDAG(x))

  # Twin data: enforce data frame
  if(is.vector(twins))
    twins = data.frame(id1 = twins[1], id2 = twins[2], code = as.integer(twins[3]))

  # (Try to) force spouse order
  if(!is.null(spouseOrder)) {
    if(!is.null(hints))
      stop2("Cannot use both `hints` and `spouseOrder` in the same call")
    hints = .spouseOrder(x, spouseOrder)
  }

  k2ped = as_kinship2_pedigree(x, twins = twins)
  align = align %||% if(straight) c(0,0) else c(1.5, 2)
  plist = kinship2::align.pedigree(k2ped, packed = packed, width = width, align = align, hints = hints)

  # Catch missing persons (kindepth bug!)
  ERR = sum(plist$n) < length(x$ID)
  if(ERR) {
    warning("Alignment failed; switching to simple DAG alignment", call. = FALSE)
    return(.alignDAG(x))
  }

  # Ad hoc fix for 3/4 siblings and similar
  if(is.null(hints))
    plist = .fix34(x, k2ped = k2ped, plist = plist, packed = packed, width = width, align = align)

  # Fix annoying rounding errors in first column of `pos`
  plist$pos[] = round(plist$pos[], 6)

  # Add further parameters
  .extendPlist(x, plist, miscarriage = miscarriage)
}


.extendPlist = function(x, plist, arrows = FALSE, miscarriage = NULL) {
  nid = plist$nid
  pos = plist$pos

  nInd = max(nid)
  maxlev = nrow(pos)

  id = as.vector(nid)
  plotord = id[id > 0]

  # Coordinates (top center)
  xall = pos[id > 0]
  yall = row(pos)[id > 0]

  xrange = range(xall)
  yrange = range(yall)

  # For completeness: Kinship2 order (1st instance of each only!)
  # Including this for completeness
  tmp = match(1:nInd, nid)
  xpos = pos[tmp]
  ypos = row(pos)[tmp]

  sex = getSex(x)
  if(length(miscarriage)) {
    idx = match(miscarriage, x$ID, nomatch = 0L)
    idx = idx[idx > 0]
    if(any(idx %in% c(x$FID, x$MID)))
      stop2("A parent cannot assigned as a miscarriage: ", miscarriage[idx %in% c(x$FID, x$MID)])
    sex[idx] = 3
  }

  list(plist = plist, x = xpos, y = ypos, nInd = nInd, sex = sex, ped = x, arrows = arrows,
       plotord = plotord, xall = xall, yall = yall, maxlev = maxlev, xrange = xrange, yrange = yrange)
}


# Annotation --------------------------------------------------------------

#' @rdname plotmethods
#' @export
.pedAnnotation = function(x, title = NULL, marker = NULL, sep = "/", missing = "-", showEmpty = FALSE,
                          labs = labels(x), foldLabs = 12, trimLabs = TRUE, col = 1, fill = NA, lty = 1, lwd = 1,
                          hatched = NULL, hatchDensity = 25, aff = NULL, carrier = NULL,
                          deceased = NULL, starred = NULL, proband = NULL, textAnnot = NULL,
                          textInside = NULL, textAbove = NULL, fouInb = "autosomal", ...) {

  res = list()
  nInd = pedsize(x)


  # Title -------------------------------------------------------------------
  if(is.function(title))
    title = title(x)
  res$title = title


  # Labels ------------------------------------------------------------------

  if(is.function(labs))
    labs = labs(x)

  if(identical(labs, "num"))
    labs = setNames(x$ID, 1:nInd)

  textu = .prepLabs(x, labs)

  # Fold
  if(isCount(foldLabs))
    textu = vapply(textu, function(s) smartfold(s, width = foldLabs), FUN.VALUE = "")
  else if(is.function(foldLabs))
    textu = vapply(textu, foldLabs, FUN.VALUE = "")

  # Add stars to labels
  if(is.function(starred))
    starred = starred(x)
  starred = internalID(x, starred, errorIfUnknown = FALSE)
  starred = starred[!is.na(starred)]
  textu[starred] = paste0(textu[starred], "*")

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

    if(is.logical(showEmpty) && length(showEmpty) == 1)
      showEmpty = if(showEmpty) x$ID else NULL
    else if (is.function(showEmpty))
      showEmpty = showEmpty(x)

    hideEmpty = match(x$ID, showEmpty, nomatch = 0L) == 0
    if (any(hideEmpty)) {
      isEmpty = rowSums(do.call(cbind, mlist)) == 0
      geno[isEmpty & hideEmpty] = ""
    }

    textu = if (!any(nzchar(textu))) geno else paste(textu, geno, sep = "\n")
  }

  if(trimLabs)
    textu = trimws(textu, which = "both", whitespace = "[\t\r\n]")

  res$textUnder = textu

  # Further text annotation

  if(!is.null(textAnnot)) {
    res$textAnnot = lapply(textAnnot, function(b) {
      if(is.atomic(b))
        b = list(b)
      b[[1]] = .prepLabs2(x, b[[1]])
      b
    })
  }


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

  res$textAbove = .prepLabs2(x, textAbove)

  # Text inside symbols ------------------------------------------------------

  if(is.function(textInside))
    textInside = textInside(x)

  res$textInside = .prepLabs2(x, textInside)


  # Affected/hathced --------------------------------------------------------

  if(is.function(aff))
    aff = aff(x)
  if(is.function(hatched))
    hatched = hatched(x)
  isaff = x$ID %in% aff
  ishatch = x$ID %in% hatched

  # filling density (-1 = fill; 25 = hatch)
  densvec = integer(nInd)
  densvec[isaff] = -1
  densvec[ishatch] = hatchDensity
  res$densvec = densvec

  # See fill color below!

  # Colours (border)----------------------------------------------------------

  res$colvec = .prepPlotarg(x, col, default = 1)

  # Fill color --------------------------------------------------------------

  affORhatch = isaff | ishatch

  # If aff/hatch given apply simple fill only to those
  if(any(affORhatch) && !is.list(fill) && is.null(names(fill)) && !identical(fill, NA)) {
    fillvec = rep(NA, length = nInd)
    fillvec[affORhatch] = fill
  }
  else
    fillvec = .prepPlotarg(x, fill, default = NA)

  # Ensure aff/hatch are filled
  fillvec[affORhatch & is.na(fillvec)] = 1

  res$fillvec = fillvec


  # Linetype ----------------------------------------------------------------
  ltyvec = .prepPlotarg(x, lty, default = 1)

  if(any(badlty <- !ltyvec %in% 0:6)) {
    ltynames = c("blank", "solid", "dashed", "dotted", "dotdash", "longdash", "twodash")
    ltyvec[badlty] = match(ltyvec[badlty], ltynames, nomatch = 2) - 1
  }
  res$ltyvec = as.numeric(ltyvec)

  # Line width ----------------------------------------------------------------

  res$lwdvec = as.numeric(.prepPlotarg(x, lwd, default = 1))

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


  # Proband -----------------------------------------------------------------

  if(is.function(proband))
    proband = proband(x)

  # Convert to T/F
  res$probandTF = x$ID %in% proband

  # Return list -------------------------------------------------------------
  res
}


# Convert `labs` to full vector with "" for unlabels indivs
.prepLabs = function(x, labs) {
  id = rep("", length(x$ID)) # Initialise

  mtch = match(x$ID, labs, nomatch = 0L)
  showIdx = mtch > 0
  showLabs = labs[mtch]

  # Use names(labs) if present
  if(!is.null(nms <- names(labs))) {
    newnames = nms[mtch]
    goodIdx = newnames != "" & !is.na(newnames)
    showLabs[goodIdx] = newnames[goodIdx]
  }

  id[showIdx] = showLabs
  id
}

# Alternative to .prepLabs (used above/inside): If vector has names, match these to x$ID.
.prepLabs2 = function(x, labs) {
  mode(labs) = "character"

  if(is.null(names(labs)) && length(labs) == length(x$ID))
    return(labs)

  ids = names(labs) %||% labs
  mtch = match(x$ID, ids, nomatch = 0L)

  txt = rep("", length(x$ID)) # Initialise
  txt[mtch > 0] = labs[mtch]
  txt
}


#--- Plot dimension and scaling parameters

#' @rdname plotmethods
#' @importFrom graphics frame strheight strwidth
#' @export
.pedScaling = function(alignment, annotation, cex = 1, symbolsize = 1, margins = 1,
                       addSpace = 0, xlim = NULL, ylim = NULL, vsep2 = FALSE,
                       autoScale = FALSE, minsize = 0.15, debug = FALSE, ...) {

  textUnder = annotation$textUnder
  textAbove = annotation$textAbove
  title = annotation$title

  maxlev = alignment$maxlev
  xrange = alignment$xrange
  yrange = alignment$yrange
  nid = alignment$plist$nid
  nid1 = nid[1, ][nid[1, ] > 0] # ids in first generation

  # Fix xrange/yrange for singletons and selfings
  if(maxlev == 1 || xrange[1] == xrange[2])
    xrange = xrange + c(-0.5, 0.5)

  if(maxlev == 1)
    yrange = yrange + c(-0.5, 0.5)

  # Margins
  mar = margins
  if(length(mar) == 1)
    mar = if(!is.null(title)) c(mar, mar, mar + 2.1, mar) else rep(mar, 4)

  # Adjust margin for proband arrows?
  if(any(prb <- annotation$probandTF)) {
    idx = prb[alignment$plotord]
    # Hard coded: All arrows are bottom left
    arrowL = any(alignment$xall[idx] == xrange[1])
    arrowB = any(alignment$yall[idx] == yrange[2])

    extraMar = c(if(arrowB) 0.5 else 0, if(arrowL) 2.5 else 0, 0, 0)
    mar = mar + extraMar
  }

  # Extra padding (e.g. for ribd::ibdDraw() and ibdsim2::haploDraw())
  if(length(addSpace) == 1)
    addSpace = rep(addSpace, 4)

  # Set margin and xpd
  oldpar = par(mar = mar, xpd = TRUE)

  # Dimensions in inches
  psize = par('pin')

  # Shortcut for finding height of a string in inches. Empty -> 0!
  hinch = function(v) {
    res = strheight(v, units = "inches", cex = cex)
    res[v == ""] = 0
    res
  }

  # Text height
  labh_in = hinch('M') # same as for "1g" used in kinship2 (`stemp2`)

  # Make room for curved duplication lines involving first generation
  # A bit hackish since curve height is only available in user coordinates.
  curvAdj = if(anyDuplicated.default(nid1)) 0.5 else if(maxlev > 1 && any(nid1 %in% nid[2, ])) 0.1225 else 0

  # Text above symbols in first generation
  # Don't adjust for text above if also for curve. (NB: Fails in extreme cases)
  abovetop_in = if(curvAdj>0) 0 else max(hinch(textAbove[nid1]))

  # Add offset: "0.5 times the width [!] of a character"
  if(abovetop_in > 0)
    abovetop_in = abovetop_in + strwidth("W", units = "inches")/2

  abovetop_in = abovetop_in + addSpace[3]

  # Separation above/below labels
  labsep1_in = 0.7*labh_in  # above text
  labsep2_in = labh_in  # below text
  labsep3_in = 0.3*labh_in  # below text in last generation

  # Max label height (except last generation)
  maxlabh_nolast_in = if(maxlev == 1) 0 else max(hinch(textUnder[nid[seq_len(maxlev - 1), ]]))

  # Max label height in last generation
  belowlast_in = max(hinch(textUnder[nid[maxlev, ]]))

  # Everything below bottom symbol (label + space above and below)
  if(belowlast_in > 0)
    belowlast_in = belowlast_in + labsep1_in + labsep3_in

  belowlast_in = belowlast_in + addSpace[1]

  # KEY TO LAYOUT: Complete y-range (psize[2]) in inches corresponds to:
  # abovetop_in + (labsep1_in + h + maxlabh_nolast_in + labsep2_in) * (maxlev - 1) + (h + belowlast_in)

  # Symbol height restriction 1 (Solve above for h)
  if(maxlev > 1) {
    sep1_in = (labsep1_in + maxlabh_nolast_in + labsep2_in)
    ht1 = (psize[2] - abovetop_in - belowlast_in - sep1_in * (maxlev - 1)) / maxlev
  }
  else
    ht1 = psize[2] - abovetop_in - belowlast_in

  # Height restriction 2
  ht2 = psize[2]/(maxlev + (maxlev-1)/2)

  # Width restriction 1: Default width = 2.5 letters (`stemp1`)
  wd1 = strwidth("ABC", units='inches', cex=cex) * 2.5/3

  # Width restriction 2
  wd2 = psize[1] * 0.8/(.8 + diff(xrange))  # = psize[1] for singletons/selfings

  # Box size in inches
  boxsize = symbolsize * min(ht1, ht2, wd1, wd2)

  # Autoscale if too small
  if(autoScale && boxsize < minsize) {
    if(minsize > ht2 | cex < 0.2)
      stop2("autoScale error: `minsize` is too large for the current window")

    trycex = round(0.95 * cex, 2)
    trysymbolsize = round(1.05 * symbolsize, 2)
    if(debug)
      message(sprintf("autoScale: cex = %g, symbolsize = %g", trycex, trysymbolsize))

    return(.pedScaling(alignment, annotation, cex = trycex, symbolsize = trysymbolsize,
                       margins = margins, addSpace = addSpace, xlim = xlim, ylim = ylim,
                       vsep2 = vsep2, autoScale = TRUE, minsize = minsize, debug = debug))
  }

  if (ht1 <= 0)
    stop2("Labels leave no room for the graph, reduce cex")

  # Horizontal scaling factor
  if(is.null(xlim)) {
    # Segment corresponding to 1 unit
    hscale = (psize[1] - boxsize - addSpace[2] - addSpace[4])/diff(xrange)
  }
  else { # override if xlim provided!
    hscale = psize[1]/diff(xlim)
  }

  # Vertical scaling factor
  if(is.null(ylim)) {
    denom = if(maxlev == 1) 1 else maxlev - 1 + curvAdj
    vscale = (psize[2] - (abovetop_in + boxsize + belowlast_in)) / denom
  }
  else {
    vscale = psize[2]/diff(ylim)
  }

  if(hscale <= 0 || vscale <= 0)
    stop2("Cannot fit the graph; please increase plot region or reduce cex and/or symbolsize")

  boxw = boxsize/hscale  # box width in user units
  boxh = boxsize/vscale  # box height in user units

  # User coordinates
  if(is.null(xlim)) {
    left   = xrange[1] - 0.5*boxw - addSpace[2]/hscale
    right  = xrange[2] + 0.5*boxw + addSpace[4]/hscale
  }
  else {
    left = xlim[1]
    right = xlim[2]
  }
  if(is.null(ylim)) {
    top    = yrange[1] - abovetop_in/vscale - curvAdj
    bottom = yrange[2] + (boxsize + belowlast_in)/vscale
  }
  else {
    top = min(ylim)
    bottom = max(ylim)
  }
  usr = c(left, right, bottom, top)

  labh = labh_in/vscale        # height of a text string
  legh = min(1/4, boxh * 1.5)  # how tall are the 'legs' up from a child

  # Return plotting/scaling parameters
  list(boxw = boxw, boxh = boxh, labh = labh, legh = legh, vsep2 = vsep2,
       cex = cex, mar = mar, usr = usr, oldpar = oldpar)
}


#' @rdname plotmethods
#' @importFrom graphics lines polygon segments
#' @export
.drawPed = function(alignment, annotation, scaling) {

  if(isTRUE(alignment$arrows))
    return(.plotDAG(alignment, annotation, scaling))

  n = alignment$nInd
  plotord = alignment$plotord
  xall = alignment$xall
  yall = alignment$yall
  maxlev = alignment$maxlev
  plist = alignment$plist
  SEX = alignment$sex

  pos = plist$pos
  nid = plist$nid

  boxh = scaling$boxh
  boxw = scaling$boxw
  legh = scaling$legh
  vsep2 = scaling$vsep2

  branch = 0.6
  pconnect = .5

  COL = annotation$colvec %||% 1
  FILL = annotation$fillvec %||% NA
  LTY = annotation$ltyvec %||% 1
  DENS = annotation$densvec %||% 0
  LWD = annotation$lwdvec %||% 1

  if (length(COL) == 1)
    COL = rep(COL, n)
  if (length(FILL) == 1)
    FILL = rep(FILL, n)
  if (length(LTY) == 1)
    LTY = rep(LTY, n)
  if (length(DENS) == 1)
    DENS = rep(DENS, n)
  if (length(LWD) == 1)
    LWD = rep(LWD, n)

  # Set user coordinates
  par(mar = scaling$mar, usr = scaling$usr, xpd = TRUE)

  # Shapes
  POLYS = list(list(x = c(0, -0.5, 0, 0.5), y = c(0, 0.5, 1, 0.5)), # diamond
               list(x = c(-1, -1, 1, 1)/2, y = c(0, 1, 1, 0)),      # square
               list(x = 0.5 * cos(seq(0, 2 * pi, length = 50)),     # circle
                    y = 0.5 * sin(seq(0, 2 * pi, length = 50)) + 0.5),
               list(x = c(0, -1.2, 1.2)/2, y = c(0, 0.6, 0.6)))     # triangle

  # Draw all the symbols
  for (k in seq_along(plotord)) {
    id = plotord[k]
    poly = POLYS[[SEX[id] + 1]]
    dens = if(DENS[id] == 0) NULL else DENS[id]
    polygon(xall[k] + poly$x * boxw,
            yall[k] + poly$y * boxh,
            border = COL[id], col = FILL[id],
            lty = LTY[id], lwd = LWD[id],
            angle = 45, density = dens)
  }

  ## Add lines between spouses (MDV: Vectorized/simplified)
  sp = plist$spouse
  cl = col(sp)[sp > 0]
  rw = row(sp)[sp > 0]
  tmpy = rw + boxh/2
  segments(pos[cbind(rw, cl)]     + boxw/2, tmpy,
           pos[cbind(rw, cl + 1)] - boxw/2, tmpy)

  # Double line for consanguineous marriage
  if(any(sp == 2)) {
    cl2 = col(sp)[sp == 2]
    rw2 = row(sp)[sp == 2]
    tmpy2 = rw2 + boxh/2 + boxh/10
    segments(pos[cbind(rw2, cl2)]     + boxw/2, tmpy2,
             pos[cbind(rw2, cl2 + 1)] - boxw/2, tmpy2)
  }

  ## Lines from offspring to parents
  ## NB: If vsep2 = T, parents are two rows above (Hack used in plot.list().)
  chRows = if(vsep2 && maxlev > 1) seq_len(maxlev-2) + 2 else seq_len(maxlev-1) + 1
  for(i in chRows) {
    parentRow = if(vsep2) i-2 else i-1
    zed = unique.default(plist$fam[i,  ]) # MDV: use unique.default
    zed = zed[zed > 0]  #list of family ids

    for(fam in zed) {
      xx = pos[parentRow, fam + 0:1]
      parentx = mean(xx)   #midpoint of parents

      # Draw the uplines
      who = (plist$fam[i,] == fam) #The kids of interest
      if (is.null(plist$twins))
        target = pos[i,who]
      else {
        twin.to.left = (c(0, plist$twins[i,who])[1:sum(who)])
        temp = cumsum(twin.to.left == 0) #increment if no twin to the left
        # 5 sibs, middle 3 are triplets gives 1,2,2,2,3
        # twin, twin, singleton gives 1,1,2,2,3
        tcount = table(temp)
        target = rep(tapply(pos[i,who], temp, mean), tcount)
      }

      yy = rep(i, sum(who))
      segments(pos[i,who], yy, target, yy-legh)

      ## Draw midpoint MZ twin line
      if (any(plist$twins[i,who] == 1)) {
        who2 = which(plist$twins[i,who] == 1)
        temp1 = (pos[i, who][who2] + target[who2])/2
        temp2 = (pos[i, who][who2+1] + target[who2])/2
        yy = rep(i, length(who2)) - legh/2
        segments(temp1, yy, temp2, yy)
      }

      # Add a question mark for those of unknown zygosity
      if (any(plist$twins[i,who] == 3)) {
        who2 = which(plist$twins[i,who] == 3)
        temp1 = (pos[i, who][who2] + target[who2])/2
        temp2 = (pos[i, who][who2+1] + target[who2])/2
        yy = rep(i, length(who2)) - legh/2
        text.default((temp1+temp2)/2, yy, '?')
      }

      # Add the horizontal line
      segments(min(target), i-legh, max(target), i-legh)

      # Draw line to parents. MDV: `pconnect` set to 0.5
      if (diff(range(target)) < 2*pconnect)
        x1 = mean(range(target))
      else
        x1 = pmax(min(target) + pconnect, pmin(max(target) - pconnect, parentx))

      # MDV: `branch` set to 0.6
      y1 = i - legh
      y2 = parentRow + boxh/2
      x2 = parentx
      ydelta = ((y2 - y1) * branch)/2
      segments(c(x1, x1, x2), c(y1, y1 + ydelta, y2 - ydelta),
               c(x1, x2, x2), c(y1 + ydelta, y2 - ydelta, y2))
    }
  } ## end of parent-child lines


  # Duplication arcs
  arcconnect = function(x, y) {
    xx = seq(x[1], x[2], length = 15)
    yy = seq(y[1], y[2], length = 15) + (seq(-7, 7))^2/98 - .5
    lines(xx, yy, lty = 2)
  }

  for (id in nid[duplicated.default(nid, incomparables = 0)]) { # faster than unique
    indx = which(nid == id)
    if (length(indx) > 1) {  # subject is a multiple
      tx = pos[indx]
      ty = row(pos)[indx]

      # MDV: Clarify code. Connect sequentially left -> right
      ord = order(tx)
      tx = tx[ord]
      ty = ty[ord]
      for (j in 1:(length(indx) - 1))
        arcconnect(tx[j + 0:1], ty[j + 0:1])
    }
  }

  ## Finish
  ckall = seq_len(n)[-nid]
  if(length(ckall>0))
    cat('Did not plot the following people:', ckall,'\n')

}


#' @rdname plotmethods
#' @importFrom graphics segments points text.default
#' @export
.annotatePed = function(alignment, annotation, scaling, font = NULL, fam = NULL,
                        col = NULL, colUnder = 1, colInside = 1, colAbove = 1,
                        cex.main = NULL, font.main = NULL, col.main = NULL, line.main = NA, ...) {

  nInd = alignment$nInd
  xall = alignment$xall
  yall = alignment$yall
  plotord = alignment$plotord

  boxh = scaling$boxh
  boxw = scaling$boxw
  labh = scaling$labh
  cex = scaling$cex

  title = annotation$title
  deceased = annotation$deceasedTF
  carrier = annotation$carrierTF
  proband = annotation$probandTF
  textUnder = annotation$textUnder
  textAnnot = annotation$textAnnot
  textInside = annotation$textInside
  textAbove = annotation$textAbove
  col = annotation$colvec

  # Add title
  if(!is.null(title)) {
    title(title, cex.main = cex.main, col.main = col.main,
          font.main = font.main, line = line.main, family = fam, xpd = NA)
  }

  # Deceased
  if(any(deceased)) {
    idx = which(deceased[plotord])
    ids = plotord[idx]
    segments(xall[idx] - .6*boxw, yall[idx] + 1.1*boxh,
             xall[idx] + .6*boxw, yall[idx] - 0.1*boxh, col = col[ids])
  }

  # Carrier dots
  if(any(carrier)) {
    idx = which(carrier[plotord])
    ids = plotord[idx]
    points(xall[idx], yall[idx] + boxh/2, pch = 16, cex = cex, col = col[ids])
  }

  # Proband arrow
  if(any(proband)) {
    pos.arrow = "bottomleft" # Hard coded for now

    mod = switch(pos.arrow,
           bottomleft = list(x = -1, y = 1),
           bottomright = list(x = 1, y = 1),
           topleft = list(x = -1, y = 0),
           topright = list(x = 1, y = 0))

    idx = which(proband[plotord])
    ids = plotord[idx]
    corner.x = xall[idx] + .5*mod$x * boxw
    corner.y = yall[idx] +    mod$y * boxh
    arrows(x0 = corner.x + 1.7*mod$x * boxw,
           y0 = corner.y + 0.9*mod$y * boxh - 0.9*(1-mod$y) * boxh,
           x1 = corner.x + 0.5*mod$x * boxw,
           y1 = corner.y,
           lwd = 1.2, length = .15)
  }

  # Colour vector
  if (length(col) == 1)
    col = rep(col, nInd)

  # Main labels
  text.default(xall, yall + boxh + labh * 0.7, textUnder[plotord], col = colUnder,
       cex = cex, adj = c(.5, 1), font = font, family = fam, xpd = NA)

  # Text inside symbols
  if(!is.null(textInside)) {
    text.default(xall, yall + boxh/2, labels = textInside[plotord], cex = cex, col = colInside,
         font = font, family = fam)
  }

  # Text above symbols
  if(!is.null(textAbove)) {
    if(is.null(font) && any(startsWith(textAbove, "f =")))
      fontAbove = 3
    else
      fontAbove = font
    text.default(xall, yall, labels = textAbove[plotord], cex = cex, col = colAbove,
         family = fam, font = fontAbove, pos = 3, offset = 0.5, xpd = NA)
  }


  if(!is.null(textAnnot)) {
    .addTxt(textAnnot[["topleft"]],     xall-boxw/2, yall,        pos = 2, plotord)
    .addTxt(textAnnot[["top"]],         xall,        yall,        pos = 3, plotord)
    .addTxt(textAnnot[["topright"]],    xall+boxw/2, yall,        pos = 4, plotord)
    .addTxt(textAnnot[["left"]],        xall-boxw/2, yall+boxh/2, pos = 2, plotord)
    .addTxt(textAnnot[["right"]],       xall+boxw/2, yall+boxh/2, pos = 4, plotord)
    .addTxt(textAnnot[["bottomleft"]],  xall-boxw/2, yall+boxh,   pos = 2, plotord)
    .addTxt(textAnnot[["bottom"]],      xall,        yall+boxh,   pos = 1, plotord)
    .addTxt(textAnnot[["bottomright"]], xall+boxw/2, yall+boxh,   pos = 4, plotord)
    .addTxt(textAnnot[["inside"]],      xall,        yall+boxh/2, pos = NULL, plotord)
  }
}


.addTxt = function(args, x, y, pos, plotord) {
  if(is.null(args))
    return()
  txt = args[[1]][plotord]
  do.call(text.default, c(list(x=x,y=y,labels=txt,pos=pos), args[-1]))
}

# Function fixing pedigree alignment of 3/4-siblings and similar
# Founders with two (or more) spouses on the same level should be placed between
.fix34 = function(x, k2ped, plist = NULL, packed = TRUE, width = 10, align = c(1.5, 2)) {

  # Large pedigrees: return unchanged
  if(length(x$ID) > 30)
    return(plist)

  fouInt = founders(x, internal = TRUE)
  nid = plist$nid

  # If no duplicated founders, return unchanged
  dups = duplicated.default(nid, incomparables = 0)
  if(!length(.myintersect(fouInt, nid[dups])))
     return(plist)

  # List of spouses
  ALLSP = vector(mode = "list", length = length(x$ID))
  ALLSP[fouInt] = lapply(fouInt, function(i) spouses(x, i, internal = TRUE))

  # Founders with multiple spouses
  fou2 = fouInt[lengths(ALLSP[fouInt]) > 1]

  # Go row by row in nid
  SP = NULL
  for(k in 2:length(plist$n)) {
    rw = nid[k, ]

    for(id in .myintersect(fou2, rw)) {
      s = .myintersect(ALLSP[[id]], rw) # spouses on that level
      if(length(s) > 1)
        SP = rbind(SP, c(s[1], id, 0), c(id, s[2], 0))
    }
  }

  # If hints added, redo alignment
  if(!is.null(SP)) {
    hints = list(order = seq_along(x$ID), spouses = SP)
    plist = kinship2::align.pedigree(k2ped, packed = packed, width = width, align = align, hints = hints)
  }

  plist
}

# Convert plot parameter (col/fill/lty/lwd) to full vector in pedigree order
.prepPlotarg = function(x, par, default) {
  nInd = length(x$ID)
  nms = names(par)

  if(!is.list(par)) {
    if(!is.null(nms)) {
      vec = rep(default, length = nInd)
      ids = intersect(x$ID, nms)
      vec[internalID(x, ids)] = par[ids]
    }
    else {
      vec = rep(par, length = nInd)
    }
  }
  else { # E.g. list(red = 1:2, "3" = males)
    vec = rep(default, nInd)
    for(cc in nms) {
      v = par[[cc]]
      if(is.function(v))
        ids = v(x)
      else
        ids = intersect(x$ID, v)

      idsInt = internalID(x, ids)
      if(length(idsInt))
        vec[idsInt] = cc
    }
  }

  vec
}


.spouseOrder = function(x, plotorder) {
  if(!is.ped(x))
    stop2("Spouse ordering is not implemented for ped lists")

  if(!is.list(plotorder))
    plotorder = list(plotorder)

  allPairs = list()
  for(ids in plotorder) {
    idsInt = internalID(x, ids)
    newPairs = lapply(2:length(idsInt), function(i) idsInt[(i-1):i])
    allPairs = c(allPairs, newPairs)
  }

  for(p in allPairs) {
    if(!p[1] %in% spouses(x, p[2], internal = TRUE))
      stop2(sprintf("'%s' is not spouse of '%s'", x$ID[p[2]], x$ID[p[1]]))
  }

  spouse = cbind(do.call(rbind, allPairs), 0)
  list(order = seq_along(x$ID), spouse = spouse)
}
