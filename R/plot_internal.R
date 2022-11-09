# Modified from `kinship2::plot.pedigree()`.
# The main changes are:
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
# Internally, the previous single `plot()` function has been refactored into the following steps:
# * .pedAlignment(): Builds on `kinship2::align.pedigree(), but also handles singletons and DAGs
# * .pedAnnotation(): Prepare and collect annotation variables
# * .pedScaling(): Calculate symbol sizes and scaling variables
# * .drawPed(): Draw symbols
# * .annotatePed(): Add labels and other annotation




#' Internal plot methods
#'
#' These functions provide the machinery for pedigree plotting. Their main
#' purpose is to be called internally by [plot.ped()], which is the recommended
#' front-end for most users. The various plot options are documented here.
#'
#' The workflow of `plot.ped(x, ...)` is approximately as follows:
#'
#' ```
#'
#' # Calculate plot parameters
#'
#' align = .pedAlignment(x, ...)
#'
#' annot = .pedAnnotation(x, ...)
#'
#' scale = .pedScaling(align, annot, ...)
#'
#' # Produce plot
#'
#' .drawPed(align, annot, scale)
#'
#' .annotatePed(align, annot, scale)
#'
#' ```
#'
#' @param x A [ped()] object.
#' @param plist Alignment list with format similar to
#'   [kinship2::align.pedigree()].
#' @param arrows A logical (default = FALSE). If TRUE, the pedigree is plotted
#'   as a DAG, i.e., with arrows connecting parent-child pairs.
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
#' @param cex Expansion factor controlling font size. This also affects symbol
#'   sizes, which by default have the width of 2.5 characters. Default: 1.
#' @param symbolsize Expansion factor for pedigree symbols. Default: 1.
#' @param margins A numeric indicating the plot margins. If a single number is
#'   given, it is recycled to length 4.
#' @param addSpace A numeric of length 4, indicating extra padding (in inches)
#'   around the pedigree inside the plot region. Default: 0.
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
#' @param packed,width,align Parameters passed on to
#'   [kinship2::align.pedigree()]. Can usually be left untouched.
#' @param hints An optional list of hints passed on to
#'   [kinship2::align.pedigree()].
#' @param fouInb Either "autosomal" (default), "x" or NULL. If "autosomal" or
#'   "x", inbreeding coefficients are added to the plot above the inbred
#'   founders. If NULL, or if no founders are inbred, nothing is added.
#' @param textInside,textAbove Character vectors of text to be printed inside or
#'   above pedigree symbols.
#' @param font,fam Arguments passed on to [text()].
#' @param colUnder,colInside,colAbove Colour vectors.
#' @param cex.main,col.main,font.main Parameters passed on to [title()].
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
#' frame()
#' drawPed(align, annot, scale)
#'
#' @name internalplot
NULL




# Alignment ---------------------------------------------------------------

#' @rdname internalplot
#' @importFrom kinship2 align.pedigree
#' @export
.pedAlignment = function(x = NULL, plist = NULL, arrows = FALSE, twins = NULL, packed = TRUE,
                         width = 10, align = c(1.5, 2), hints = NULL, ...) {

  if(hasSelfing(x) && !arrows) {
    message("Pedigree has selfing, switching to DAG mode. Use `arrows = TRUE` to avoid this message.")
    arrows = TRUE
  }

  if(!is.null(plist))
    return(.extendPlist(x, plist, arrows))

  # Singleton
  if(is.singleton(x)) {
    plist = list(n = 1, nid = cbind(1), pos = cbind(0), fam = cbind(0), spouse = cbind(0))
    return(.extendPlist(x, plist))
  }

  if(arrows)
    return(.alignDAG(x))

  # Twin data: enforce data frame
  if(is.vector(twins))
    twins = data.frame(id1 = twins[1], id2 = twins[2], code = as.integer(twins[3]))

  k2ped = as_kinship2_pedigree(x, twins = twins)
  plist = kinship2::align.pedigree(k2ped, packed = packed, width = width, align = align, hints = hints)

  # Ad hoc fix for 3/4 siblings and similar
  if(is.null(hints))
    plist = .fix34(x, k2ped = k2ped, plist = plist, packed = packed, width = width, align = align)

  # Fix annoying rounding errors in first column of `pos`
  plist$pos[] = round(plist$pos[], 6)

  # Add further parameters
  .extendPlist(x, plist)
}


.extendPlist = function(x, plist, arrows = FALSE) {
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

  list(plist = plist, x = xpos, y = ypos, nInd = nInd, sex = getSex(x), ped = x, arrows = arrows,
       plotord = plotord, xall = xall, yall = yall, maxlev = maxlev, xrange = xrange, yrange = yrange)
}


# Annotation --------------------------------------------------------------

#' @rdname internalplot
#' @export
.pedAnnotation = function(x, title = NULL, marker = NULL, sep = "/", missing = "-",
                          showEmpty = FALSE, labs = labels(x), col = 1, aff = NULL, carrier = NULL,
                          hatched = NULL, deceased = NULL, starred = NULL,
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
    if (!showEmpty)
      geno[rowSums(do.call(cbind, mlist)) == 0] = ""

    textu = if (!any(nzchar(textu))) geno else paste(textu, geno, sep = "\n")
  }

  res$textUnder = trimws(textu, "right", whitespace = "[\t\r\n]")


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


#--- Plot dimension and scaling parameters

#' @rdname internalplot
#' @importFrom graphics frame strheight strwidth
#' @export
.pedScaling = function(alignment, annotation, cex = 1, symbolsize = 1, margins = 1, addSpace = 0, ...) {

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

  if (ht1 <= 0)
    stop2("Labels leave no room for the graph, reduce cex")

  # Height restriction 2
  ht2 = psize[2]/(maxlev + (maxlev-1)/2)

  # Width restriction 1: Default width = 2.5 letters (`stemp1`)
  wd1 = strwidth("ABC", units='inches', cex=cex) * 2.5/3

  # Width restriction 2
  wd2 = psize[1] * 0.8/(.8 + diff(xrange))  # = psize[1] for singletons/selfings

  # Box size in inches
  boxsize = symbolsize * min(ht1, ht2, wd1, wd2)

  # Segments of length 1 inch
  hscale = (psize[1] - boxsize - addSpace[2] - addSpace[4])/diff(xrange)
  denom = if(maxlev == 1) 1 else maxlev - 1 + curvAdj
  vscale = (psize[2] - (abovetop_in + boxsize + belowlast_in)) / denom

  if(hscale <= 0 || vscale <= 0)
    stop2("Cannot fit the graph; please increase plot region or reduce cex and/or symbolsize")

  boxw = boxsize/hscale  # box width in user units
  boxh = boxsize/vscale  # box height in user units

  # User coordinates
  left   = xrange[1] - 0.5*boxw - addSpace[2]/hscale
  right  = xrange[2] + 0.5*boxw + addSpace[4]/hscale
  top    = yrange[1] - abovetop_in/vscale - curvAdj
  bottom = yrange[2] + (boxsize + belowlast_in)/vscale
  usr = c(left, right, bottom, top)

  labh = labh_in/vscale        # height of a text string
  legh = min(1/4, boxh * 1.5)  # how tall are the 'legs' up from a child

  # Return plotting/scaling parameters
  list(boxw = boxw, boxh = boxh, labh = labh, legh = legh,
       cex = cex, mar = mar, usr = usr, oldpar = oldpar)
}


#' @rdname internalplot
#' @importFrom graphics lines polygon segments
#' @export
.drawPed = function(alignment, annotation, scaling) {

  if(isTRUE(alignment$arrows))
    return(.plotDAG(alignment, annotation, scaling))

  AFF = annotation$aff01 %||% 0
  COL = annotation$colvec %||% 1
  density = annotation$density %||% -1
  angle = annotation$angle %||% 90

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

  branch = 0.6
  pconnect = .5

  # Colour vector
  if (length(COL) == 1)
    COL = rep(COL, n)

  # Aff vector
  if (length(AFF) == 1)
    AFF = rep(AFF, n)

  # Set user coordinates
  par(mar = scaling$mar, usr = scaling$usr, xpd = TRUE)

  # Shapes
  POLYS = list(list(x = c(0, -0.5, 0, 0.5), y = c(0, 0.5, 1, 0.5)), # diamond
               list(x = c(-1, -1, 1, 1)/2, y = c(0, 1, 1, 0)),      # square
               list(x = 0.5 * cos(seq(0, 2 * pi, length = 50)),     # circle
                    y = 0.5 * sin(seq(0, 2 * pi, length = 50)) + 0.5))

  # Function for drawing a single symbol
  drawbox = function(xpos, ypos, sex, aff, col) {
    poly = POLYS[[sex + 1]]

    polygon(xpos + poly$x * boxw, ypos + poly$y * boxh, border = col,
            col = if(aff == 1) col else NA, angle = angle,
            density = if(aff == 1) density else NULL)
  }

  # Draw all the symbols
  for (k in seq_along(plotord)) {
    id = plotord[k]
    drawbox(xpos = xall[k], ypos = yall[k],
            sex = SEX[id], aff = AFF[id],
            col = COL[id])
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
  for(i in seq_len(maxlev-1) + 1) {       # MDV: safer than 2:maxlev
    zed = unique.default(plist$fam[i,  ]) # MDV: use unique.default
    zed = zed[zed > 0]  #list of family ids

    for(fam in zed) {
      xx = pos[i - 1, fam + 0:1]
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
        text((temp1+temp2)/2, yy, '?')
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
      y2 = (i-1) + boxh/2
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


#' @rdname internalplot
#' @importFrom graphics segments points text
#' @export
.annotatePed = function(alignment, annotation, scaling, font = NULL, fam = NULL,
                        col = NULL, colUnder = 1, colInside = 1, colAbove = 1,
                        cex.main = NULL, font.main = NULL, col.main = NULL, ...) {

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
  textUnder = annotation$textUnder
  textInside = annotation$textInside
  textAbove = annotation$textAbove
  col = annotation$colvec

  # Add title
  if(!is.null(title)) {
    title(title, cex.main = cex.main, col.main = col.main,
          font.main = font.main, family = fam, xpd = NA)
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

  # Colour vector
  if (length(col) == 1)
    col = rep(col, nInd)

  # Main labels
  text(xall, yall + boxh + labh * 0.7, textUnder[plotord], col = colUnder,
       cex = cex, adj = c(.5, 1), font = font, family = fam, xpd = NA)

  # Text inside symbols
  if(!is.null(textInside)) {
    text(xall, yall + boxh/2, labels = textInside[plotord], cex = cex, col = colInside,
         font = font, family = fam)
  }

  # Text above symbols
  if(!is.null(textAbove)) {
    if(is.null(font) && any(startsWith(textAbove, "f =")))
      fontAbove = 3
    else
      fontAbove = font
    text(xall, yall, labels = textAbove[plotord], cex = cex, col = colAbove,
         family = fam, font = fontAbove, pos = 3, offset = 0.5, xpd = NA)
  }
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
    plist = align.pedigree(k2ped, packed = packed, width = width, align = align, hints = hints)
  }

  plist
}


