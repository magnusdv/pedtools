# Modified from `kinship2::plot.pedigree()` (GPL >=2 licence).
# The main changes are:
#
# * Separate setup (plot dimensions and scaling) from the actual drawing
# * Fixed scaling bugs, documented here: https://github.com/mayoverse/kinship2/pull/10
# * Adjust scaling to account for duplication arrows
# * Fix unwanted duplications, e.g. in 3/4 siblings
#
# * Removed unused features, especially i) multiple phenotypes and ii) subregions
# * Fixated some parameters at their default values:
#     - pconnect = 0.5
#     - branch = 0.6
#     - align = c(1.5, 2)
#     - packed = TRUE



# Alignment ---------------------------------------------------------------

#' @importFrom kinship2 align.pedigree
.alignPed = function(x, dag = FALSE, packed = TRUE, width = 10, align = c(1.5, 2),
                     hints = NULL, twins = NULL, ...) {

  if(!is.ped(x))
    stop2("First argument to `alignPed()` must be a `ped` object")
  nInd = pedsize(x)

  # Singleton
  if(nInd == 1) {
    plist = list(n = 1, nid = cbind(1), pos = cbind(0), fam = cbind(0), spouse = cbind(0))
    return(list(plist = plist, x = 0, y = 1, nInd = 1, sex = x$SEX, ped = x,
                plotord = 1, xall = 0, yall = 1, maxlev = 1, xrange = c(0,0)))
  }

  # Alignment for DAG presentation (with arrows)
  if(dag) {
    # Generation number of each
    gvec = generations(x, maxOnly = FALSE)
    gvec = as.integer(gvec) # remove names

    # Founders: Place with earliest spouse
    fou = founders(x, internal = TRUE)
    for(f in fou) {
      sp = spouses(x, f, internal = TRUE)
      if(length(sp))
        gvec[f] = min(gvec[sp])
    }

    maxlev = max(gvec)
    glist = lapply(1:maxlev, function(i) which(gvec == i))

    n = lengths(glist)
    maxnum = max(n)

    # Transform glist -> nid
    nid = pos = matrix(0, nrow = maxlev, ncol = maxnum)
    for(i in seq_along(glist)) {
      v = glist[[i]]
      lenv = length(v)
      nid[i, 1:lenv] = v
      pos[i, 1:lenv] = seq.default(from = 0.5*(maxnum - lenv), length.out = lenv)
    }

    plist = list(n = n, nid = nid, pos = pos, fam = NULL, spouse = NULL)
  }
  else {
    # Otherwise: align with kinship2 alignment
    k2ped = as_kinship2_pedigree(x, twins = twins)
    plist = align.pedigree(k2ped, packed = packed, width = width, align = align, hints = hints)

    # Fix annoying rounding error in first column of `pos`
    plist$pos[, 1] = round(plist$pos[, 1], 6)

    # Ad hoc fix for 3/4 siblings and similar
    if(is.null(hints))
      plist = .fix34(x, k2ped, plist, packed, width, align)
  }

  nid = plist$nid
  pos = plist$pos

  xrange = range(pos[nid > 0])
  maxlev = nrow(pos)

  id = as.vector(nid)
  plotord = id[id > 0]

  # Coordinates (top center)
  xall = pos[id > 0]
  yall = row(pos)[id > 0]

  # Kinship2 order (1st instance of each only!)
  # Including this for completeness
  tmp = match(1:nInd, nid)
  xpos = pos[tmp]
  ypos = row(pos)[tmp]

  list(plist = plist, x = xpos, y = ypos, nInd = nInd, sex = x$SEX, ped = x,
       plotord = plotord, xall = xall, yall = yall, maxlev = maxlev, xrange = xrange)
}


#--- Plot dimension and scaling parameters

#' @importFrom graphics frame strheight strwidth
plotSetup = function(pdat0, textUnder = NULL, textAbove = NULL,
                     hasTitle = FALSE, cex = 1, symbolsize = 1, mar = c(1,1,1,1), ...) {

  maxlev = pdat0$maxlev
  xrange = pdat0$xrange
  nid = pdat0$plist$nid

  # Create or advance frame
  frame()

  # Margins
  if(length(mar) == 1)
    mar = if(hasTitle) c(mar, mar, mar + 2.1, mar) else rep(mar, 4)

  # Set basic parameters (usr comes later!)
  oldpar = par(mar = mar, xpd = TRUE)

  psize = par('pin')  # plot region in inches
  stemp1 = strwidth("ABC", units='inches', cex=cex) * 2.5/3
  stemp2 = strheight('1g', units='inches', cex=cex)
  stemp3 = max(strheight(textUnder, units='inches', cex=cex))

  # Text above 1st-generation?
  stemp4 = max(strheight(textAbove[nid[1, ]], units='inches', cex=cex))

  ht1 = psize[2]/maxlev - (stemp2 + stemp3 + stemp4)
  if (ht1 <= 0)
    stop2("Labels leave no room for the graph, reduce cex")

  # Singleton and selfing towers
  if(xrange[1] == xrange[2]) {
    wd2 = psize[1] * 0.5
    boxsize = symbolsize * min(stemp1, wd2) # don't use ht1 here
    hscale = psize[1]
    vscale = (psize[2] - (stemp2 + stemp3 + stemp4 + boxsize)) / (max(1, maxlev - 1))

    yrange = if(maxlev == 1) c(0.5, 1.5) else c(1, maxlev)
    top    = yrange[1] - stemp4/vscale
    bottom = yrange[2] + (stemp2 + stemp3 + boxsize)/vscale

    left  = xrange[1] - 0.5
    right = xrange[2] + 0.5

  }
  else {
    ht2 = psize[2]/(maxlev + (maxlev-1)/2)
    wd2 = psize[1] * 0.8/(.8 + diff(xrange))
    boxsize = symbolsize * min(ht1, ht2, stemp1, wd2) # box size in inches
    hscale = (psize[1] - boxsize)/diff(xrange)  #horizontal scale from user-> inch

    # MDV: Adjust top for curved duplication lines
    nid1 = nid[1, ]
    nid1 = nid1[nid1 > 0]
    curvAdj = if(anyDuplicated.default(nid1)) 0.5 else if(any(nid1 %in% nid[2, ])) 0.1225 else 0

    # Don't adjust for text above if also for curve. (NB: Fails in extreme cases)
    if(curvAdj > 0)
      stemp4 = 0

    denom = maxlev - 1 + curvAdj
    vscale = (psize[2] - (stemp2 + stemp3 + boxsize + stemp4)) / denom

    left   = xrange[1] - 0.5*boxsize/hscale
    right  = xrange[2] + 0.5*boxsize/hscale

    top    = 1 - stemp4/vscale - curvAdj
    bottom = maxlev + (stemp2 + stemp3 + boxsize)/vscale
  }

  usr = c(left, right, bottom, top)

  boxw = boxsize/hscale  # box width in user units
  boxh = boxsize/vscale  # box height in user units
  labh = stemp2/vscale   # height of a text string
  legh = min(1/4, boxh * 1.5)  # how tall are the 'legs' up from a child

  # Return plotting/scaling parameters
  params = list(boxw = boxw, boxh = boxh, labh = labh, legh = legh, usr = usr, oldpar = oldpar)

  # Add params to alignment spec
  c(pdat0, params)
}


#' @importFrom graphics lines polygon segments
.plotPed = function(plotdata, dag = FALSE, aff = 0L, density = NULL, angle = NULL, col = 1, ...) {

  if(dag)
    return(.plotDAG(plotdata, aff = aff, density = density, angle = angle, col = col, ...))

  branch = 0.6
  pconnect = .5

  n = plotdata$nInd
  sex = plotdata$sex
  plist = plotdata$plist
  pos = plist$pos
  nid = plist$nid
  boxh = plotdata$boxh
  boxw = plotdata$boxw
  legh = plotdata$legh
  plotord = plotdata$plotord
  xall = plotdata$xall
  yall = plotdata$yall
  maxlev = nrow(pos)

  # Colour vector
  if (length(col) == 1)
    col = rep(col, n)

  # Aff vector
  if (length(aff) == 1)
    aff = rep(aff, n)

  # Set user coordinates
  par(usr = plotdata$usr, xpd = TRUE)

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
            sex = sex[id], aff = aff[id],
            col = col[id])
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
      if (is.null(plist$twins)) target = pos[i,who]
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


  ## Arcs for multiple instances of subj
  arcconnect = function(x, y) {
    xx = seq(x[1], x[2], length = 15)
    yy = seq(y[1], y[2], length = 15) + (seq(-7, 7))^2/98 - .5
    lines(xx, yy, lty = 2)
  }

  #for (id in unique.default(nid[nid>0])) {
  for (id in nid[duplicated.default(nid, incomparables = 0)]) { # MDV: faster
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

  ## Finish/Final
  ckall = seq_len(n)[-nid]
  if(length(ckall>0))
    cat('Did not plot the following people:', ckall,'\n')

}

#' @importFrom graphics segments
.annotatePed = function(plotdat, deceased = FALSE, carrier = FALSE, title = NULL,
                         textUnder = NULL, textInside = NULL, textAbove = NULL,
                         col = 1, cex = 1, font = NULL, fam = NULL, colUnder = 1, colInside = 1,
                         colAbove = 1, cex.main = cex, font.main = NULL, col.main = NULL, ...) {

  plotord = plotdat$plotord
  xall = plotdat$xall
  yall = plotdat$yall
  boxh = plotdat$boxh
  boxw = plotdat$boxw

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
    col = rep(col, plotdat$nInd)

  # Main labels
  colUnder = rep(colUnder, )
  text(xall, yall + plotdat$boxh + plotdat$labh * 0.7, textUnder[plotord], col = colUnder,
       cex = cex, adj = c(.5, 1), font = font, family = fam, xpd = NA)

  # Add title
  if(!is.null(title)) {
    title(title, cex.main = cex.main, col.main = col.main,
          font.main = font.main, family = fam, xpd = NA)
  }

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
