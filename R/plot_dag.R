# Alignment (corresp to kinship2::plist) for DAG plots (with arrows)
.alignDAG = function(x) {

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

  .extendPlist(x, plist, arrows = TRUE)
}


.plotDAG = function(alignment, annotation, scaling) {

  n = alignment$nInd
  sex = alignment$sex
  plotord = alignment$plotord
  xall = alignment$xall
  yall = alignment$yall

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

  # Symbol radii
  boxw = scaling$boxw
  boxh = scaling$boxh
  rx = boxw/2
  ry = boxh/2

  branch = 0.6
  pconnect = .5

  # Shapes
  POLYS = list(list(x = c(0, -0.5, 0, 0.5), y = c(0, 0.5, 1, 0.5)), # diamond
               list(x = c(-1, -1, 1, 1)/2, y = c(0, 1, 1, 0)),      # square
               list(x = 0.5 * cos(seq(0, 2 * pi, length = 50)),     # circle
                    y = 0.5 * sin(seq(0, 2 * pi, length = 50)) + 0.5))

  # Draw symbols
  for(k in seq_along(plotord)) {
    id = plotord[k]
    poly = POLYS[[sex[id] + 1]]
    dens = if(DENS[id] == 0) NULL else DENS[id]
    polygon(xall[k] + poly$x * boxw,
            yall[k] + poly$y * boxh,
            border = COL[id], col = FILL[id],
            lty = LTY[id], lwd = LWD[id],
            angle = 45, density = dens)
  }

  # Arrows
  x = alignment$ped
  nonf = nonfounders(x, internal = TRUE)
  ch = match(nonf, plotord)
  fa = match(x$FIDX[nonf], plotord)
  mo = match(x$MIDX[nonf], plotord)
  shortArrows(xall[fa],yall[fa], xall[ch],yall[ch],sex0=1,sex1=sex[nonf],rx=rx,ry=ry)
  shortArrows(xall[mo],yall[mo], xall[ch],yall[ch],sex0=2,sex1=sex[nonf],rx=rx,ry=ry)
}


#' @importFrom graphics arrows
shortArrows = function(x0, y0, x1, y1, sex0, sex1, rx, ry, sex, length = 0.1, ...) {
  # NB: Points are given for top center of symbol!
  # Adjusting to mid center
  y0 = y0 + ry
  y1 = y1 + ry

  n = length(x0)
  sex0 = rep_len(sex0, n)
  sex1 = rep_len(sex1, n)

  # Slope
  a = (y1-y0)/(x1-x0) # slopes

  # Males
  xfirst = abs(a) <= ry/rx
  yfirst = !xfirst

  # Females: Params at intersection (rx*cos(t), ry*sin(t))
  t0 = atan(a * rx / ry)
  t0[t0<0] = t0[t0<0] + pi  # bottom intersection
  t1 = t0 - pi

  # Shrink
  fromMx = sex0 == 1 & xfirst
  if(any(fromMx)) {
    x0[fromMx] = x0[fromMx] + rx * sign(a)[fromMx]
    y0[fromMx] = y0[fromMx] + rx * abs(a)[fromMx]
  }

  fromMy = sex0 == 1 & yfirst
  if(any(fromMy)) {
    x0[fromMy] = x0[fromMy] + ifelse(a[fromMy] < Inf, ry/a[fromMy], 0)
    y0[fromMy] = y0[fromMy] + ry
  }

  toMx = sex1 == 1 & xfirst
  if(any(toMx)) {
    x1[toMx] = x1[toMx] - rx * sign(a)[toMx]
    y1[toMx] = y1[toMx] - rx * abs(a)[toMx]
  }

  toMy = sex1 == 1 & yfirst
  if(any(toMy)) {
    x1[toMy] = x1[toMy] - ifelse(a[toMy] < Inf, ry/a[toMy], 0)
    y1[toMy] = y1[toMy] - ry
  }

  # Shrink females (and diamonds)
  if(any(fromF <- sex0 != 1)) {
    x0[fromF] = x0[fromF] + rx * cos(t0[fromF])
    y0[fromF] = y0[fromF] + ry * sin(t0[fromF])
  }
  if(any(toF <- sex1 != 1)) {
    x1[toF] = x1[toF] + rx * cos(t1[toF])
    y1[toF] = y1[toF] + ry * sin(t1[toF])
  }

  arrows(x0, y0, x1, y1, length = length, ...)
}
