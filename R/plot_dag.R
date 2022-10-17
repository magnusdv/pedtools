

.plotDAG = function(plotdata, aff = NULL, density = NULL, angle = NULL, col = 1, ...) {

  # Set user coordinates
  par(usr = plotdata$usr, xpd = TRUE)

  # Symbol radii
  boxw = plotdata$boxw
  boxh = plotdata$boxh
  rx = boxw/2
  ry = boxh/2

  xall = plotdata$xall
  yall = plotdata$yall

  # Shapes
  POLYS = list(list(x = c(0, -0.5, 0, 0.5), y = c(0, 0.5, 1, 0.5)), # diamond
               list(x = c(-1, -1, 1, 1)/2, y = c(0, 1, 1, 0)),      # square
               list(x = 0.5 * cos(seq(0, 2 * pi, length = 50)),     # circle
                    y = 0.5 * sin(seq(0, 2 * pi, length = 50)) + 0.5))

  sex = plotdata$sex
  plotord = plotdata$plotord

  # Draw symbols
  for(k in seq_along(plotord)) {
    id = plotord[k]
    poly = POLYS[[sex[id] + 1]]
    polygon(xall[k] + poly$x * boxw, yall[k] + poly$y * boxh,
            col = if(aff[id]) col[id] else NA, border = col[id],
            density = if(aff[id]) density else NULL, angle = angle)
  }

  # Arrows
  x = plotdata$ped
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
