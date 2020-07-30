
#' @export
format.marker = function(x, sep = "/", missing = "-", ...) {
  als = c(missing, alleles(x))
  al1 = als[x[, 1] + 1]
  al2 = als[x[, 2] + 1]
  gt = paste(al1, al2, sep = sep)
  if (isXmarker(x)) {
    male = attr(x, 'sex') == 1

    # Male 'homozygosity': Just show first allele
    maleHom = male & (al1 == al2)
    gt[maleHom] = al1[maleHom]
  }
  gt
}

#' @export
print.marker = function(x, sep = "/", missing = "-", ...) {
  ids = attr(x, 'pedmembers')
  gt = format(x, sep = sep, missing = missing)

  df = data.frame(id = ids, gt = gt)

  # Use marker name
  mname = name(x)
  names(df)[2] = if(is.na(mname)) "<NA>" else mname

  # Extra space between columns
  df = commentAndRealign(df, 1, rep(TRUE, nrow(df)), " ")

  # If X: add question mark for heterozygous males
  if(isXmarker(x)) {
    maleHet = attr(x, "sex") == 1 & !gt %in% c(alleles(x), missing)
    df = commentAndRealign(df, 2, maleHet, "?")
  }

  print(df, row.names = FALSE)

  ### Other attributes

  # Mutation model
  mut = attr(x, "mutmod")
  muttxt = if(is.null(mut)) "none" else toString(mut)

  # Position
  chr = chrom(x)
  pos = posMb(x)
  postxt = if(is.na(chr) && is.na(pos)) NA else sprintf("chr = %s, Mb = %g", chr, pos)

  # Print info
  cat(strrep("* ", (max(nchar(ids)) + max(nchar(gt)) + 6)/2), "\n")
  cat("Position:", postxt, "\n")
  cat("Mutation:", muttxt, "\n")
  cat("Frequencies:\n")

  # Hack to get one space indentation
  print(data.frame(as.list(afreq(x)), check.names = FALSE), row.names = FALSE)

  invisible(x)
}



