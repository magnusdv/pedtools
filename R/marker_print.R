
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
  pedlabels = attr(x, 'pedmembers')
  gt = format(x, sep = sep, missing = missing)
  df = data.frame(pedlabels, gt, stringsAsFactors = FALSE)

  # Use marker name as column header, or "<NA>"
  mname = name(x)
  if(is.na(mname))
    mname = "<NA>"
  names(df) = c("", mname)

  if(isXmarker(x)) {
    maleHet = attr(x, "sex") == 1 & !gt %in% c(alleles(x), missing)
    df = commentAndRealign(df, 2, maleHet, "?")
  }
  print(df, row.names = FALSE)

  cat(strrep("-", max(nchar(pedlabels)) + max(4, nchar(mname)) +3), "\n")

  # Chromosome - position
  chr = chrom(x)
  mb = posMb(x)
  cat(sprintf("Chrom %s: %s (Mb)\n", chr, mb))

  # Mutations
  mut = attr(x, "mutmod")
  if(!is.null(mut))
    muttxt = if(trivialMut(mut)) "Yes (trivial)" else "Yes"
  else muttxt = "No"
  cat("Mutation model:", muttxt, "\n")

  # Allele freqs - use hack to get one space indentation
  cat("Allele frequencies:\n")
  afr = afreq(x)
  print(data.frame(as.list(afr), check.names = FALSE), row.names = FALSE)
}



