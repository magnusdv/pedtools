
format.marker = function(x, sep = "/", missing = "-", ...) {
  als = c(missing, alleles(x))
  al1 = als[x[, 1] + 1]
  al2 = als[x[, 2] + 1]
  gt = paste(al1, al2, sep=sep)
  if (is_Xmarker(x)) {
    sex = attr(x, 'sex')
    gt[sex == 1] = al1[sex == 1]
  }
  gt
}

print.marker = function(x, sep = "/", missing = "-", ...) {
  gt = format(x, sep=sep, missing=missing)
  df = data.frame(id = attr(x, 'pedmembers'), geno=gt)
  print(df, row.names=FALSE)

  cat("---------\n")
  # Name
  cat("Name:", attr(x, 'name'), "\n")

  # Chromosome - position
  chr = attr(x, 'chrom')
  mb = attr(x, 'posMb')
  cm = attr(x, 'posCm')
  if(!is.na(mb) && !is.na(cm))
    pos = sprintf("%f (mb); %f (cm)", mb, cm)
  else if(!is.na(mb))
    pos = sprintf("%f (Mb)", mb)
  else if(!is.na(cm))
    pos = sprintf("%f (cM)", cm)
  else
    pos = NA
  chr_pos = if(is.na(chr) && is.na(pos)) NA else sprintf("%s - %s", chr, pos)
  cat(sprintf("Position: %s\n", chr_pos))

  # Mutations
  mut = attr(x, "mutmat")
  possible = if(is.null(mut)) "No" else "Yes"
  cat("Mutations possible:", possible, "\n")

  # Allele freqs
  cat("Allele frequencies:\n")
  afr = attr(x, 'afreq')
  names(afr) = attr(x, 'alleles')
  print(afr)
}

.prettyMarkers = function(m, alleles = NULL, sep = "/", missing = NULL, singleCol = FALSE, sex) {
  if (is.null(m))
    return(m)
  if (is.matrix(m))
    m = list(m)
  if ((n <- length(m)) == 0)
    return(m)
  if (is.null(alleles))
    alleles = lapply(m, attr, "alleles")
  else {
    if (!is.atomic(alleles))
      stop("The parameter 'alleles' must be NULL, or an atomic vector.")
    if (length(alleles) < max(unlist(lapply(m, attr, "nalleles"))))
      stop("The indicated 'alleles' vector has too few alleles for some markers.")
    alleles = rep(list(alleles), n)
    }
  if (is.null(missing)) missing = 0

  mNames = unlist(lapply(m, attr, "name"))
  mNames[is.na(mNames)] = ""
  pretty_m = do.call(c, lapply(seq_len(n), function(i)
    c(missing, alleles[[i]])[m[[i]] + 1]))
  dim(pretty_m) = c(length(pretty_m)/(2 * n), 2 * n)

  if (singleCol) {
    al1 = pretty_m[, 2 * seq_len(n) - 1, drop = FALSE]
    al2 = pretty_m[, 2 * seq_len(n), drop = FALSE]
    m.matrix = matrix(paste(al1, al2, sep = sep), ncol = n)
    chrom_X = unlist(lapply(m, function(mm) identical(23L, as.integer(attr(mm, "chrom")))))
    if (any(chrom_X)) {
      males = (sex == 1)
      if (!all(hh <- al1[males, chrom_X] == al2[males, chrom_X]))
        warning("Male heterozygosity at X-linked marker detected.")
      m.matrix[males, chrom_X] = al1[males, chrom_X]
    }
    colnames(m.matrix) = mNames
    return(m.matrix)
  } else {
    nam = rep(mNames, each = 2)
    nam[nzchar(nam)] = paste(nam[nzchar(nam)], 1:2, sep = "_")
    colnames(pretty_m) = nam
    return(pretty_m)
  }
}


