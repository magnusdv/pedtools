#' Marker objects
#'
#' Creating a marker object associated with a pedigree
#'
#'
#' @param x a [`ped`] object
#' @param ... one or more expressions of the form `id = genotype`, where `id` is
#'   the ID label of a member of `x`, and genotype is a numeric or character
#'   vector of length 1 or 2 (see Examples).
#' @param alleles a character (or coercible to character) containing allele
#'   names. If not given, the default is to take the sorted vector of distinct
#'   alleles occuring in `...`.
#' @param afreq a numeric of the same length as `alleles`, indicating the
#'   population frequency of each allele. A warning is issued if the frequencies
#'   don't sum to 1 after rounding to 3 decimals. If `afreq` is not specified,
#'   all alleles are given equal frequencies.
#' @param chrom a single integer: the chromosome number. Default: NA.
#' @param posMb a nonnegative real number: the phsyical position of the marker,
#'   in megabases. Default: NA.
#' @param posCm a nonnegative real number: the centiMorgan position of the
#'   marker. Default: NA.
#' @param name a character string: the name of the marker. Default: NA.
#' @param mutmat a mutation matrix, or a list of two such matrices named
#'   'female' and 'male'. If given, the mutation matrices must be square, with
#'   the allele labels as dimnames, and each row must sum to 1 (after rounding
#'   to 3 decimals). Default: NULL.
#'
#' @return An object of class `marker`: This is an integer matrix with 2 columns
#'   and one row per individual, and attributes 'alleles' (a character vector
#'   with allele labels), 'afreq` (allele frequencies), 'chrom' (chromosome
#'   number), 'posMb' (physical location in megabases),'posCm' (position in
#'   centiMorgan), 'name' (marker identifier) and 'mutmat' (a list of two (male
#'   and female) mutation matrices).
#'
#' @seealso [marker_attach]
#'
#' @examples
#' x = nuclearPed(father="fa", mother="mo", children="child")
#'
#' # A rare SNP marker for which the child is heterozygous
#' m = marker(x, child = 1:2, alleles = 1:2, afreq = c(0.01, 0.99))
#'
#' # Sometimes it is useful to attach the marker to the pedigree
#' x = setMarkers(x, m)
#'
#' @export
marker = function(x, ...,  alleles = NULL, afreq = NULL, chrom = NA,
                  posMb = NA, posCm = NA, name = NA, mutmat = NULL) {

  # Initalize empty allele matrix
  m = matrix(0, ncol = 2, nrow = pedsize(x))

  # Capture genotypes given in dots
  dots = eval(substitute(alist(...)))
  if(length(dots) > 0 && is.null(names(dots)))
    stop2("Genotype assignments in `...` must be named. See ?marker")

  ids_int = internalID(x, names(dots))

  genos = lapply(dots, eval.parent)

  for(i in seq_along(dots)) {
    g = genos[[i]]
    if(!is.vector(g) || !length(g) %in% 1:2) # TODO: deal with NA and '' inputs
      stop2("Genotype must be a vector of length 1 or 2: ", deparse(g))
    m[ids_int[i], ] = g
  }

  .createMarkerObject(m, alleles = alleles, afreq = afreq, chrom = chrom,
      posMb = posMb, posCm=posCm, name = name, mutmat = mutmat,
      pedmembers = x$LABELS, sex = x$SEX)
}


# TODO: rename to new_marker() and move checks to validate_marker()
.createMarkerObject = function(matr, alleles = NULL, afreq = NULL, chrom = NA,
                               posMb = NA, posCm = NA, name = NA, mutmat = NULL,
                               pedmembers = NULL, sex = NULL, check_input = TRUE) {
  NA_allele_ = c(0, "", NA, "-")

  ### Checks
  if(check_input) {

    # Check that the observed alleles are OK
    if(!is.null(alleles)) {
      if(any(alleles %in% NA_allele_))
        stop2("Invalid entry in `alleles`: ", intersect(alleles, NA_allele_))
      if(!all(matr %in% c(NA_allele_, alleles))) {
        notfound = setdiff(matr, c(NA_allele_, alleles))
        stop2("Allele used in genotype but not included in `alleles` argument: ",
              notfound)
      }
    }

    # Check frequencies
    if(!is.null(afreq)) {
      if(is.null(alleles))
         stop2("Argument `alleles` cannot be NULL if `afreq` is non-NULL")
      if (length(afreq) != length(alleles))
        stop2("Number of alleles doesn't match length of frequency vector")
      if (round(sum(afreq), 3) != 1)
        stop2("Allele frequencies do not sum to 1: ", afreq)
    }

    # Check name
    if(length(name) != 1)
      stop2("Length of `name` must be 1: ", name)
    if (isTRUE(suppressWarnings(name == as.integer(name))))
      stop2("Attribute `name` cannot consist entirely of digits: ", name)

    # Check chrom
    if(length(chrom) != 1)
      stop2("Length of `chrom` must be 1: ", chrom)

    # Check pedmembers and sex
    if(length(pedmembers) != nrow(matr))
      stop2("Number of rows in the allele matrix must equal the length of `pedmembers`")
    if(length(sex) != nrow(matr))
      stop2("Number of rows in the allele matrix must equal the length of `sex`")

    # Check mutation matrices
    if (!is.null(mutmat)) {
      if (is.matrix(mutmat))
        mutmat = .checkMutationMatrix(mutmat, alleles)
      else if(is.list(mutmat)) {
        nms = names(mutmat)
        if(!setequal(nms, c("female", "male")))
          stop2("List of mutation matrices must have names 'male' and 'female': ", nms)
        mutmat$female = .checkMutationMatrix(mutmat$female, alleles)
        mutmat$male = .checkMutationMatrix(mutmat$male, alleles)
      }
      else
        stop2("Argument `mutmat` must be either a matrix or a list of two matrices
              (named `male` and `female`)")
    }
  }

  ### Alleles
  if (is.null(alleles)) {
    alleles = .mysetdiff(matr, NA_allele_)
    if (!length(alleles)) alleles = 1
  }
  else if (!is.numeric(alleles) && !anyNA(suppressWarnings(as.numeric(alleles))))
    alleles = as.numeric(alleles) # ensure numerical sorting if appropriate

  # Sort (same order used below for frequencies)
  allele_order = order(alleles)
  alleles = as.character(alleles[allele_order])

  ### Frequencies
  if (is.null(afreq)) {
    nall = length(alleles)
    afreq = rep.int(1/nall, nall)
  }
  else
    afreq = afreq[allele_order]

  ### Name
  name = as.character(name)

  ### Chromomsome
  chrom = as.character(chrom)

  ### Mutation matrices
  if (is.matrix(mutmat)) {
    mutmat = list(male = mutmat, female = mutmat)
  }

  ### Create the internal allele matrix
  m = match(matr, alleles, nomatch = 0)

  ### Add everything else as attributes, including class = "marker".
  attributes(m) = list(dim = dim(matr), alleles = alleles, afreq = afreq,
                       chrom = chrom, posMb = posMb, posCm = posCm, name = name,
                       mutmat = mutmat, pedmembers = pedmembers, sex = sex,
                       class = "marker")
  m
}

.checkMutationMatrix = function(mutmat, alleles) {
    # Check that mutation matrix is compatible with allele number / names.  Sort matrix
    # according to the given allele order (this is important since dimnames are NOT used in
    # calculations).
  alleles = as.character(alleles)
  N = length(alleles)

  if (any((dm <- dim(mutmat)) != N)) {
    print(mutmat)
    stop2(sprintf("Dimension of mutation matrix (%d x %d)
                     inconsistent with number of alleles (%d)", dm[1], dm[2], N))
  }

  if (any(round(.rowSums(mutmat, N, N), 3) != 1)) {
    print(mutmat)
    stop2("Row sums of mutation matrix are not 1")
  }

  row_nms = rownames(mutmat)
  col_nms = colnames(mutmat)
  if (!setequal(row_nms, alleles) || !setequal(col_nms, alleles)) {
    print(mutmat)
    stop2("Dimnames of mutation do not match allele names")
  }

  m = mutmat[alleles, alleles]

  # Lumbability: always lumpable (any alleles) if rows are identical (except diagonal)
  attr(m, "lumpability") = NA
  if (N > 2) {
    if (all(vapply(1:N, function(i) diff(range(m[-i, i])) == 0, logical(1))))
      attr(m, "lumpability") = "always"  # If each column has 0 range
  }
  m
}

.setSNPfreqs = function(x, newfreqs) { # TODO: review function
    stopifnot(all(vapply(x$markerdata, function(m) attr(m, "nalleles"), numeric(1)) == 2))
    newfreqs = rep(newfreqs, length = x$nMark)
    for (i in seq_len(x$nMark)) attr(x$markerdata[[i]], "afreq") = c(newfreqs[i], 1 - newfreqs[i])
    x
}
