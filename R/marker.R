#' Marker objects
#'
#' Creating a marker object associated with a pedigree.
#'
#'
#' @param x a [`ped`] object
#' @param ... one or more expressions of the form `id = genotype`, where `id` is
#'   the ID label of a member of `x`, and `genotype` is a numeric or character
#'   vector of length 1 or 2 (see Examples).
#' @param geno a character vector of length `pedsize(x)`, with genotypes
#'   written in the format "a/b".
#' @param allelematrix a matrix with 2 columns and `pedsize(x)` rows. If this is
#'   non-NULL, then `...` must be empty.
#' @param alleles a character (or coercible to character) containing allele
#'   names. If not given, and `afreq` is named, `names(afreq)` is used. The
#'   default action is to take the sorted vector of distinct alleles occurring
#'   in `allelematrix`, `geno` or `...`.
#' @param afreq a numeric of the same length as `alleles`, indicating the
#'   population frequency of each allele. A warning is issued if the frequencies
#'   don't sum to 1 after rounding to 3 decimals. If the vector is named, and
#'   `alleles` is not NULL, an error is raised if `setequal(names(afreq),
#'   alleles)` is not TRUE. If `afreq` is not specified, all alleles are given
#'   equal frequencies.
#' @param chrom a single integer: the chromosome number. Default: NA.
#' @param posMb a nonnegative real number: the physical position of the marker,
#'   in megabases. Default: NA.
#' @param name a character string: the name of the marker. Default: NA.
#' @param NAstrings A character vector containing strings to be treated as
#'   missing alleles. Default: `c("", "0", NA, "-")`.
#' @param mutmod,rate mutation model parameters. These are passed directly to
#'   [pedmut::mutationModel()]; see there for details. Note: `mutmod`
#'   corresponds to the `model` parameter. Default: NULL (no mutation model).
#' @param validate if TRUE, the validity of the created `marker` object is
#'   checked.
#'
#' @return An object of class `marker`. This is an integer matrix with 2 columns
#'   and one row per individual, and the following attributes:
#'
#'   * `alleles` (a character vector with allele labels)
#'
#'   * `afreq` (allele frequencies; default `rep.int(1/length(alleles),
#'   length(alleles))`)
#'
#'   * `chrom` (chromosome number; default = NA)
#'
#'   * `posMb` (physical location in megabases; default = NA)
#'
#'   * `name` (marker identifier; default = NA)
#'
#'   * `mutmod` (a list of two (male and female) mutation matrices; default =
#'   NULL)
#'
#' @seealso [marker_attach]
#'
#' @examples
#' x = nuclearPed(father = "fa", mother = "mo", children = "child")
#'
#' # An empty SNP with alleles "A" and "B"
#' marker(x, alleles = c("A", "B"))
#'
#' # Alleles/frequencies can be given jointly or separately
#' stopifnot(identical(
#'   marker(x, afreq = c(A = 0.01, B = 0.99)),
#'   marker(x, alleles = c("A", "B"), afreq = c(0.01, 0.99)),
#'   ))
#'
#' # Genotypes can be assigned individually ...
#' marker(x, fa = "1/1", mo = "1/2")
#'
#' # ... or using the `geno` vector (all members in order)
#' marker(x, geno = c("1/1", "1/2", NA))
#'
#' # For homozygous genotypes, a single allele suffices
#' marker(x, fa = 1)
#'
#' # Attaching a marker to the pedigree
#' m = marker(x) # By default a SNP with alleles 1,2
#' x = setMarkers(x, m)
#'
#' # A marker with a "proportional" mutation model,
#' # with different rates for males and females
#' mutrates = list(female = 0.1, male = 0.2)
#' marker(x, alleles = 1:2, mutmod = "prop", rate = mutrates)
#'
#' @export
marker = function(x, ...,  geno = NULL, allelematrix = NULL,
                  alleles = NULL, afreq = NULL,
                  chrom = NA, posMb = NA, name = NA,
                  NAstrings = c(0, "", NA, "-"),
                  mutmod = NULL, rate = NULL,
                  validate = TRUE) {

  # Some parameters cannot have length 0 or be ""
  if(length(chrom) == 0) chrom = NA
  if(length(posMb) == 0) posMb = NA
  if(length(name) == 0 || identical(name, "")) name = NA

  pedN = pedsize(x)

  if (is.null(geno) && is.null(allelematrix)) {
    # Initialise empty allele matrix
    m = matrix(0, ncol = 2, nrow = pedN)

    # Capture genotypes given in dots
    dots = eval(substitute(alist(...)))
    if((ld <- length(dots)) > pedN)
      stop2("Too many genotype assignments")

    # Genotypes (may be empty)
    genos = lapply(dots, eval.parent)

    # Internal ID of each genotype
    # (If no names, take pedigree members in sequence)
    dotnames = names(dots)
    if(is.null(dotnames))
      ids_int = seq_len(ld)
    else
      ids_int = internalID(x, dotnames)

    for(i in seq_len(ld)) {
      g = genos[[i]]
      lg = length(g)

      if(!is.vector(g) || !lg %in% 1:2)
        stop2("Genotype must be a vector of length 1 or 2: ", deparse(g))

      # Split compound genotypes, e.g., "a/b"
      if(lg == 1 && is.character(g))
        g = strsplit(g, "/", fixed = TRUE)[[1]]

      # Insert in `m`
      m[ids_int[i], ] = g
    }
  }
  else if(!is.null(geno)) {
    if(!is.null(allelematrix))
      stop2("At least one of `geno` and `allelematrix` must be NULL")
    geno = as.character(geno)
    if(length(geno) != pedN)
      stop2("`geno` incompatible with pedigree")
    s = strsplit(geno, "/")
    s[lengths(s) < 2] = lapply(s[lengths(s) < 2], rep, length.out = 2)
    m = matrix(unlist(s), ncol = 2, byrow = TRUE)
  }
  else
    m = allelematrix

  ### Alleles and frequencies
  if(!is.null(alleles) && !is.null(names(afreq)))
    stop2("Argument `alleles` should not be used when `afreq` has names")
  if(is.null(alleles) && !is.null(afreq) && is.null(names(afreq)))
    stop2("When `alleles` is NULL, `afreq` must be named")


  # If alleles are NULL, take from afreq names, otherwise from supplied genos
  als = alleles %||% names(afreq) %||% .mysetdiff(m, NAstrings)
  if(length(als) == 0)
    als = 1:2

  ### Frequencies
  afreq = afreq %||% {rep_len(1, length(als))/length(als)}
  names(afreq) = names(afreq) %||% als

  # Sort alleles and frequencies (numerical sorting if appropriate)
  if (!is.numeric(als) && !anyNA(suppressWarnings(as.numeric(als))))
    ord = order(as.numeric(als))
  else
    ord = order(als)

  # Final ordered objects
  AFR = afreq[ord]
  ALS = names(AFR)

  ### Mutation model
  if(!is.null(mutmod)) {
    if (!requireNamespace("pedmut", quietly = TRUE))
      stop2("Package `pedmut` must be installed in order to include mutation models")
    mutmod = pedmut::mutationModel(mutmod, alleles = ALS, afreq = AFR, rate = rate)
  }

  ### Internal allele matrix
  if(!all(m %in% c(NAstrings, ALS)))
    stop2("Invalid allele", if(!is.na(name)) sprintf(" for marker `%s`", name), ": ",
          setdiff(m, c(NAstrings, ALS)))

  m_int = match(m, ALS, nomatch = 0, incomparables = NAstrings)
  dim(m_int) = dim(m)

  # Create marker object
  ma = newMarker(m_int, alleles = ALS, afreq = unname(AFR),
            name = as.character(name), chrom = as.character(chrom),
            posMb = as.numeric(posMb), mutmod = mutmod,
            pedmembers = labels(x), sex = x$SEX)

  if(validate)
    validateMarker(ma)
  ma
}


newMarkerOLD = function(allelematrix_int, alleles, afreq, name = NA_character_,
                     chrom = NA_character_, posMb = NA_real_,
                     mutmod = NULL, pedmembers, sex) {

  stopifnot2(is.matrix(allelematrix_int),
            ncol(allelematrix_int) == 2,
            is.integer(allelematrix_int),
            is.character(alleles),
            is.numeric(afreq),
            is.character(name),
            is.character(chrom),
            is.numeric(posMb),
            is.null(mutmod) || is.list(mutmod),
            is.character(pedmembers),
            is.integer(sex))

  structure(allelematrix_int, alleles = alleles, afreq = afreq, name = name,
            chrom = chrom, posMb = posMb, mutmod = mutmod,
            pedmembers = pedmembers, sex = sex, class = "marker")
}

#' Internal marker constructor
#'
#' This is the internal constructor of `marker` objects. It does not do any
#' input validation and should only be used in programming scenarios, and only
#' if you know what you are doing. Most users are recommended to use the regular
#' constructor [marker()].
#'
#' See [marker()] for more details about the marker attributes.
#'
#' @param alleleMatrixInt An integer matrix.
#' @param alleles A character vector.
#' @param afreq A numeric vector.
#' @param name A character of length 1.
#' @param chrom A character of length 1.
#' @param posMb A numeric of length 1.
#' @param mutmod A mutation model.
#' @param pedmembers A character vector.
#' @param sex An integer vector.
#'
#' @return A `marker` object.
#'
#' @examples
#'
#' newMarker(matrix(c(1L, 0L, 1L, 1L, 0L, 2L), ncol = 2),
#'           alleles = c("A", "B"), afreq = c(0.1, 0.9), name = "M",
#'           pedmembers = c("1", "2", "3"), sex = c(1L, 2L, 1L))
#'
#' @export
newMarker = function(alleleMatrixInt, alleles, afreq, name = NA_character_,
                     chrom = NA_character_, posMb = NA_real_,
                     mutmod = NULL, pedmembers, sex) {

  x = alleleMatrixInt
  attributes(x) = list(dim = dim(alleleMatrixInt),
                       alleles = alleles, afreq = afreq, name = name,
                       chrom = chrom, posMb = posMb, mutmod = mutmod,
                       pedmembers = pedmembers, sex = sex)
  class(x) = "marker"
  x
}

validateMarker = function(x) {
  attrs = attributes(x)

  ## alleles
  alleles = attrs$alleles
  NA_allele_ = c(0, "", NA)
  if(any(alleles %in% NA_allele_))
    stop2("Invalid entry in `alleles`: ", intersect(alleles, NA_allele_))

  if(dup <- anyDuplicated(alleles))
    stop2("Duplicated allele label: ", alleles[dup])
  ## afreq
  afreq = attrs$afreq
  if (length(afreq) != length(alleles))
    stop2("Frequency vector doesn't match the number of alleles")
  if (round(sum(afreq), 3) != 1)
    stop2("Allele frequencies do not sum to 1 (after rounding to 3 decimal places): ", afreq)

  # name
  name = attrs$name
  if(length(name) != 1)
    stop2("Length of `name` must be 1: ", name)
  if (isTRUE(suppressWarnings(name == as.integer(name))))
    stop2("Attribute `name` cannot consist entirely of digits: ", name)

  # chrom
  chrom = attrs$chrom
  if(length(chrom) != 1)
    stop2("Length of `chrom` must be 1: ", chrom)

  # pedmembers and sex
  pedmembers = attrs$pedmembers
  if(length(pedmembers) != nrow(x))
    stop2("`pedmembers` attribute must have same length as nrows of the allele matrix")
  sex = attrs$sex
  if(length(sex) != nrow(x))
    stop2("`sex` attribute must have same length as nrows of the allele matrix")

  # mutation model
  mutmod = attrs$mutmod
  if(!is.null(mutmod)) {
    if (!requireNamespace("pedmut", quietly = TRUE))
      stop2("Package `pedmut` is not installed")

    pedmut::validateMutationModel(mutmod)
  }

  x
}

checkMutationMatrix = function(mutmat, alleles, identifier = NULL) {
  prefix = if(!is.null(identifier)) sprintf("%s mutation matrix: ", identifier) else ""
  N = length(alleles)

  if(!is.numeric(mutmat))
    stop2(prefix, "Type must be numeric, not ", typeof(mutmat))

  if (!identical(dm <- dim(mutmat), c(N,N)))
    stop2(prefix,
          sprintf("Dimensions (%d x %d) incompatible with number of alleles (%d)",
                  dm[1], dm[2], N))

  if(!identical(alleles, rownames(mutmat)) ||
     !identical(alleles, colnames(mutmat)))
    stop2(prefix, "Dimnames differ from allele names")

  if (any(round(rowSums(mutmat), 3) != 1))
    stop2(prefix, "Row sums are not 1")
}



