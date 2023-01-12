#' Marker objects
#'
#' Creating a marker object associated with a pedigree. The function `marker()`
#' returns a marker object, while `addMarker()` first creates the marker and
#' then attaches it to `x`.
#'
#' @param x A `ped` object.
#' @param ... One or more expressions of the form `id = genotype`, where `id` is
#'   the ID label of a member of `x`, and `genotype` is a numeric or character
#'   vector of length 1 or 2 (see Examples).
#' @param geno A character vector of length `pedsize(x)`, with genotypes written
#'   in the format "a/b".
#' @param allelematrix A matrix with 2 columns and `pedsize(x)` rows. If this is
#'   non-NULL, then `...` must be empty.
#' @param alleles A character containing allele names. If not given, and `afreq`
#'   is named, `names(afreq)` is used. The default action is to take the sorted
#'   vector of distinct alleles occurring in `allelematrix`, `geno` or `...`.
#' @param afreq A numeric of the same length as `alleles`, indicating the
#'   population frequency of each allele. A warning is issued if the frequencies
#'   don't sum to 1 after rounding to 3 decimals. If the vector is named, and
#'   `alleles` is not NULL, an error is raised if `setequal(names(afreq),
#'   alleles)` is not TRUE. If `afreq` is not specified, all alleles are given
#'   equal frequencies.
#' @param chrom A single integer: the chromosome number. Default: NA.
#' @param posMb A nonnegative real number: the physical position of the marker,
#'   in megabases. Default: NA.
#' @param name A character string: the name of the marker. Default: NA.
#' @param mutmod,rate Mutation model parameters to be passed on to
#'   [pedmut::mutationModel()]; see there for details. Note: `mutmod`
#'   corresponds to the `model` parameter. Default: NULL (no mutation model).
#' @param locusAttr A list with names `alleles`, `afreq`, `chrom`, `name`,
#'   `posMb`, `mutmod`, `rate` (or a subset of these). This can be used as an
#'    alternative to entering the arguments as function parameters.
#' @param NAstrings A character vector containing strings to be treated as
#'   missing alleles. Default: `c("", "0", NA, "-")`.
#' @param validate A logical indicating if the validity of the marker object
#'   should be checked. Default: TRUE.
#' @param validateMut A logical indicating if the mutation model (if present)
#'   should be checked.
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
#' @seealso Get/set marker attributes: [marker_getattr], [marker_setattr].
#'
#' Retrieve various marker properties: [marker_prop], [nMarkers()],
#'
#' Add alleles to an existing marker: [addAllele()]
#'
#' Attach multiple markers: [marker_attach]
#'
#'
#'
#' @examples
#' x = nuclearPed(father = "fa", mother = "mo", children = "child")
#'
#' # An empty SNP with alleles "A" and "B"
#' marker(x, alleles = c("A", "B"))
#'
#' # Creating and attaching to `x`
#' addMarker(x, alleles = c("A", "B"))
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
                  mutmod = NULL, rate = NULL,
                  NAstrings = c(0, "", NA, "-"),
                  validate = TRUE, validateMut = validate) {

  # Some parameters cannot have length 0 or be ""
  if(length(chrom) == 0) chrom = NA
  if(length(posMb) == 0) posMb = NA
  if(length(name) == 0 || identical(name, "")) name = NA

  labs = labels(x)

  # Parse genotypes into list of 2-vectors
  glist = parseGeno(geno) %||% parseDots(...)

  if(length(glist) && !is.null(allelematrix))
    stop2("When specifying genotypes, `allelematrix` must be NULL")

  # Allele matrix (N x 2) of actual alleles
  m = allelematrix %||% glist2amat(glist, labs)

  # Trim whitespace of alleles
  # Assumption to avoid trimws (slow!): At most 1 space, no linebreaks
  m[] = sub(" ", "", m, fixed = TRUE)

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
    # Create model (postpone validation)
    mutmod = pedmut::mutationModel(mutmod, alleles = ALS, afreq = AFR, rate = rate, validate = FALSE)
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

# Convert genotypes given in marker(x, ...) into named list
parseDots = function(...) {
  dots = eval(substitute(alist(...)))
  if(!length(dots))
    return(NULL)

  glist = lapply(dots, eval.parent)

  lg = lengths(glist)

  # Split compound genotypes, e.g., "a/b"
  g1 = as.character(unlist(glist[lg == 1]))
  glist[lg == 1] = strsplit(g1, split = "/", fixed = TRUE)

  glist
}

# Convert `geno` argument into named list
parseGeno = function(geno) {
  if(is.null(geno))
    return(NULL)

  storage.mode(geno) = "character"
  strsplit(geno, split = "/", fixed = TRUE)
}

# Convert glist to allele matrix (actual alleles)
glist2amat = function(glist, labs) {
  # Empty allele matrix
  m = matrix(0, ncol = 2, nrow = length(labs))

  if(!length(glist))
    return(m)

  if(is.null(names(glist)) && length(glist) == length(labs))
    names(glist) = labs

  nms = names(glist)
  if(is.null(nms) || any(nms == ""))
    stop2("Genotypes must be named")

  unkn = setdiff(nms, unlist(labs))
  if(length(unkn))
    stop2("Unknown ID label: ", unkn)

  lg = lengths(glist)

  if(any(bad <- (lg < 1 | lg > 2))) {
    i = which(bad)[1]
    stop2(sprintf("Unclear genotype for `%s`: ", nms[i]), glist[[i]])
  }

  # Recycle single alleles
  glist[lg == 1] = lapply(glist[lg == 1], rep, length.out = 2)

  if(d <- anyDuplicated(nms))
    stop2(sprintf("Multiple genotypes given for individual `%s`", nms[d]))

  # Fill matrix
  idx = match(nms, labs)
  m[idx, ] = do.call(rbind, glist)

  m
}

#' @rdname marker
#' @export
addMarker = function(x, ..., geno = NULL, allelematrix = NULL, alleles = NULL,
                     afreq = NULL, chrom = NA, posMb = NA, name = NA,
                     mutmod = NULL, rate = NULL, locusAttr = NULL,
                     NAstrings = c(0, "", NA, "-"), validate = TRUE) {

  if(is.pedList(x)) {
    if(!is.null(allelematrix))
      stop2("The argument `allelematrix` cannot be used when `x` is a list of pedigrees")

    # If attributes given as list, use these
    if(!is.null(locusAttr)) {
      locusAttr = checkLocusAttribs(locusAttr)[[1]] # NB: returns list of lists
      alleles = locusAttr$alleles %||% alleles
      afreq = locusAttr$afreq %||% afreq
      chrom = locusAttr$chrom %||% chrom
      name = locusAttr$name %||% name
      posMb = locusAttr$posMb %||% posMb
      mutmod = locusAttr$mutmod %||% mutmod
      rate = locusAttr$rate %||% rate
    }

    if(is.null(alleles) && is.null(afreq))
      stop2("Either `alleles` or `afreq` must be specified when `x` is a list of pedigrees")

    glist = parseGeno(geno) %||% parseDots(...)
    if(hasGeno <- !is.null(glist)) {
      if(is.null(nms <- names(glist)))
        stop2("Genotypes must be named when `x` is a ped list")

      unkn = setdiff(nms, unlist(labels(x)))
      if(length(unkn))
        stop2("Unknown ID label: ", unkn)
    }

    y = lapply(x, function(comp) {
      labsi = labels(comp)
      mi = if(hasGeno) glist2amat(glist[intersect(nms, labsi)], labsi) else NULL
      addMarker(comp, allelematrix = mi, alleles = alleles,
                afreq = afreq, chrom = chrom, posMb = posMb, name = name,
                mutmod = mutmod, rate = rate, locusAttr = locusAttr,
                NAstrings = NAstrings, validate = validate)
    })
    return(y)
  }

  if(!is.ped(x))
    stop2("Input to `addMarker()` must be a `ped` object or a list of such")

  # If attributes given as list, use these
  if(!is.null(locusAttr)) {
    locusAttr = checkLocusAttribs(locusAttr)
    alleles = locusAttr$alleles %||% alleles
    afreq = locusAttr$afreq %||% afreq
    chrom = locusAttr$chrom %||% chrom
    name = locusAttr$name %||% name
    posMb = locusAttr$posMb %||% posMb
    mutmod = locusAttr$mutmod %||% mutmod
    rate = locusAttr$rate %||% rate
  }

  m = marker(x, ..., geno = geno, allelematrix = allelematrix,
             alleles = alleles, afreq = afreq,
             chrom = chrom, posMb = posMb, name = name,
             mutmod = mutmod, rate = rate, NAstrings = NAstrings,
             validate = validate)
  addMarkers(x, m)
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

validateMarker = function(x, validateMut = TRUE) {
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
  if (all(strsplit(name, "", fixed = TRUE)[[1]] %in% 0:9))
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
  if(validateMut && !is.null(mutmod)) {
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



