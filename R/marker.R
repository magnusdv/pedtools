#' Marker objects
#'
#' Creating a marker object associated with a pedigree
#'
#'
#' @param x a [`ped`] object
#' @param ... one or more expressions of the form `id = genotype`, where `id` is
#'   the ID label of a member of `x`, and genotype is a numeric or character
#'   vector of length 2 (see Examples).
#' @param alleles a character (or coercible to character) containing allele
#'   names. If not given, the default is to take the sorted vector of distinct
#'   alleles occuring in `...`.
#' @param afreq a numeric of the same length as `alleles`, indicating the
#'   population frequency of each allele. A warning is issued if the frequencies
#'   don't sum to 1 after rounding to 3 decimals. If `afreq` is not specified,
#'   all alleles are given equal frequencies.
#' @param chrom a single integer: the chromosome number. Default: NA.
#' @param posMb a single numeric: the phsyical position of the marker, in
#'   megabases. Default: NA.
#' @param posCm a single numeric: the centiMorgan position of the marker.
#'   Default: NA.
#' @param name a single character: the name of the marker. Default: NA.
#' @param mutmat a mutation matrix, or a list of two such matrices named
#'   'female' and 'male'. If given, the mutation matrices must be square, with
#'   the allele labels as dimnames, and each row must sum to 1 (after rounding
#'   to 3 decimals). Default: NULL.
#'
#' @return An object of class `marker`: This is a numerical 2-column matrix with
#'   one row per individual, and attributes 'alleles' (a character vector with
#'   allele labels), 'afreq` (allele frequencies), 'chrom' (chromosome number),
#'   'posMb' (physical location in megabases),'posCm' (position in centiMorgan),
#'   'name' (marker identifier) and 'mutmat' (a list of two (male and female)
#'   mutation matrices).
#'
#' @author Magnus Dehli Vigeland
#'
#' @examples
#' x = nuclearPed(father="fa", mother="mo", children="child")
#'
#' # Create a SNP marker with alleles 1 and 2.
#' # Father is homozygous 1/1; mother is 2/2; child is heterozygous.
#' marker(x, fa=1, mo=2, child=1:2)
#'
#' @export
marker = function(x, ...,  alleles = NULL, afreq = NULL, chrom = NA,
                  posMb = NA, posCm = NA, name = NA, mutmat = NULL) {

  # Initalize empty allele matrix
  m = matrix(0, ncol = 2, nrow = pedsize(x))

  # Capture genotypes given in dots
  dots = eval(substitute(alist(...)))

  ids_int = internalID(x, names(dots))
  genos = lapply(dots, eval.parent)

  for(i in seq_along(dots)) {
    g = genos[[i]]
    if(!is.vector(g) || !length(g) %in% 1:2) # TODO: deal with NA and '' inputs
      stop("Genotype is not a vector of length 1 or 2: ", deparse(g), call.=FALSE)
    m[ids_int[i], ] = g
  }

  .createMarkerObject(m, alleles = alleles, afreq = afreq, chrom = chrom,
      posMb = posMb, posCm=posCm, name = name, mutmat = mutmat,
      pedmembers = x$LABELS, sex = x$SEX)
}


as.marker.matrix = function(allelematrix, x, alleles = NULL, afreq = NULL, chrom = NA,
                            posMb = NA, posCm = NA, name = NA, mutmat = NULL)
NULL



.createMarkerObject = function(matr, alleles = NULL, afreq = NULL, chrom = NA,
                               posMb = NA, posCm = NA, name = NA, mutmat = NULL,
                               pedmembers = NULL, sex = NULL) {
  if(!is.null(alleles) && !all(matr %in% c(0, alleles))) {
    notfound = setdiff(matr, c(0, alleles))
    stop("Alleles used in genotypes but not included in `alleles` argument: ",
         paste(notfound, collapse=", "), call.=FALSE)
  }
  if(!is.null(afreq)) {
    if(is.null(alleles))
       stop("Argument `alleles` cannot be NULL if `afreq` is non-NULL", call.=FALSE)
    if (length(afreq) != length(alleles))
      stop("Number of alleles doesn't match length of frequency vector", call.=FALSE)
    if (round(sum(afreq), 3) != 1)
      stop("Allele frequencies do not sum to 1: ", paste(afreq, collapse = ", "), call.=FALSE)
  }
  if(length(pedmembers) != nrow(matr))
    stop("Number of allele matrix must equal the length of `pedmembers`", call.=FALSE)
  if(length(sex) != nrow(matr))
    stop("Number of allele matrix must equal the length of `sex`", call.=FALSE)

  ### Alleles
  if (is.null(alleles)) {
    vec = as.vector(matr)
    alleles = unique.default(vec[vec != 0])
    if (!length(alleles)) alleles = 1
  }
  else if (!anyNA(suppressWarnings(as.numeric(alleles))))
    alleles = as.numeric(alleles) # ensure numerical sorting if appropriate

  # Sort (same order used below for frequencies)
  allele_order = order(alleles)
  alleles = alleles[allele_order]

  ### Frequencies
  if (is.null(afreq)) {
    nalleles = length(alleles)
    afreq = rep.int(1, nalleles)/nalleles
  }
  else
    afreq = afreq[allele_order]

  ### Mutation matrices
  if (!is.null(mutmat)) {
    stopifnot(is.list(mutmat) || is.matrix(mutmat))
    # If single matrix given: make sex specific list
    if (is.matrix(mutmat)) {
      mutmat = .checkMutationMatrix(mutmat, alleles)
      mutmat = list(male = mutmat, female = mutmat)
    } else {
      stopifnot(length(mutmat) == 2, setequal(names(mutmat), c("female", "male")))
      mutmat$female = .checkMutationMatrix(mutmat$female, alleles, label = "female")
      mutmat$male = .checkMutationMatrix(mutmat$male, alleles, label = "male")
    }
  }

  ### Create the internal allele matrix
  m = match(matr, alleles, nomatch = 0)

  ### Add everything else as attributes, including class = "marker".
  attributes(m) = list(dim = dim(matr), alleles = as.character(alleles),
                           afreq = afreq, chrom = chrom, posMb = posMb,
                           posCm=posCm, name = name, mutmat = mutmat,
                           pedmembers = pedmembers, sex = sex, class = "marker")
  m
}

.checkMutationMatrix = function(mutmat, alleles, label = "") {
    # Check that mutation matrix is compatible with allele number / names.  Sort matrix
    # according to the given allele order (this is important since dimnames are NOT used in
    # calculations).
    N = length(alleles)
    if (label != "")
        label = sprintf("%s ", label)
    if (any((dm <- dim(mutmat)) != N))
        stop(sprintf("Dimension of %smutation matrix (%d x %d) inconsistent with number of alleles (%d).",
            label, dm[1], dm[2], N))
    if (any(round(.rowSums(mutmat, N, N), 3) != 1))
        stop(sprintf("Row sums of %smutation matrix are not 1.", label))
    alleles.char = as.character(alleles)
    if (!setequal(rownames(mutmat), alleles) || !setequal(colnames(mutmat), alleles))
        stop(sprintf("Dimnames of %smutation do not match allele names.", label))
    m = mutmat[alleles.char, alleles.char]

    # lumbability: always lumpable (any alleles) if rows are identical (except diagonal)
    attr(m, "lumpability") = NA
    if (N > 2) {
        if (all(vapply(1:N, function(i) diff(range(m[-i, i])) == 0, logical(1))))
            attr(m, "lumpability") = "always"  # If each column has 0 range
    }
    m
}


.setSNPfreqs = function(x, newfreqs) {
    stopifnot(all(vapply(x$markerdata, function(m) attr(m, "nalleles"), numeric(1)) == 2))
    newfreqs = rep(newfreqs, length = x$nMark)
    for (i in seq_len(x$nMark)) attr(x$markerdata[[i]], "afreq") = c(newfreqs[i], 1 - newfreqs[i])
    x
}
