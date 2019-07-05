#' Attach genotype data to pedigrees
#'
#' Reads genotype data from a reference dataframe consisting of four columns:
#' Id, Marker, Allele 1 and Allele 2. Creates one marker per row and attaches
#' them to the pedigree in the correct format. Allele denomination and frequency
#' data are read from the frequency dataframe which must contain a table with
#' one allele per row and one marker per column, where a cell (allele, marker)
#' contains the frequency of that allele in the given population. Only positive
#' numerical values are passed to the pedigree (alleles with frequencies 0 or NA
#' are ignored). Optionally, the user can choose to only attach some of the
#' markers in the reference file by passing a list of which ones to attach.
#'
#' If the pedigree has genotype information already attached to it, that info is
#' left intact.
#'
#' @param ped a `ped` object
#' @param references a [data.frame()] with 4 columns (Id, Marker, Allele 1,
#'   Allele 2)
#' @param frequencies a [data.frame()] with one column per marker and one row
#'   per allele where a cell (allele, marker) contains a value for the frequency
#' @param markers a vector of strings with the marker names that should be
#'   attached to the pedigree. Defaults to `NULL` in which case all information
#'   available in `references` will be attached to the pedigree
#' @return a `ped` object with the attached markers
#'
#' @author Elias Hernandis <eliashernandis@gmail.com>
#' @seealso [as.ped()]
#' @export
attachGenotypeToPed = function(ped, references, frequencies, markers = NULL) {
  # if no markers were specified, load all markers
  if (is.null(markers)) {
    markers = unique(as.vector(references[,2]))
  }

  for (markerName in markers) {
    # get allele denominations
    als = as.vector(frequencies[!is.na(frequencies[markerName]),1], mode="character")

    # get frequency data
    freqs = frequencies[!is.na(frequencies[markerName]),markerName]

    # get assignment of alleles to persons in the pedigree
    ret = apply(ref[ref[,2] == markerName,], 1, function(row) {
      as.vector(row[c(3,4)], mode="character")
    })
    args = split(ret, rep(1:ncol(ret)))
    args = setNames(args, ref[ref[,2] == markerName,1])
    args$x = ped
    args$alleles = als
    args$afreq = freqs
    args$name = markerName
    print(args)

    # build marker
    m = do.call(marker, args)

    # attach marker to pedigree
    ped = addMarkers(ped, m)
  }

  ped
}
