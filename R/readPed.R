#' Read a pedigree from file
#'
#' @param pedfile A file name
#' @param header A logical. If NA, the program will interpret the first line as
#'   a header line if the first entry contains "id" AND the word "sex" is an
#'   entry.
#' @param ... Further parameters passed on to [read.table()], e.g. `sep`,
#'   `comment.char` and `quote`.
#'
#' @return A [ped] object or a list of such.
#'
#' @examples
#'
#' tf = tempfile()
#'
#' ### Write and read a trio
#' trio = data.frame(id = 1:3, fid = c(0,0,1), mid = c(0,0,2), sex = c(1,2,1))
#' write.table(trio, file = tf, row.names = FALSE)
#' readPed(tf)
#'
#' # With marker data in one column
#' trio.marker = cbind(trio, M = c("1/1", "2/2", "1/2"))
#' write.table(trio.marker, file = tf, row.names = FALSE)
#' readPed(tf, sep = "/")
#'
#' # With marker data in two allele columns
#' trio.marker2 = cbind(trio, A1 = c(1,2,1), A2 = c(1,2,2))
#' write.table(trio.marker2, file = tf, row.names = FALSE)
#' readPed(tf)
#'
#' ### Two singletons in the same file
#' singles = data.frame(id = c("S1", "S2"),
#'                      fid = c(0,0), mid = c(0,0), sex = c(2,1),
#'                      M = c("9/14.2", "9/9"))
#' write.table(singles, file = tf, row.names = FALSE)
#' readPed(tf, sep = "/")
#'
#' ### Two trios in the same file
#' trio2 = cbind(famid = rep(c("trio1", "trio2"), each = 3), rbind(trio, trio))
#'
#' # Without column names
#' write.table(trio2, file = tf, row.names = FALSE)
#' readPed(tf)
#'
#' # With column names
#' write.table(trio2, file = tf, col.names = FALSE, row.names = FALSE)
#' readPed(tf, famid = 1, id = 2, fid = 3, mid = 4, sex = 5)
#'
#' # Cleanup
#' unlink(tf)
#'
#' @inheritParams as.ped.data.frame
#' @importFrom utils read.table
#' @export
readPed = function(pedfile, header = NA, famid_col = NA, id_col = NA, fid_col = NA,
                   mid_col = NA, sex_col = NA, marker_col = NA,
                   locusAttributes = NULL, missing = 0,
                   sep = NULL, validate = TRUE, ...) {

  # Check for deprecated arguments
  dots = list(...)
  for(arg in names(dots)) {
    if(!is.na(pmatch(arg, "locus_annotations"))) {
      warning("Argument `locus_annotations` is deprecated; use `locusAttributes` instead")
      locusAttributes = dots[[arg]]
    }
    if(!is.na(pmatch(arg, "allele_sep"))) {
      warning("Argument `allele_sep` is deprecated; use `sep` instead")
      sep = dots[[arg]]
    }
  }

  # If header = NA, check first line
  if(is.na(header)) {
    # Read first line
    first = tolower(scan(pedfile, what = "", nlines = 1, quiet = TRUE, ...))

    # Interpret as header line if 1) first element contains "id" and 2) "sex" is an entry
    header = grepl("id", first[1], fixed = TRUE) && "sex" %in% first
  }

  ped.df = read.table(pedfile, header = header, colClasses = "character", check.names = FALSE, ...)

  as.ped(ped.df, famid_col = famid_col, id_col = id_col, fid_col = fid_col,
         mid_col = mid_col, sex_col = sex_col, marker_col = marker_col,
         locusAttributes = locusAttributes, missing = missing,
         sep = sep, validate = validate)

}
