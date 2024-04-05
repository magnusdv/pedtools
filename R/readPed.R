#' Read a pedigree from file
#'
#' Reads a text file in pedigree format, or something fairly close to it.
#'
#' If there are no headers, and no column information is provided by the user,
#' the program assumes the following column order:
#'
#' * family ID (optional; guessed from the data)
#'
#' * individual ID
#'
#' * father's ID
#'
#' * mother's ID
#'
#' * sex
#'
#' * marker data (remaining columns)
#'
#' #### Reading SNP data
#'
#' Adding the argument `locusAttributes = "snp-AB"`, sets all markers to be
#' equifrequent SNPs with alleles A and B. Moreover, the letters A and B may be
#' replaced by other single-character letters or numbers, e.g., "snp-12" gives
#' alleles 1 and 2.
#'
#' @inheritParams as.ped.data.frame
#' @param pedfile A file name
#' @param colSep A column separator character, passed on as the `sep` argument
#'   of [read.table()]. The default is to separate on white space, that is, one
#'   or more spaces, tabs, newlines or carriage returns. (Note: the parameter
#'   `sep` is used to indicate allele separation in genotypes.)
#' @param header A logical. If NA, the program will interpret the first line as
#'   a header line it contains both "id" and "sex" as part of some entries
#'   (ignoring case).
#' @param colSkip Columns to skip, given as a vector of indices or columns
#'   names. If given, these columns are removed directly after `read.table()`,
#'   before any other processing.
#' @param ... Further parameters passed on to [read.table()], e.g.
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
#' readPed(tf)
#'
#' # With marker data in two allele columns
#' trio.marker2 = cbind(trio, M.1 = c(1,2,1), M.2 = c(1,2,2))
#' write.table(trio.marker2, file = tf, row.names = FALSE)
#' readPed(tf)
#'
#' ### Two singletons in the same file
#' singles = data.frame(id = c("S1", "S2"),
#'                      fid = c(0,0), mid = c(0,0), sex = c(2,1),
#'                      M = c("9/14.2", "9/9"))
#' write.table(singles, file = tf, row.names = FALSE)
#' readPed(tf)
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
#' ### With non-standard `sex` codes
#' trio3 = data.frame(id = 1:3, fid = c(0,0,1), mid = c(0,0,2),
#'                    sex = c("male","female","?"))
#' write.table(trio3, file = tf, row.names = FALSE)
#' readPed(tf, sexCodes = list(male = "male", female = "female", unknown = "?"))
#'
#' # Cleanup
#' unlink(tf)
#'
#' @importFrom utils read.table
#' @export
readPed = function(pedfile, colSep = "", header = NA,
                   famid_col = NA, id_col = NA, fid_col = NA,
                   mid_col = NA, sex_col = NA, marker_col = NA,
                   locusAttributes = NULL, missing = 0,
                   sep = NULL, colSkip = NULL, sexCodes = NULL,
                   addMissingFounders = FALSE, validate = TRUE, ...) {

  # If header = NA, check first line
  if(is.na(header)) {
    # Read first line
    first = tolower(scan(pedfile, what = "", nlines = 1, quiet = TRUE, ...))

    # Interpret as header line if it contains both "id" and "sex"
    header = any(grepl("id", first, fixed = TRUE)) && any(grepl("sex", first, fixed = TRUE))
  }

  ped.df = read.table(pedfile, sep = colSep, header = header,
                      colClasses = "character", check.names = FALSE, ...)
  if(!is.null(colSkip)) {
    if(is.character(colSkip))
      colSkip = match(colSkip, names(ped.df), nomatch = 0)
    if(is.numeric(colSkip) && any(colSkip > 0))
      ped.df = ped.df[, -colSkip[colSkip > 0], drop = FALSE]
  }

  nc = ncol(ped.df)

  # TODO: Clean this up and make smarter
  # guess columns if no header info
  if(!header && isTRUE(all(is.na(c(famid_col, id_col, fid_col, mid_col, sex_col, marker_col))))) {
    hasFamid = anyDuplicated(ped.df[,1]) || (nc >=4 && isTRUE(all(0 == ped.df[,3:4])))
    if(hasFamid) {
      famid_col = 1
      id_col = 2
      fid_col = 3
      mid_col = 4
      sex_col = if(nc > 4) 5 else NA
    }
    else {
      id_col = 1
      fid_col = 2
      mid_col = 3
      sex_col = if(nc > 3) 4 else NA
    }
  }

  as.ped(ped.df, famid_col = famid_col, id_col = id_col, fid_col = fid_col,
         mid_col = mid_col, sex_col = sex_col, marker_col = marker_col,
         locusAttributes = locusAttributes, missing = missing,
         sep = sep, sexCodes = sexCodes, addMissingFounders = addMissingFounders,
         validate = validate)

}
