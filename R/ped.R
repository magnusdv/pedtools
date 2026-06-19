#' Pedigree construction
#'
#' This is the basic constructor of `ped` objects. Utility functions for
#' creating many common pedigree structures are described in [ped_basic]. See
#' also [as.ped()] and [readPed()], which are more liberal regarding the input
#' format.
#'
#' Each individual must have either both parents specified, or no parents.
#' Missing parents are indicated with entries "0", "" or NA in `fid` and `mid`.
#' Note that `id`,`fid`,`mid` are all converted to character vectors before
#' matching to establish the parent connections.
#'
#' If the pedigree is disconnected, it is split into its connected components
#' and returned as a list of `ped` objects.
#'
#' A singleton is a special `ped` object whose pedigree contains 1 individual.
#' The class attribute of a singleton is `c('singleton', 'ped')`.
#'
#' `singletons()` creates a list of singletons with the indicated labels and
#' sexes.
#'
#' Selfing, i.e. the presence of pedigree members whose father and mother are
#' the same individual, is allowed in `ped` objects. Any such "self-fertilizing"
#' parent must have undecided sex (`sex = 0`).
#'
#' @param id A character (or coercible to character) of individual ID labels.
#' @param fid,mid Vectors of the same length as `id`, naming each individual's
#'   father and mother. Missing parents (of founders) may be entered as "0", ""
#'   or NA.
#' @param sex A numeric of the same length as `id`, describing the genders of
#'   the individuals (in the same order as `id`.) Each entry must be either 1
#'   (=male), 2 (=female) or 0 (=unknown).
#' @param famid A character string. Default: An empty string.
#' @param reorder A logical indicating if the pedigree should be reordered so
#'   that all parents precede their children. Default: TRUE.
#' @param validate A logical indicating if a validation of the pedigree
#'   structure should be performed. Default: TRUE.
#' @param isConnected A logical indicating if the input is known to be a
#'   connected pedigree. Setting this to TRUE speeds up the processing. Default:
#'   FALSE.
#' @param detectLoops This argument is deprecated; loops are now always detected.
#' @param verbose A logical.
#'
#' @return A `ped` object, which is essentially a list with the following
#'   entries:
#'
#'   * `ID`: A character vector of ID labels. Unless the pedigree is reordered
#'   during creation, this equals `as.character(id)`
#'
#'   * `FIDX`: An integer vector with paternal indices: For each \eqn{j =
#'   1,2,...}, `FIDX[j]` is 0 if `ID[j]` has no father; otherwise `ID[FIDX[j]]`
#'   is the father of `ID[j]`.
#'
#'   * `MIDX`: An integer vector with maternal indices: For each \eqn{j =
#'   1,2,...}, `MIDX[j]` is 0 if `ID[j]` has no mother; otherwise `ID[MIDX[j]]`
#'   is the mother of `ID[j]`.
#'
#'   * `SEX`: An integer vector with gender codes. Unless the pedigree is
#'   reordered, this equals `as.integer(sex)`.
#'
#'   * `FAMID`: The family ID.
#'
#'   * `UNBROKEN_LOOPS`: A logical indicating if the pedigree has unbroken
#'   loops, or NA if the status is currently unknown.
#'
#'   * `LOOP_BREAKERS`: A matrix with loop breaker ID's in the first column and
#'   their duplicates in the second column. All entries refer to the internal
#'   IDs. This is usually set by [breakLoops()].
#'
#'   * `FOUNDER_INBREEDING`: A list of two potential entries, "autosomal" and
#'   "x"; both numeric vectors with the same length as `founders(x)`.
#'   `FOUNDER_INBREEDING` is always NULL when a new `ped` is created. See
#'   [founderInbreeding()].
#'
#'   * `MARKERS`: A list of `marker` objects, or NULL.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [newPed()], [ped_basic], [ped_modify], [ped_subgroups], [relabel()]
#'
#' @examples
#' # Trio
#' x = ped(id = 1:3, fid = c(0,0,1), mid = c(0,0,2), sex = c(1,2,1))
#'
#' # Female singleton
#' y = singleton('NN', sex = 2)
#'
#' # Selfing
#' z = ped(id = 1:2, fid = 0:1, mid = 0:1, sex = 0:1)
#' stopifnot(hasSelfing(z))
#'
#' # Disconnected pedigree: Trio + singleton
#' ped(id = 1:4, fid = c(2,0,0,0), mid = c(3,0,0,0), sex = c(1,1,2,1))
#'
#' # List of singletons
#' singletons(1:2)
#'
#' @export
ped = function(id, fid, mid, sex, famid = "", reorder = TRUE, validate = TRUE,
               detectLoops = TRUE, isConnected = FALSE, verbose = FALSE) {

  # Check input
  n = length(id)

  if(n == 0)
    stop2("`id` vector has length 0")
  if(length(fid) != n)
    stop2(sprintf("Incompatible input: length(id) = %d, but length(fid) = %d", n, length(fid)))
  if(length(mid) != n)
    stop2(sprintf("Incompatible input: length(id) = %d, but length(mid) = %d", n, length(mid)))
  if(length(sex) != n)
    stop2(sprintf("Incompatible input: length(id) = %d, but length(sex) = %d", n, length(sex)))

  # Coerce
  id = as.character(id)
  fid = as.character(fid)
  mid = as.character(mid)
  famid = as.character(famid)

  # Duplicated IDs
  if(anyDuplicated.default(id) > 0)
    stop2("Duplicated entry in `id` vector: ", id[duplicated(id)])

  # Parental index vectors (integer).
  missing = c("", "0", NA)
  FIDX = match(fid, id)
  FIDX[fid %in% missing] = 0L

  MIDX = match(mid, id)
  MIDX[mid %in% missing] = 0L

  if(any(is.na(FIDX)))
    stop2("`fid` entry does not appear in `id` vector: ", fid[is.na(FIDX)])
  if(any(is.na(MIDX)))
    stop2("`mid` entry does not appear in `id` vector: ", mid[is.na(MIDX)])

  if(all(FIDX + MIDX > 0))
    stop2("Pedigree has no founders")

  if(length(famid) != 1)
    stop2("`famid` must be a character string: ", famid)

  # Check for illegal entries in `sex``
  if(anyNA(match(sex, 0:2)))
    stop2("Illegal sex: ", .mysetdiff(sex, 0:2))
  sex = as.integer(sex)

  # Connected components
  if(!isConnected) {

    # Identify components
    comps = connectedComponents(id, fidx = FIDX, midx = MIDX)

    if(length(comps) > 1) {
      famids = paste0(famid, "_comp", seq_along(comps))

      pedlist = lapply(seq_along(comps), function(i) {
        idx = match(comps[[i]], id)
        ped(id = id[idx], fid = fid[idx], mid = mid[idx],
            sex = sex[idx], famid = famids[i], reorder = reorder,
            validate = validate, isConnected = TRUE, verbose = verbose)
      })

      return(structure(pedlist, names = famids, class = c("pedList", "list")))
    }
  }

  if(validate)
    validatePed(id = id, fidx = FIDX, midx = MIDX, sex = sex, famid = famid)

  # Initialise ped object
  x = newPed(id, FIDX, MIDX, sex, famid)

  # reorder so that parents precede their children
  if(reorder)
    x = parentsBeforeChildren(x)

  x
}


#' @export
#' @rdname ped
singleton = function(id = 1, sex = 1, famid = "") {
  if (length(id) != 1)
    stop2("`id` must have length 1")
  sex = validate_sex(sex, nInd = 1)
  newPed(ID = as.character(id), FIDX = 0L, MIDX = 0L, SEX = sex, FAMID = famid)
}


#' @export
#' @rdname ped
singletons = function(id, sex = 1) {
  n = length(id)
  id = as.character(id)
  sex = validate_sex(sex, nInd = n)

  lapply(seq_len(n), function(i)
    newPed(ID = id[i], FIDX = 0L, MIDX = 0L, SEX = sex[i], FAMID = ""))
}


#' Internal ped constructor
#'
#' This is the internal constructor of `ped` objects. It does not do any
#' validation of input other than simple type checking. In particular it should
#' only be used in programming scenarios where it is known that the input is a
#' valid, connected pedigree. End users are recommended to use the regular
#' constructor [ped()].
#'
#' See [ped()] for details about the input parameters.
#'
#' @param ID A character vector.
#' @param FIDX An integer vector.
#' @param MIDX An integer vector.
#' @param SEX An integer vector.
#' @param FAMID A string.
#' @param detectLoops (Deprecated).
#'
#' @return A `ped` object.
#'
#' @examples
#'
#' newPed("a", 0L, 0L, 1L, "")
#'
#' @export
newPed = function(ID, FIDX, MIDX, SEX, FAMID, detectLoops = NULL) {
  if(!all(is.character(ID), is.integer(FIDX), is.integer(MIDX),
          is.integer(SEX), is.character(FAMID)))
    stop2("Type error in the creation of `ped` object")

  LOOPS = hasLoop(fidx = FIDX, midx = MIDX)

  # Initialise ped object
  x = list(ID = ID,
           FIDX = FIDX,
           MIDX = MIDX,
           SEX = SEX,
           FAMID = FAMID,
           UNBROKEN_LOOPS = LOOPS,
           LOOP_BREAKERS = NULL,
           FOUNDER_INBREEDING = NULL,
           MARKERS = NULL)

  class(x) = if(length(ID) == 1) c("singleton", "ped") else "ped"

  x
}


#' Pedigree errors
#'
#' Validate the internal pedigree structure. The input may be either a (possibly
#' malformed) [ped()] object, or its defining vectors `id`, `fid`, `mid`, `sex`.
#'
#' @param x A `ped` object.
#' @param fidx,midx Integer vectors of parental indices (for internal use).
#' @inheritParams ped
#'
#' @return If no errors are detected, the function returns NULL invisibly.
#'   Otherwise, messages describing the errors are printed to the screen and an
#'   error is raised.
#'
#' @examples
#' x = nuclearPed()
#' validatePed(x)
#'
#' # Various errors
#' # validatePed(id = c(1,2), fid = c(2,0), mid = c(0,1), sex = c(1,2))
#'
#' @export
#'
validatePed = function(x = NULL, id = x$ID, fid = NULL, mid = NULL, sex = x$SEX,
                       famid = x$FAMID, fidx = x$FIDX, midx = x$MIDX) {
  id = as.character(id)

  if(!is.null(fid))
    fidx = match(fid, id, nomatch = 0L)
  if(!is.null(mid))
    midx = match(mid, id, nomatch = 0L);

  n = length(id)

  # Type verification (mainly for developer)
  stopifnot2(is.integer(fidx), is.integer(midx), is.numeric(sex))

  # Other verifications that don't need friendly messages at this point
  # (since they should be caught earlier during construction)
  stopifnot2(n > 0, length(fidx) == n, length(midx) == n, length(sex) == n,
             min(fidx) >= 0, min(midx) >= 0, max(fidx) <= n, max(midx) <= n)

  if(!is.null(famid))
    stopifnot2(is.character(famid), length(famid) == 1, !is.na(famid))

  errs = character(0)

  # Either 0 or 2 parents
  has1parent = (fidx > 0L) != (midx > 0L)
  if (any(has1parent))
    errs = c(errs, paste("Individual", id[has1parent], "has exactly 1 parent; this is not allowed"))

  # Sex
  if(anyNA(sex) || min(sex) < 0L || max(sex) > 2L)
    errs = c(errs, paste("Illegal sex:", unique(setdiff(sex, 0:2))))

  # Self ancestry
  self_anc = any_self_ancestry(id, fidx, midx)
  if(length(self_anc) > 0L)
    errs = c(errs, paste("Individual", self_anc, "is their own ancestor"))

  # Duplicated IDs
  if(anyDuplicated.default(id) > 0L)
    errs = c(errs, paste("Duplicated ID label:", id[duplicated.default(id)]))

  # Female fathers
  if(any(sex[fidx] == 2L, na.rm = TRUE)) {
    female_fathers_int = intersect(which(sex == 2L), fidx) # note: zeroes in FIDX disappear
    first_child = id[match(female_fathers_int, fidx)]
    errs = c(errs, paste("Individual", id[female_fathers_int],
                         "is female, but appear as the father of", first_child))
  }

  # Male mothers
  if(any(sex[midx] == 1L, na.rm = TRUE)) {
    male_mothers_int = intersect(which(sex == 1L), midx) # note: zeroes in MIDX disappear
    first_child = id[match(male_mothers_int, midx)]
    errs = c(errs, paste("Individual", id[male_mothers_int],
                         "is male, but appear as the mother of", first_child))
  }

  # Connected?
  #if (all(c(FIDX, MIDX) == 0))
  #    message("Pedigree is not connected.")

  if(length(errs) > 0) {
    errs = c("Malformed pedigree.", errs)
    stop2(paste0(errs, collapse = "\n "))
  }

  invisible(NULL)
}


any_self_ancestry = function(id, fidx, midx) {
  n = length(id)
  nseq = seq_len(n)

  # Quick check if anyone is their own parent
  self_parent = nseq == fidx | nseq == midx
  if(any(self_parent))
    return(id[self_parent])

  # OK means all known ancestors have been resolved.
  # Start with founders OK
  OK = fidx == 0L

  for(i in nseq) {

    # Ensure parent index 0 is treated as OK.
    parOK = c(TRUE, OK)
    OKnew = OK | parOK[fidx + 1L] & parOK[midx + 1L]

    # No change -> remaining FALSE are cyclic/unresolvable
    if(identical(OKnew, OK))
      break

    OK = OKnew
  }

  id[!OK]
}
