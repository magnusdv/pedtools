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
#' @param detectLoops A logical indicating if the presence of loops should be
#'   detected. Setting this to FALSE may speed up the processing of large
#'   pedigrees. Default: TRUE.
#' @param validate A logical indicating if a validation of the pedigree
#'   structure should be performed. Default: TRUE.
#' @param isConnected A logical indicating if the input is known to be a
#'   connected pedigree. Setting this to TRUE speeds up the processing. Default:
#'   FALSE.
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
  if(!all(sex %in% 0:2))
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
            validate = validate, detectLoops = detectLoops,
            isConnected = TRUE, verbose = verbose)
      })

      return(structure(pedlist, names = famids, class = c("pedList", "list")))
    }
  }

  # Initialise ped object
  x = newPed(id, FIDX, MIDX, sex, famid, detectLoops = FALSE) # TODO

  # Detect loops (by trying to find a peeling order)
  if(detectLoops)
    x$UNBROKEN_LOOPS = hasUnbrokenLoops(x)

  if(validate)
    validatePed(x)

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
  newPed(ID = as.character(id), FIDX = 0L, MIDX = 0L, SEX = sex,
         FAMID = famid, detectLoops = FALSE)
}


#' @export
#' @rdname ped
singletons = function(id, sex = 1) {
  n = length(id)
  id = as.character(id)
  sex = validate_sex(sex, nInd = n)

  lapply(seq_len(n), function(i)
    newPed(ID = id[i], FIDX = 0L, MIDX = 0L, SEX = sex[i],
           FAMID = "", detectLoops = FALSE))
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
#' @param detectLoops A logical.
#'
#' @return A `ped` object.
#'
#' @examples
#'
#' newPed("a", 0L, 0L, 1L, "")
#'
#' @export
newPed = function(ID, FIDX, MIDX, SEX, FAMID, detectLoops = TRUE) {
  if(!all(is.character(ID), is.integer(FIDX), is.integer(MIDX),
          is.integer(SEX), is.character(FAMID)))
    stop2("Type error in the creation of `ped` object")

  # Initialise ped object
  x = list(ID = ID,
           FIDX = FIDX,
           MIDX = MIDX,
           SEX = SEX,
           FAMID = FAMID,
           UNBROKEN_LOOPS = NA,
           LOOP_BREAKERS = NULL,
           FOUNDER_INBREEDING = NULL,
           MARKERS = NULL)

  if(length(ID) == 1) {
    class(x) = c("singleton", "ped")
    x$UNBROKEN_LOOPS = FALSE
    return(x)
  }

  class(x) = "ped"
  if(detectLoops)
    x$UNBROKEN_LOOPS = hasUnbrokenLoops(x)

  x
}


#' Pedigree errors
#'
#' Validate the internal pedigree structure. The input may be either a (possibly
#' malformed) [ped()] object, or its defining vectors `id`, `fid`, `mid`, `sex`.
#'
#' @param x A `ped` object.
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
validatePed = function(x = NULL, id = NULL, fid = NULL, mid = NULL, sex = NULL) {
  if(!is.null(x)) {
    ID = x$ID; FIDX = x$FIDX; MIDX = x$MIDX; SEX = x$SEX; FAMID = x$FAMID
  }
  else {
    ID = as.character(id); FIDX = match(fid, id, nomatch = 0L); MIDX = match(mid, id, nomatch = 0L);
    SEX = as.integer(sex); FAMID = ""
  }

  n = length(ID)

  # Type verification (mainly for developer)
  stopifnot2(is.character(ID), is.integer(FIDX), is.integer(MIDX), is.integer(SEX),
            is.character(FAMID), is.singleton(x) == (n == 1))

  # Other verifications that don't need friendly messages at this point
  # (since they should be caught earlier during construction)
  stopifnot2(n > 0, length(FIDX) == n, length(MIDX) == n, length(SEX) == n,
            all(FIDX >= 0), all(MIDX >= 0), all(FIDX <= n), all(MIDX <= n),
            length(FAMID) == 1)

  errs = character(0)

  # Either 0 or 2 parents
  has1parent = (FIDX > 0) != (MIDX > 0)
  if (any(has1parent))
    errs = c(errs, paste("Individual", ID[has1parent], "has exactly 1 parent; this is not allowed"))

  # Sex
  if (!all(SEX %in% 0:2))
    errs = c(errs, paste("Illegal sex:", unique(setdiff(SEX, 0:2))))

  # Self ancestry
  self_anc = any_self_ancestry(list(ID = ID, FIDX = FIDX, MIDX = MIDX))
  if(length(self_anc) > 0)
    errs = c(errs, paste("Individual", self_anc, "is their own ancestor"))

  # If singleton: return here
  # if(n == 1) return()

  # Duplicated IDs
  if(anyDuplicated.default(ID) > 0)
    errs = c(errs, paste("Duplicated ID label:", ID[duplicated(ID)]))

  # Female fathers
  if(any(SEX[FIDX] == 2)) {
    female_fathers_int = intersect(which(SEX == 2), FIDX) # note: zeroes in FIDX disappear
    first_child = ID[match(female_fathers_int, FIDX)]
    errs = c(errs, paste("Individual", ID[female_fathers_int],
                         "is female, but appear as the father of", first_child))
  }

  # Male mothers
  if(any(SEX[MIDX] == 1)) {
    male_mothers_int = intersect(which(SEX == 1), MIDX) # note: zeroes in MIDX disappear
    first_child = ID[match(male_mothers_int, MIDX)]
    errs = c(errs, paste("Individual", ID[male_mothers_int],
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


any_self_ancestry = function(x) {
  ID = x$ID
  FIDX = x$FIDX
  MIDX = x$MIDX

  n = length(ID)
  nseq = seq_len(n)

  # Quick check if anyone is their own parent
  self_parent = (nseq == FIDX) | (nseq == MIDX)
  if(any(self_parent))
    return(ID[self_parent])

  fou_int = which(FIDX == 0)
  OK = rep(FALSE, n)
  OK[fou_int] = TRUE

  # TODO: works, but not optimised for speed
  for(i in nseq) { # note that i is not used
    parents = which(OK)
    children = which(FIDX %in% parents | MIDX %in% parents)

    fatherOK = OK[FIDX[children]]
    motherOK = OK[MIDX[children]]
    childrenOK = children[fatherOK & motherOK]

    # If these were already ok, there is nothing more to do
    if(all(OK[childrenOK]))
      break

    OK[childrenOK] = TRUE
  }
  ID[!OK]
}
