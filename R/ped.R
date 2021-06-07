#' Pedigree construction
#'
#' This is the basic constructor of `ped` objects. Utility functions for
#' creating many common pedigree structures are described in [ped_basic].
#'
#' A singleton is a special `ped` object whose pedigree contains 1 individual.
#' The class attribute of a singleton is `c('singleton', 'ped')`.
#'
#' Selfing, i.e. the presence of pedigree members whose father and mother are
#' the same individual, is allowed in `ped` objects. Any such "self-fertilizing"
#' parent must have undecided sex (`sex = 0`).
#'
#' If the pedigree is disconnected, it is split into its connected components
#' and returned as a list of `ped` objects.
#'
#' @param id a vector (numeric or character) of individual ID labels.
#' @param fid a vector of the same length as `id`, containing the labels of the
#'   fathers. In other words `fid[i]` is the father of `id[i]`, or 0 if `id[i]`
#'   is a founder.
#' @param mid a vector of the same length as `id`, containing the labels of the
#'   mothers. In other words `mid[i]` is the mother of `id[i]`, or 0 if `id[i]`
#'   is a founder.
#' @param sex a numeric of the same length as `id`, describing the genders of
#'   the individuals (in the same order as `id`.) Each entry must be either 1
#'   (=male), 2 (=female) or 0 (=unknown).
#' @param famid a character string. Default: An empty string.
#' @param reorder a logical. If TRUE, the pedigree is reordered so that all
#'   parents precede their children.
#' @param validate a logical. If TRUE, [validatePed()] is run before returning
#'   the pedigree.
#' @param isConnected a logical, by default FALSE. If it is known that the input
#'   is connected, setting this to TRUE speeds up the processing.
#' @param verbose a logical.
#'
#' @return A `ped` object, which is essentially a list with the following
#'   entries:
#'
#'   * `ID` : A character vector of ID labels. Unless the pedigree is reordered
#'   during creation, this equals `as.character(id)`
#'
#'   * `FIDX` : An integer vector with paternal indices: For each \eqn{j =
#'   1,2,...}, the entry `FIDX[j]` is 0 if `ID[j]` has no father within the
#'   pedigree; otherwise `ID[FIDX[j]]` is the father of `ID[j]`.
#'
#'   * `MIDX` : An integer vector with maternal indices: For each \eqn{j =
#'   1,2,...}, the entry `MIDX[j]` is 0 if `ID[j]` has no mother within the
#'   pedigree; otherwise `ID[MIDX[j]]` is the mother of `ID[j]`.
#'
#'   * `SEX` : An integer vector with gender codes. Unless the pedigree is
#'   reordered, this equals `as.integer(sex)`.
#'
#'   * `FAMID` : The family ID.
#'
#'   * `UNBROKEN_LOOPS` : A logical: TRUE if the pedigree is inbred.
#'
#'   * `LOOP_BREAKERS` : A matrix with loop breaker ID's in the first column and
#'   their duplicates in the second column. All entries refer to the internal
#'   IDs. This is usually set by [breakLoops()].
#'
#'   * `FOUNDER_INBREEDING` : A list of two potential entries, "autosomal" and
#'   "x"; both numeric vectors with the same length as `founders(x)`.
#'   `FOUNDER_INBREEDING` is always NULL when a new `ped` is created. See
#'   [founderInbreeding()].
#'
#'   * `MARKERS` : A list of `marker` objects, or NULL.
#'
#' @author Magnus Dehli Vigeland
#' @seealso [ped_basic], [ped_modify], [ped_subgroups], [relabel()]
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
#' w = ped(id = 1:4, fid = c(2,0,0,0), mid = c(3,0,0,0), sex = c(1,1,2,1))
#' stopifnot(is.pedList(w), length(w) == 2)
#'
#' @export
ped = function(id, fid, mid, sex, famid = "", reorder = TRUE, validate = TRUE, isConnected = FALSE, verbose = FALSE) {

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
  sex = as.integer(sex)
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

  # Initialise ped object
  x = newPed(id, FIDX, MIDX, sex, famid)

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
  ped(id = id, fid = 0, mid = 0, sex = sex, famid = famid)
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
#'
#' @return A `ped` object.
#'
#' @examples
#'
#' newPed("a", 0L, 0L, 1L, "")
#'
#' @export
newPed = function(ID, FIDX, MIDX, SEX, FAMID) {
  if(!all(is.character(ID), is.integer(FIDX), is.integer(MIDX),
          is.integer(SEX), is.character(FAMID)))
    stop2("Type error in the creation of `ped` object")

  # Initialise ped object
  x = list(ID = ID,
           FIDX = FIDX,
           MIDX = MIDX,
           SEX = SEX,
           FAMID = FAMID,
           UNBROKEN_LOOPS = FALSE,
           LOOP_BREAKERS = NULL,
           FOUNDER_INBREEDING = NULL,
           MARKERS = NULL)

  if(length(ID) == 1) {
    class(x) = c("singleton", "ped")
    return(x)
  }

  class(x) = "ped"

  # Detect loops (by trying to find a peeling order)
  nucs = peelingOrder(x)
  lastnuc_link = nucs[[length(nucs)]]$link
  x$UNBROKEN_LOOPS = is.null(lastnuc_link)

  x
}


#' Pedigree errors
#'
#' Validate the internal structure of a `ped` object.
#'
#' @param x object of class `ped`.
#'
#' @return If no errors are detected, the function returns NULL invisibly.
#'   Otherwise, messages describing the errors are printed to the screen and an
#'   error is raised.
#'
#' @export
validatePed = function(x) {
  ID = x$ID; FIDX = x$FIDX; MIDX = x$MIDX; SEX = x$SEX; FAMID = x$FAMID
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
  self_anc = any_self_ancestry(x)
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
  n = pedsize(x)
  nseq = 1:n
  FIDX = x$FIDX
  MIDX = x$MIDX


  # Quick check if anyone is their own parent
  self_parent = (nseq == FIDX) | (nseq == MIDX)
  if(any(self_parent))
    return(labels(x)[self_parent])

  fou_int = founders(x, internal = TRUE)
  OK = rep(FALSE, n)
  OK[fou_int] = TRUE

  # TODO: works, but not optimised for speed
  for(i in 1:n) {
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
  labels(x)[!OK]
}
