#' Pedigree construction
#'
#' Basic construction of `ped` objects. Utility functions for creating many
#' common pedigree structures are described in [ped_basic].
#'
#' Internally, this happens: ... #TODO
#'
#' A singleton is a special `ped` object whose pedigree contains 1
#' individual. The class attribute of a singleton is `c('singleton',
#' 'ped')`
#'
#' @param id a vector (numeric or character) of individual ID labels.
#' @param fid a vector of the same length as `id`, containing the labels of the fathers.
#'   In other words `fid[i]` is the father of `id[i]`, or 0 if `id[i]` is a founder.
#' @param mid a vector of the same length as `id`, containing the labels of the mothers.
#'   In other words `mid[i]` is the mother of `id[i]`, or 0 if `id[i]` is a founder.
#' @param sex a numeric of the same length as `id`, describing the genders of the individuals
#'   (in the same order as `id`.) Each entry must be either 1 (=male), 2 (=female) or 0 (=unknown).
#' @param famid a character string. Default: An empty string.
#' @param reorder a logical. If TRUE, the pedigree is reordered so that all
#'   parents precede their children.
#' @param check a logical. If TRUE, [checkped()] is run on the
#'   pedigree before it is returned.
#' @param verbose a logical.
#' @param ... further arguments
#'
#' @return A `ped` object, which is essentially a list with the following
#'   entries:
#'
#'   * `ID` : A character vector of ID labels. Unless the pedigree is reordered during creation, this equals `as.character(id)`
#'   * `FIDX` : An integer vector with paternal indices: For each $j=1,2,...$, `ID[FIDX[j]]` is the father of `ID[j]`, or 0 if `ID[j]` has no father within the pedigree.
#'   * `MIDX` : An integer vector with maternal indices: For each $j=1,2,...$, `ID[MIDX[j]]` is the mother of `ID[j]`, or 0 if `ID[j]` has no mother within the pedigree.
#'   * `SEX` : An integer vector with gender codes. Unless the pedigree is reordered, this equals `as.integer(sex)`.
#'   * `FAMID` : The family ID.
#'   * `UNBROKEN_LOOPS` : A logical: TRUE if the pedigree is inbred.
#'   * `LOOP_BREAKERS` : A matrix with loop breaker ID's in the first
#'   column and their duplicates in the second column. All entries refer to the internal IDs.
#'   This is usually set by [breakLoops()].
#'   * `FOUNDER_INBREEDING` : A numeric vector with the same length as `founders(x)`, or NULL. This is always NULL when a new `ped` is created. See [founder_inbreeding()].
#'   * `MARKERS` : A list of `marker` objects.
#' @author Magnus Dehli Vigeland
#' @seealso [ped_basic], [ped_add], [ped_modify], [ped_subsets]
#'
#' @examples
#' x = ped(id=1:3, fid=c(0,0,1), mid=c(0,0,2), sex=c(1,2,1))
#' y = singleton('NN', sex=2, famid="SINGLETON GIRL")
#'
#' @export
ped = function(id, fid, mid, sex, famid=NULL, reorder = TRUE, check = TRUE, verbose = FALSE) {
  n = length(id)
  if(n ==0) stop2("`id` vector has length 0")
  if(length(fid) != n)
    stop2(sprintf("Incompatible input: length(id) = %d and length(fid) = %d", n, length(fid)))
  if(length(mid) != n)
    stop2(sprintf("Incompatible input: length(id) = %d and length(mid) = %d", n, length(mid)))
  if(length(fid) != n)
    stop2(sprintf("Incompatible input: length(id) = %d and length(sex) = %d", n, length(sex)))
  if (!all(sex %in% 0:2))
    stop2("Illegal gender code: ", setdiff(sex, 0:2))
  if(n == 1 && (fid!=0 || mid!=0))
    stop2("Singleton error: Parent IDs must be 0")

  # Parental index vectors (integer).
  FIDX = match(fid, id, nomatch=0L)
  MIDX = match(mid, id, nomatch=0L)

  # Initialise ped object
  x = list(ID = as.character(id), FIDX = FIDX, MIDX = MIDX, SEX = as.integer(sex),
           FAMID = if(is.null(famid)) "" else as.character(famid),
           UNBROKEN_LOOPS = FALSE, LOOP_BREAKERS = NULL, FOUNDER_INBREEDING = NULL,
           markerdata = NULL)

  class(x) = "ped"

  if(n == 1) {
    class(x) = c("singleton", class(x))
    return(x)
  }

  if (check) checkped(x)

  # Detect loops (by trying to find a peeling order)
  nucs = peelingOrder(x)
  lastnuc_link = nucs[[length(nucs)]]$link
  x$UNBROKEN_LOOPS = is.null(lastnuc_link)

  # reorder so that parents precede their children
  if(reorder) x = parents_before_children(x)

  x
}

#' @export
#' @rdname ped
singleton = function(id, sex = 1, famid = NULL) {
  if (length(id) != 1)
    stop2("Parameter `id` must have length 1")
  if (length(sex) != 1 || !sex %in% 0:2)
    stop2("Parameter `sex` must be either 0 (unknown), 1 (male) or 2 (female)")
  ped(id=id, fid=0, mid=0, sex=sex, famid=famid)
}

#' Pedigree errors
#'
#' Check a `ped` object for pedigree errors.
#'
#' @param x object of class `ped`.
#'
#' @return If no errors are detected, the function returns NULL invisibly.
#'   Otherwise, messages describing the errors are printed to the screen and an
#'   error is raised.
#'
#' @export
checkped = function(x) { #TODO: Some of this code doesn't make sense now - move to ped()
  ID = x$ID; FIDX = x$FIDX; MIDX = x$MIDX; SEX = x$SEX
  idx = seq_along(ID)
  if (length(idx) < 2) return()

  if (all(c(FIDX, MIDX) == 0))
      message("Pedigree is not connected.")

  fatherErr = !FIDX %in% c(0, idx)
  motherErr = !MIDX %in% c(0, idx)
  self_ancest = rep(F, length(idx)) # TODO: fix! sapply(seq_along(idx), function(i) idx[i] %in% ancestors(p, idx[i]))
  quick.check <- all(SEX %in% 0:2) &&
    all((FIDX > 0) == (MIDX > 0)) &&
    !any(duplicated(idx)) &&
    !any(fatherErr) &&
    !any(motherErr) &&
    all(SEX[match(FIDX[FIDX != 0], idx)] %in% c(0,1)) &&
    all(SEX[match(MIDX[MIDX != 0], idx)] %in% c(0,2)) &&
    !any(self_ancest)

  if (quick.check)
      return(invisible())  #if all tests are passed

  for (i in seq_along(idx)) {
    intro = sprintf("Individual %d: ", idx[i])
    if (!SEX[i] %in% 0:2)
      message(intro, "SEX must be either 0 (unknown), 1 (male) or 2 (female).")
    if ((FIDX[i] > 0) != (MIDX[i] > 0))
      message(intro, "Either both parents or none of them must be included.")
    if (duplicated(idx)[i]) {
      message(intro, "ID not unique.")
      next
    }
    if (fatherErr[i])
      message(intro, "Father's ID (", FIDX[i], ") does not appear in ID column.")
    else if (isTRUE(SEX[FIDX[i]] == 2))
      message(intro, "Father is female.")
    if (motherErr[i])
      message(intro, "Mother's ID (", MIDX[i], ") does not appear in ID column.")
    else if (isTRUE(SEX[MIDX[i]] == 1))
      message(intro, "Mother is male.")
    if (self_ancest[i])
      message(intro, "Is", switch(SEX[i]+1, "its", "his", "her"), "own ancestor.")
  }
  stop2("Pedigree errors detected")
}

# TODO
any_self_ancest = function(id, fid, mid) {
  is_self_anc = numeric(length(id))
  for(i in seq_along(id)) {
    anci = numeric()

  }
}
