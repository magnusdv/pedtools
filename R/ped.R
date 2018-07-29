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
#' @param famid a character of length 1
#' @param reorder a logical. If TRUE, the pedigree is reordered so that all
#'   parents preceeede their children.
#' @param check a logical. If TRUE, [checkped()] is run on the
#'   pedigree before it is returned.
#' @param verbose a logical.
#' @param ... further arguments
#'
#' @return A `ped` object, which is essentially a list with the following
#'   entries:
#'   \describe{
#'   \item{ID}{A numerical vector with the internal IDs, which are always 1,2,3,...,N.}
#'   \item{FID}{A numerical vector indicating the internal IDs of the fathers.}
#'   \item{MID}{A numerical vector indicating the internal IDs of the mothers.}
#'   \item{SEX}{A numerical vector with gender codes. Unless the pedigree is reordered,
#'     this equals the input argument `sex`.}
#'   \item{FAMID}{The family ID.}
#'   \item{LABELS}{A character vector containing the original id labels.
#'     Unless the pedigree has been reordered, this equals the input argument `id`.}
#'   \item{UNBROKEN_LOOPS}{A logical: TRUE if the pedigree is inbred.}
#'   \item{LOOP_BREAKERS}{A matrix with loop breaker ID's in the first
#'   column and their duplicates in the second column. All entries refer to the internal IDs.
#'   This is usually set by [breakLoops()].}
#'   }
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
  assert_that(n>0, length(fid)==n, length(mid)==n, length(sex)==n)
  if(n == 1 && (fid!=0 || mid!=0))
    stop2("Singleton error: Parent IDs must be 0")

  # Internal order 1,2,...
  ID = 1:n
  FID = match(fid, id, nomatch=0)
  MID = match(mid, id, nomatch=0)

  # Initialise ped object
  x = list(ID = ID, FID = FID, MID = MID, SEX = sex,
           LABELS = NULL, FAMID = NULL,
           UNBROKEN_LOOPS = FALSE, LOOP_BREAKERS = NULL,
           markerdata = NULL)

  class(x) = "ped"
  x = setFamid(x, famid)
  x = setLabels(x, id)

  if(n == 1) {
    class(x) = c("singleton", class(x))
    return(x)
  }

  if (check) checkped(x)

  # Detect loops (by trying to find a peeling order)
  nucs = peelingOrder(x)
  lastnuc_link = nucs[[length(nucs)]]$link
  x$UNBROKEN_LOOPS = is.null(lastnuc_link)

  # reorder so that parents preceede their children
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
checkped = function(x) {
  ID=x$ID; FID=x$FID; MID=x$MID; SEX=x$SEX
  if (length(ID) < 2) return()

  if (all(c(FID, MID) == 0))
      message("Pedigree is not connected.")

  fatherErr = !FID %in% c(0, ID)
  motherErr = !MID %in% c(0, ID)
  self_ancest = rep(F, length(ID)) # TODO: fix! sapply(seq_along(ID), function(i) ID[i] %in% ancestors(p, ID[i]))
  quick.check <- all(SEX %in% 0:2) &&
    all((FID > 0) == (MID > 0)) &&
    !any(duplicated(ID)) &&
    !any(fatherErr) &&
    !any(motherErr) &&
    all(SEX[match(FID[FID != 0], ID)] %in% c(0,1)) &&
    all(SEX[match(MID[MID != 0], ID)] %in% c(0,2)) &&
    !any(self_ancest)

  if (quick.check)
      return(invisible())  #if all tests are passed

  for (i in seq_along(ID)) {
    intro = sprintf("Individual %d: ", ID[i])
    if (!SEX[i] %in% 0:2)
      message(intro, "SEX must be either 0 (unknown), 1 (male) or 2 (female).")
    if ((FID[i] > 0) != (MID[i] > 0))
      message(intro, "Either both parents or none of them must be included.")
    if (duplicated(ID)[i]) {
      message(intro, "ID not unique.")
      next
    }
    if (fatherErr[i])
      message(intro, "Father's ID (", FID[i], ") does not appear in ID column.")
    else if (isTRUE(SEX[FID[i]] == 2))
      message(intro, "Father is female.")
    if (motherErr[i])
      message(intro, "Mother's ID (", MID[i], ") does not appear in ID column.")
    else if (isTRUE(SEX[MID[i]] == 1))
      message(intro, "Mother is male.")
    if (self_ancest[i])
      message(intro, "Is", switch(SEX[i]+1, "its", "his", "her"), "own ancestor.")
  }
  stop2("Pedigree errors detected")
}

#TODO
as_ped.matrix = function(m) ped(id=m[,1], fid=m[,2], mid=m[,3], sex=m[,4])


# TODO
any_self_ancest = function(id, fid, mid) {
  is_self_anc = numeric(length(id))
  for(i in seq_along(id)) {
    anci = numeric()

  }
}
