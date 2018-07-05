#' Modify the pedigree of 'ped' objects
#'
#' Functions to modify the pedigree of a 'ped' object.
#'
#' When removing an individual, all descendants are also removed as well as
#' founders remaining without offspring.
#'
#' The `branch()` function extracts the pedigree subset consisting of all
#' descendants of `id`, including `id` itself and all relevant
#' spouses.
#'
#' @param x A [ped()] object
#' @param id,ids Individual ID label(s). In `addOffspring` the (optional)
#' `ids` argument is used to specify ID labels for the offspring to be
#' created.
#' @param father,mother Integers indicating the IDs of parents. If missing, a
#' new founder individual is created (whose ID will be 1+the largest ID already
#' in the pedigree).
#' @param nch A single integer indicating the number of children to be
#' created. Default: 1.
#' @param sex Integer vectors indicating the genders
#' of the offspring to be created (recycled if less than `nch`
#' elements).
#' @param verbose A logical: Verbose output or not.
#' @param parent Integer ID of any pedigree member, which will be the father or
#' mother (depending on its gender) of the new child.
#' @param new a numeric containing new labels to replace those in `old`.
#' @param old a numeric containing ID labels to be replaced by those in
#' `new`. If missing, `old` is set to `x$LABELS`.
#' @return The modified `ped` object.
#' @author Magnus Dehli Vigeland
#' @seealso [ped()], [nuclearPed()]
#'
#' @examples
#'
#' x = nuclearPed(1)
#'
#' # To see the effect of each command below, use plot(x) in between.
#' x = addSon(x, 3)
#' x = addParents(x, id=4, father=6, mother=7)
#' x = swapSex(x, c(1,4))
#' x = removeIndividuals(x, 4)
#'
#' @name pedModify
NULL

#' @rdname pedModify
#' @export
swapSex = function(x, ids, verbose = TRUE) {
  assert_that(is.ped(x))
  if(!length(ids)) return(x)
  ids = internalID(x, ids)
  FID = x$FID
  MID = x$MID
  spouses = c(MID[FID %in% ids], FID[MID %in% ids])

  if (!all(spouses %in% ids)) {
    if (verbose) {
      extra = setdiff(spouses, ids)
      message("Changing sex of spouses as well: ", catLabels(x, extra))
    }
    return(swapSex(x, x$LABELS[union(ids, spouses)]))
  }

  # Swap sex
  x$SEX[ids] = 3 - x$SEX[ids]

  # # Swap parents wherever any of the 'ids' occur as parents
  ids_as_parents = FID %in% ids # same with MID!
  FID[ids_as_parents] = x$MID[ids_as_parents]
  MID[ids_as_parents] = x$FID[ids_as_parents]
  x$FID = FID
  x$MID = MID

  x
}

# CHANGE: addOffspring has changed name to addChildren.

#' @rdname pedModify
#' @export
addChildren = function(x, father=NULL, mother=NULL, nch = 1, sex = 1, ids = NULL, verbose = TRUE) {
  assert_that(is.ped(x), is.count(nch), all(sex %in% 0:2))
  father_exists = isTRUE(father %in% x$LABELS)
  mother_exists = isTRUE(mother %in% x$LABELS)
  if (!father_exists && !mother_exists)
    stop("At least one parent must be an existing pedigree member.", call. = FALSE)
  if (!is.null(ids) && length(ids) != nch)
    stop("Length of 'ids' must equal the number of children.", call. = FALSE)
  if (any(ids %in% x$LABELS))
    stop("Individuals already exist: ", catLabels(x, ids[ids %in% x$LABELS]), call.=FALSE)

  sex = rep_len(sex, nch)
  n = pedSize(x)

  # Check if labels are coercible to integers
  if(!has_numlabs(x))
    stop("non-integer labels not implemented", call. = FALSE)

  p = as.matrix(x)
  attrs = attributes(p)

  numlabs = as.numeric(x$LABELS)

  father_int = if(father_exists) internalID(x, father) else n + 1
  mother_int = if(mother_exists) internalID(x, mother) else n + 1

  if(!father_exists) {
    if(is.null(father)) {father = max(numlabs) + 1}
    if (verbose) message("Father: Creating new individual with ID ", father)
    p = rbind(p, c(father_int, 0, 0, 1))
    numlabs = c(numlabs, as.numeric(father))
  }
  if(!mother_exists) {
    if(is.null(mother)) {mother = max(numlabs) + 1}
    if (verbose) message("Mother: Creating new individual with ID ", mother)
    p = rbind(p, c(mother_int, 0, 0, 2))
    numlabs = c(numlabs, as.numeric(mother))
  }

  # Children
  numids = if(is.null(ids)) max(numlabs) + (1:nch) else suppressWarnings(as.numeric(ids))
  numlabs = c(numlabs, numids)

  if (isFALSE(all(numids == as.integer(numids))))
    stop("non-integer labels not implemented", call. = F)
  if (anyDuplicated(numlabs)) stop("Duplicated ID labels", call. = FALSE)

  p = rbind(p, cbind(nrow(p) + (1:nch), father_int, mother_int, sex))
  attrs$labels = as.character(numlabs)
  restore_ped(p, attrs = attrs)
}

#' @rdname pedModify
#' @export
addOffspring = addChildren


#' @rdname pedModify
#' @export
addSon = function(x, parent, id = NULL, verbose = TRUE) {
  parent_sex = getSex(x, parent)
  if (parent_sex == 1)
    addChildren(x, father = parent, nch = 1, sex = 1, ids = id, verbose = verbose)
  else if (parent_sex == 2)
    addChildren(x, mother = parent, nch = 1, sex = 1, ids = id, verbose = verbose)
  else stop("Not implemented for parents of unknown sex")
}

#' @rdname pedModify
#' @export
addDaughter = function(x, parent, id = NULL, verbose = TRUE) {
  parent_sex = getSex(x, parent)
  if (parent_sex == 1)
    addChildren(x, father = parent, nch = 1, sex = 2, ids = id, verbose = verbose)
  else if (parent_sex == 2)
    addChildren(x, mother = parent, nch = 1, sex = 2, ids = id, verbose = verbose)
  else stop("Not implemented for parents of unknown sex")
}

#' @rdname pedModify
#' @export
addParents = function(x, id, father=NULL, mother=NULL, verbose = TRUE) {
  if (length(id) > 1)
      stop("Only one individual at the time, please")

  id_int = internalID(x, id)
  if (id_int %in% x$NONFOUNDERS)
    stop("Individual ", id, " already has parents in the pedigree", call. = FALSE)

  # Check that assigned parents are OK
  desc = descendants(x, id)
  if (!is.null(father)) {
    if (father == id) stop("Pedigree error: Father and child are equal", call. = FALSE)
    if (father %in% desc) stop("Pedigree error: Assigned father is a descendant of ", id, call. = FALSE)
    if (father %in% x$LABELS && getSex(x, father) == 2) stop("Pedigree error: Assigned father is female", call. = FALSE)
  }
  if (!is.null(mother)) {
    if (mother == id) stop("Pedigree error: Mother and child are equal", call. = FALSE)
    if (mother %in% desc) stop("Pedigree error: Assigned mother is a descendant of ", id, call. = FALSE)
    if (mother %in% x$LABELS && getSex(x, mother) == 1) stop("Pedigree error: Assigned mother is male", call. = FALSE)
  }

  # If no labels are given, create them
  if (is.null(father) || is.null(mother)) {
    if(has_numlabs(x)) {
      mx = max(as.numeric(x$LABELS))
      if (is.null(father)) {father = mx + 1; mx = mx + 1}
      if (is.null(mother)) {mother = mx + 1}
    }
    else {
      if (is.null(father)) {father = nextNN(x$LABELS)}
      if (is.null(mother)) {mother = nextNN(c(x$LABELS, father))}
    }
  }

  p = as.matrix(x)
  attrs = attributes(p)
  labels = attrs$labels

  new.father = !father %in% x$LABELS
  new.mother = !mother %in% x$LABELS

  if (new.father) {
    if(verbose) message("Father: Creating new individual with ID ", father)
    fath_int = nrow(p) + 1
    p = rbind(p, c(fath_int, 0, 0, 1))
    labels = c(labels, father)
  }
  else fath_int = internalID(x, father)

  if (new.mother) {
    if(verbose) message("Mother: Creating new individual with ID ", mother)
    moth_int = nrow(p) + 1
    p = rbind(p, c(moth_int, 0, 0, 2))
    labels = c(labels, mother)
  }
  else moth_int = internalID(x, mother)

  p[id_int, 2:3] = c(fath_int, moth_int)
  attrs$labels = labels

  restore_ped(p, attrs = attrs)
}



#' @rdname pedModify
#' @export
removeIndividuals = function(x, ids, verbose = TRUE) {
  assert_that(is.ped(x))
  # Remove individuals 'ids' and all their descendants.
  # Spouse-founders are removed as well.
  if(!length(ids))
    return(x)

  ids_int = internalID(x, ids)

  # Founders without children after 'id' and 'desc' indivs are removed.
  # The redundancy in 'desc' does not matter.
  desc = numeric(0)
  for (id in ids_int) {
    dd = descendants(x, id, internal=T)
    desc = c(desc, dd)

    if (verbose) {
      if(length(dd)) {
        hisher = switch(x$SEX[id] + 1, "its", "his", "her")
        message("Removing ", x$LABELS[id], " and ", hisher, " descendants: ", catLabels(x, dd))
      }
      else message("Removing ", x$LABELS[id], " (no descendants)")
    }
  }

  # Keep: parents of remaining indivs (includes harmless zeroes)
  parents_of_remain = c(x$FID[-c(ids_int, desc)], x$MID[-c(ids_int, desc)])

  # But remove founders that are NOT among the above
  leftover_spouses = setdiff(x$FOUNDERS, c(ids_int, parents_of_remain))

  if (verbose && length(leftover_spouses))
    message("Removing leftover spouses: ", catLabels(x, leftover_spouses))

  # These are the ones to be removed (redundancy harmless)
  remov = c(ids_int, desc, leftover_spouses)

  # The actual reduction. Using the as.matrix trick anticipating marker data a.s.o.
  xmatr = as.matrix(x)
  new = xmatr[-remov, , drop=F]

  if(nrow(new) == 0) {if(verbose) message("Remaining pedigree is empty!"); return(NULL)}

  attrs = attributes(xmatr)
  attrs$labels = attrs$labels[-remov]

  restore_ped(new, attrs = attrs)
}

#' @rdname pedModify
#' @export
branch = function(x, id) {
  assert_that(is.ped(x))
  desc = descendants(x, id)
  spous = unlist(lapply(c(id, desc), spouses, x = x))
  subset(x, subset = c(id, desc, spous))
}


#' @rdname pedModify
#' @export
relabel = function(x, new, old=x$LABELS) {
  assert_that(is.ped(x), all(old %in% x$LABELS), length(new)==length(old))
  lab = x$LABELS
  lab[match(old, lab)] = new
  x$LABELS = lab
  x
}


#' Pedigree member labels
#'
#' Set labels of pedigree members
#'
#' @param x a `ped` object
#' @param labels a character (or coercible to character) of length `x$NIND`
#'
#' @export
setLabels = function(x, labels) {
  labels = as.character(labels)
  assert_that(is.ped(x), length(labels) == pedSize(x), !anyDuplicated(labels))
  x$LABELS = labels
  x
}

#' Family ID
#'
#' Set a family identifier of a `ped` object.
#'
#' @param x a `ped` object
#' @param famid a character of length 1. If missing, an emtpy string is used.
#'
#' @export
setFamid = function(x, famid) {
  famid = as.character(famid)
  if(!length(famid)) famid=""
  assert_that(is.ped(x), length(famid) == 1)
  x$FAMID = famid
  x
}
