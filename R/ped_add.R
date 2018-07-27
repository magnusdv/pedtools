#' Add/remove pedigree members
#'
#' Functions for adding or removing individuals in a 'ped' object.
#'
#' In `addChildren()` and `addParents()`, labels of added individuals are
#' created automatically if they are not specified by the user. In the automatic
#' case, the labelling depends on whether the existing labels are integer-like
#' or not. (To be precise, the program checks if `x$LABELS` equals
#' `as.character(as.integer(x$LABELS))`.) If so, the new labels are integers
#' subsequent to the largest of the existing labels. If not, the new labels are
#' "NN_1", "NN_2", ... If any such label already exists, the numbers are
#' adjusted accordingly.
#'
#' `addSon()` and `addDaughter()` are wrappers for a common use of
#' `addChildren()`, namely adding a single child to a pedigree member. Note that
#' its argument `parent` is gender-neutral, unlike in addChildren where you have
#' to know the parental genders. Also note that the other parent is always
#' created as a new individual. Thus, applying `addDaughter()` twice with the
#' same parent will create half sisters.
#'
#' In `removeIndividuals()` all descendants of `ids` are also removed. Any
#' individuals (spouses) left unconnected to the remaining pedigree are also
#' removed.
#'
#'
#' @param x A `ped` object.
#' @param ids A character vector (or coercible to such) with ID labels. In
#'   `addChildren` the (optional) `ids` argument is used to specify labels for
#'   the created children. If given, its length must equal `nch`. If not given,
#'   labels are assigned automatically as explained in Details.
#' @param id The ID label of some existing pedigree member.
#' @param father,mother Single ID labels. At least one of these must belong to
#'   an existing pedigree member. The other label may either: 1) belong to an
#'   existing member, 2) not belong to any existing member, or 3) be missing
#'   (i.e. not included in the function call). In cases 2 and 3 a new founder is
#'   added to the pedigree. In case 2 its label is the one given, while in case
#'   3 a suitable label is created by the program (see Details).
#' @param nch A single integer indicating the number of children to be created.
#'   Default: 1.
#' @param sex Gender codes of the created children (recycled if needed).
#' @param verbose A logical: Verbose output or not.
#' @param parent The ID label (coercible to character) of a single pedigree
#'   member, which will be the father or mother (depending on its gender) of the
#'   new child.
#'
#' @return The modified `ped` object.
#' @author Magnus Dehli Vigeland
#' @seealso [ped()], [ped_modify]
#'
#' @examples
#'
#' x = nuclearPed(1)
#'
#' # To see the effect of each command below, use plot(x) in between.
#' x = addSon(x, 3)
#' x = addParents(x, id=4, father=6, mother=7)
#' x = removeIndividuals(x, 4)
#'
#' @name ped_add
NULL

#' @rdname ped_add
#' @export
addChildren = function(x, father=NULL, mother=NULL, nch = 1, sex = 1, ids = NULL, verbose = TRUE) {
  assert_that(is.ped(x), is.count(nch), all(sex %in% 0:2))
  father_exists = isTRUE(father %in% x$LABELS)
  mother_exists = isTRUE(mother %in% x$LABELS)
  if (!father_exists && !mother_exists)
    stop2("At least one parent must be an existing pedigree member")
  if (!is.null(ids) && length(ids) != nch)
    stop2("Length of 'ids' must equal the number of children")
  if (any(ids %in% x$LABELS))
    stop2("Individuals already exist: ", intersect(ids, x$LABELS))

  # Recycle `sex` if needed
  sex = rep_len(sex, nch)

  # Prepare manipulation of pedigree matrix
  p = as.matrix(x)
  attrs = attributes(p)
  nmark = nMarkers(x)

  # Utility for creating new labels, depending on existing labels being numeric
  nextlabs = function(labs, len) {
    if(has_numlabs(x)) {
      mx = max(as.numeric(labs))
      seq.int(mx + 1, length.out=len)
    }
    else {
      res = character(0)
      for(i in seq_len(len)) res = c(res, nextNN(c(labs, res)))
      res
    }
  }

  # These variables will change as new members are created
  n = pedsize(x)
  labs = x$LABELS

  if(!father_exists) {
    if(is.null(father))
      father = nextlabs(labs, 1)
    if (verbose)
      message("Father: Creating new individual with ID = ", father)

    labs = c(labs, father)
    father_int = n = n + 1
    p = rbind(p, c(father_int, 0, 0, 1, rep(0, 2*nmark)))
  }
  else {
    father_int = internalID(x, father)
  }
  if(!mother_exists) {
    if(is.null(mother))
      mother = nextlabs(labs, 1)
    if (verbose)
      message("Mother: Creating new individual with ID = ", mother)

    labs = c(labs, mother)
    mother_int = n = n + 1
    p = rbind(p, c(mother_int, 0, 0, 2, rep(0, 2*nmark)))
  }
  else {
    mother_int = internalID(x, mother)
  }
  # Children
  if(is.null(ids)) ids = nextlabs(labs, len=nch)
  labs = c(labs, ids)

  if (anyDuplicated(labs)) stop2("Duplicated ID labels")

  children_pedcols = cbind(nrow(p) + (1:nch), father_int, mother_int, sex)
  children_markers = matrix(0, nrow=nch, ncol = nmark*2)
  p = rbind(p, cbind(children_pedcols, children_markers))

  attrs$labels = as.character(labs)
  restore_ped(p, attrs = attrs)
}


#' @rdname ped_add
#' @export
addSon = function(x, parent, id = NULL, verbose = TRUE) {
  parent_sex = getSex(x, parent)
  if (parent_sex == 1)
    addChildren(x, father = parent, nch = 1, sex = 1, ids = id, verbose = verbose)
  else if (parent_sex == 2)
    addChildren(x, mother = parent, nch = 1, sex = 1, ids = id, verbose = verbose)
  else stop2("Not implemented for parents of unknown sex: ", parent)
}

#' @rdname ped_add
#' @export
addDaughter = function(x, parent, id = NULL, verbose = TRUE) {
  parent_sex = getSex(x, parent)
  if (parent_sex == 1)
    addChildren(x, father = parent, nch = 1, sex = 2, ids = id, verbose = verbose)
  else if (parent_sex == 2)
    addChildren(x, mother = parent, nch = 1, sex = 2, ids = id, verbose = verbose)
  else stop2("Not implemented for parents of unknown sex: ", parent)
}

#' @rdname ped_add
#' @export
addParents = function(x, id, father=NULL, mother=NULL, verbose = TRUE) {
  if (length(id) > 1)
      stop2("Parents cannot be added to multiple individuals at once: ", id)
  if (id %in% nonfounders(x))
    stop2("Individual ", id, " already has parents in the pedigree")

  id_int = internalID(x, id)
  
  # Check that assigned parents are OK
  desc = descendants(x, id)
  if (!is.null(father)) {
    if (father == id) stop2("Father and child are equal")
    if (father %in% desc) stop2("Assigned father is a descendant of ", id)
    if (father %in% x$LABELS && getSex(x, father) == 2) stop2("Assigned father is female")
  }
  if (!is.null(mother)) {
    if (mother == id) stop2("Mother and child are equal")
    if (mother %in% desc) stop2("Assigned mother is a descendant of ", id)
    if (mother %in% x$LABELS && getSex(x, mother) == 1) stop2("Assigned mother is male")
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
  nmark = nMarkers(x)

  new.father = !father %in% x$LABELS
  new.mother = !mother %in% x$LABELS

  if (new.father) {
    if(verbose) message("Father: Creating new individual with ID = ", father)
    fath_int = nrow(p) + 1
    p = rbind(p, c(fath_int, 0, 0, 1, rep(0, 2*nmark)))
    labels = c(labels, father)
  }
  else fath_int = internalID(x, father)

  if (new.mother) {
    if(verbose) message("Mother: Creating new individual with ID = ", mother)
    moth_int = nrow(p) + 1
    p = rbind(p, c(moth_int, 0, 0, 2, rep(0, 2*nmark)))
    labels = c(labels, mother)
  }
  else moth_int = internalID(x, mother)

  p[id_int, 2:3] = c(fath_int, moth_int)
  attrs$labels = labels

  restore_ped(p, attrs = attrs)
}



#' @rdname ped_add
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
        message("Removing ", x$LABELS[id], " and ", hisher, " descendants: ", toString(x$LABELS[dd]))
      }
      else message("Removing ", x$LABELS[id], " (no descendants)")
    }
  }

  # Keep: parents of remaining indivs (includes harmless zeroes)
  parents_of_remain = c(x$FID[-c(ids_int, desc)], x$MID[-c(ids_int, desc)])

  # But remove founders that are NOT among the above
  FOU = founders(x, internal=T)
  leftover_spouses = setdiff(FOU, c(ids_int, parents_of_remain))

  if (verbose && length(leftover_spouses))
    message("Removing leftover spouses: ", toString(x$LABELS[leftover_spouses]))

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

# #' @rdname ped_add
# #' @export
branch = function(x, id) {
  assert_that(is.ped(x))
  desc = descendants(x, id)
  spous = unlist(lapply(c(id, desc), spouses, x = x))
  subset(x, subset = c(id, desc, spous))
}

