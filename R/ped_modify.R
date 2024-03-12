#' Add/remove pedigree members
#'
#' Functions for adding or removing individuals in a 'ped' object.
#'
#' In `addChildren()` and `addParents()`, labels of added individuals are
#' generated automatically if they are not specified by the user. If the
#' existing labels are all integer-like, the generated labels are the smallest
#' available integers. Otherwise, the new labels are the first available in the
#' sequence "a1", "a2", ...
#'
#' `addChild()`, `addSon()` and `addDaughter()` are convenient wrappers for the
#' most common use of `addChildren()`, namely adding a single child to a
#' pedigree. Note that the parents can be given in any order. If only one parent
#' is supplied, the other is created as a new individual.
#'
#' `removeIndividuals()` removes the individuals indicated with `ids` along with
#' all of their ancestors OR descendants, depending on the `remove` argument.
#' Leftover spouses disconnected to the remaining pedigree are also removed. An
#' error is raised if result is a disconnected pedigree.
#'
#' The `branch()` function extracts the sub-pedigree formed by `id` and all
#' his/her spouses and descendants.
#'
#' Finally, `subset()` can be used to extract any connected sub-pedigree. (Note
#' that in the current implementation, the function does not actually check that
#' the indicated subset forms a connected pedigree; failing to comply with this
#' may lead to obscure errors.)
#'
#' @param x A `ped` object, or a list of such.
#' @param ids A vector of ID labels. In `addChildren()` these are the children
#'   to be created. If NULL (default) given, automatic labels are generated.
#' @param remove Either "ancestors" or "descendants" (default), dictating the
#'   method of removing pedigree members. Abbreviations are allowed.
#' @param returnLabs A logical, by default FALSE. If TRUE, `removeIndividuals()`
#'   returns only the labels of all members to be removed, instead of actually
#'   removing them.
#' @param id The ID label of a pedigree member.
#' @param father,mother Single ID labels. At least one of these must be an
#'   existing member of `x`. The other may be (i) another existing member, (ii)
#'   a new founder to be created, or (iii) missing (i.e., NULL), in which case
#'   the other parent is created and given a suitable name.
#' @param nch A positive integer indicating the number of children to be
#'   created. Default: 1.
#' @param sex Gender codes of the created children (recycled if needed).
#' @param verbose A logical: Verbose output or not.
#' @param parents A vector of 1 or 2 ID labels, of which at least one must be an
#'   existing member of `x`.
#'
#' @return The modified `ped` object.
#' @seealso [ped()], [relabel()], [swapSex()]
#'
#' @examples
#'
#' x = nuclearPed(1) |>
#'   addSon(3) |>
#'   addParents(4, father = 6, mother = 7) |>
#'   addChildren(father = 6, mother = 7, nch = 3, sex = c(2,1,2))
#'
#' # Remove 6 and 7 and their descendants
#' y1 = removeIndividuals(x, 6:7)
#'
#' # Remove 8-10 and their parents
#' y2 = removeIndividuals(x, 8:10, remove = "ancestors")
#'
#' # Adding a child across components
#' z = singletons(1:2, sex = 1:2) |> addDaughter(1:2)
#'
#'
#' @name ped_modify
NULL

#' @rdname ped_modify
#' @export
addChildren = function(x, father = NULL, mother = NULL, nch = NULL, sex = 1L, ids = NULL, verbose = TRUE) {
  islist = is.pedList(x)
  if(!is.ped(x) && !islist)
    stop2("Input is not a `ped` object or a list of such")

  # NB! Labels will change as new members are created
  labs = labels(x)
  isnum = .isIntegral(labs)

  # Check input
  if(length(father) > 1)
    stop2("More than one father indicated: ", father)
  if(length(mother) > 1)
    stop2("More than one mother indicated: ", mother)

  father = father %||% generateLabs(labs, n=1, num = isnum, avoid = c(mother, ids))
  if(!(father_exists <- father %in% labs))
    labs = c(labs, father)

  mother = mother %||% generateLabs(labs, n=1, num = isnum, avoid = ids)
  if(!(mother_exists <- mother %in% labs))
    labs = c(labs, mother)

  if(!father_exists && !mother_exists)
    stop2("At least one parent must be an existing pedigree member")

  if(father_exists && getSex(x, father) == 2)
    stop2("Assigned father is female: ", father)
  if(mother_exists && getSex(x, mother) == 1)
    stop2("Assigned mother is male: ", mother)

  # Number of children
  nch = nch %||% if(!is.null(ids)) length(ids) else length(sex)
  if(!isCount(nch))
    stop2("Argument `nch` must be a positive integer: ", nch)

  # Children IDs
  ids = ids %||% generateLabs(labs, n = nch, num = isnum)
  if(length(ids) != nch)
    stop2("Length of `ids` must equal the number of children")
  if(any(ids %in% labs))
    stop2("Individual already exists: ", intersect(ids, labs))
  if (anyDuplicated.default(labs))
    stop2("Duplicated ID label: ", labs[duplicated(labs)])

  labs = c(labs, ids)

  # Check `sex` and recycle if needed
  if(!is.numeric(sex) || !all(sex %in% 0:2))
    stop2("Illegal value of `sex`: ", .mysetdiff(sex, 0:2))
  sex = rep_len(as.integer(sex), nch)


  if(islist) {
    comp = getComponent(x, c(father, mother), checkUnique = TRUE)
    compSet = unique.default(comp[!is.na(comp)])
    if(length(compSet) == 1) {
      newcomp = addChildren(x[[compSet]], father, mother, nch = nch, sex = sex,
                            ids = ids, verbose = verbose)
      x[[compSet]] = newcomp
    }
    else if(length(compSet) == 2) {
      facomp = compSet[1]
      mocomp = compSet[2]
      newcomp = x[[facomp]] |>
        addChildren(father, mother, nch = nch, sex = sex, ids = ids, verbose = FALSE) |>
        mergePed(x[[mocomp]], by = mother, relabel = FALSE)
      x[[facomp]] = newcomp
      x[[mocomp]] = NULL
      if(length(x) == 1) x = x[[1]]
    }
    return(x)
  }

  n = pedsize(x)

  # Prepare manipulation of pedigree matrix
  p = as.matrix(x)
  attrs = attributes(p)
  nmark = nMarkers(x)

  if(!father_exists) {
    if(verbose) message("Creating new father: ", father)
    father_int = n = n + 1L
    p = rbind(p, c(father_int, 0L, 0L, 1L, rep(0L, 2*nmark)))
  }
  else {
    father_int = internalID(x, father)
  }

  if(!mother_exists) {
    if (verbose) message("Creating new mother: ", mother)
    mother_int = n = n + 1L
    p = rbind(p, c(mother_int, 0L, 0L, 2L, rep(0L, 2*nmark)))
  }
  else {
    mother_int = internalID(x, mother)
  }

  children_pedcols = cbind(nrow(p) + (1:nch), father_int, mother_int, sex)
  children_markers = matrix(0L, nrow = nch, ncol = nmark*2)
  p = rbind(p, cbind(children_pedcols, children_markers))

  attrs$LABELS = as.character(labs)
  restorePed(p, attrs = attrs)
}


#' @rdname ped_modify
#' @export
addChild = function(x, parents, id = NULL, sex = 1, verbose = TRUE) {
  # Convenience function for a single child of sex 0, 1 or 2.
  # Note: The implementation actually allows multiple children

  npar = length(parents)
  if(npar == 0)
    stop2("No parents indicated")

  if(npar > 2)
    stop2("Too many parents indicated: ", parents)

  if(npar == 2 && parents[1] == parents[2])
    stop2("Duplicated parent: ", parents[1])

  parents = as.character(parents) # remove potential names etc
  existing = parents %in% labels(x)
  if(!any(existing))
    stop2("At least one parent must be an existing pedigree member: ", parents)

  par1 = if(existing[1]) parents[1] else parents[2]
  par2 = if(npar == 2) setdiff(parents, par1) else NULL

  sex1 = getSex(x, par1)
  if(sex1 == 1)      { fa = par1; mo = par2 }
  else if(sex1 == 2) { mo = par1; fa = par2 }
  else
    stop2("Not implemented for parents of unknown sex: ", par1)

  nch = length(id %||% sex) # default 1
  addChildren(x, father = fa, mother = mo, nch = nch, sex = sex, ids = id, verbose = verbose)
}

#' @rdname ped_modify
#' @export
addSon = function(x, parents, id = NULL, verbose = TRUE) {
  addChild(x, parents = parents, id = id, sex = 1, verbose = verbose)
}

#' @rdname ped_modify
#' @export
addDaughter = function(x, parents, id = NULL, verbose = TRUE) {
  addChild(x, parents = parents, id = id, sex = 2, verbose = verbose)
}

#' @rdname ped_modify
#' @export
addParents = function(x, id, father = NULL, mother = NULL, verbose = TRUE) {
  if (length(id) > 1)
      stop2("Cannot add parents to multiple individuals: ", id)

  labs = labels(x)
  isnum = .isIntegral(labs)

  father = father %||% generateLabs(labs, n = 1, avoid = mother, num = isnum)
  mother = mother %||% generateLabs(labs, n = 1, avoid = father, num = isnum)

  fatherExists = father %in% labs
  motherExists = mother %in% labs

  if(is.pedList(x)) {
    if(fatherExists || motherExists)
      stop2("The current implementation doesn't support existing parents in ped lists")

    comp = getComponent(x, id, checkUnique = TRUE, errorIfUnknown = TRUE)
    x[[comp]] = addParents(x[[comp]], id, father, mother, verbose = verbose)
    return(x)
  }

  ### By now, x is a connected pedigree

  if (id %in% nonfounders(x))
    stop2(sprintf("Individual '%s' already has parents in the pedigree", id))

  id_int = internalID(x, id)

  # Check that assigned parents are OK
  desc = descendants(x, id)
  if (fatherExists) {
    if (father == id) stop2("Father and child are equal")
    if (father %in% desc) stop2("Assigned father is a descendant of ", id)
    if (getSex(x, father) == 2) stop2("Assigned father is female")
  }
  if (motherExists) {
    if (mother == id) stop2("Mother and child are equal")
    if (mother %in% desc) stop2("Assigned mother is a descendant of ", id)
    if (getSex(x, mother) == 1) stop2("Assigned mother is male")
  }


  p = as.matrix(x)
  attrs = attributes(p)
  newlabs = attrs$LABELS
  nmark = nMarkers(x)

  if (!fatherExists) {
    if(verbose) message("Creating new father: ", father)
    fath_int = nrow(p) + 1
    p = rbind(p, c(fath_int, 0, 0, 1, rep(0, 2*nmark)))
    newlabs = c(newlabs, father)
  }
  else fath_int = internalID(x, father)

  if (!motherExists) {
    if(verbose) message("Creating new mother: ", mother)
    moth_int = nrow(p) + 1
    p = rbind(p, c(moth_int, 0, 0, 2, rep(0, 2*nmark)))
    newlabs = c(newlabs, mother)
  }
  else moth_int = internalID(x, mother)

  p[id_int, 2:3] = c(fath_int, moth_int)
  attrs$LABELS = newlabs

  y = restorePed(p, attrs = attrs)
  neworder = c(seq_len(id_int - 1),
               if(!fatherExists) fath_int,
               if(!motherExists) moth_int,
               id_int:pedsize(x))
  reorderPed(y, neworder)
}



#' @rdname ped_modify
#' @export
removeIndividuals = function(x, ids, remove = c("descendants", "ancestors"),
                             returnLabs = FALSE, verbose = TRUE) {
  if(is.pedList(x)) {
    y = lapply(x, function(comp)
      removeIndividuals(comp, .myintersect(comp$ID, ids), remove = remove,
                        returnLabs = returnLabs, verbose = FALSE))
    # Remove NULLs
    y = y[!sapply(y, is.null)]
    return(y)
  }

  if(!is.ped(x))
    stop2("Input is not a `ped` object or a list of such")

  if(!length(ids))
    return(x)

  if(anyDuplicated.default(ids))
    ids = unique.default(ids)

  ids_int = internalID(x, ids, errorIfUnknown = TRUE)
  labs = x$ID

  remove = match.arg(remove)

  if(verbose)
    message(sprintf("Removing individual(s) and their %s: %s", remove, toString(ids)))

  # Descendants OR ancestors.
  remov = switch(remove,
                ancestors = ancestors(x, ids_int, inclusive = TRUE, internal = TRUE),
                descendants = descendants(x, ids_int, inclusive = TRUE, internal = TRUE))

  makeFounderIdx = integer(0)

  if(remove == "descendants") {
    # Remove founders that are NOT parents of any remaining
    FOU = founders(x, internal = TRUE)
    leftovers = .mysetdiff(FOU, c(remov, x$FIDX[-remov], x$MIDX[-remov]))

    if(length(leftovers))
      remov = unique.default(c(remov, leftovers))
  }
  else if(remove == "ancestors") {
    while(TRUE) {
      # Remaining parents (idx) - includes zeroes
      remainFidx = x$FIDX[-remov]
      remainMidx = x$MIDX[-remov]

      # Remaining individuals for which at least 1 parent is removed
      makeFounder = remainFidx %in% remov | remainMidx %in% remov
      makeFounderIdx = seq_along(x$ID)[-remov][makeFounder]

      # Remove spouses who become unattached
      sp = .mysetdiff(c(remainFidx[makeFounder], remainMidx[makeFounder]), remov)
      isParent = sp %in% c(remainFidx[!makeFounder], remainMidx[!makeFounder])
      isChild = !x$FIDX[sp] %in% c(0, remov) & !x$MIDX[sp] %in% c(0, remov)
      leftovers = sp[!isParent & !isChild]

      if(length(leftovers))
        remov = unique.default(c(remov, leftovers))
      else
        break
    }
  }
  if(returnLabs)
    return(labs[.mysortInt(remov)])

  if(verbose) {
    message(sprintf("Total number to be removed: %d. (Remaining: %d)",
                    length(remov), length(labs) - length(remov)))
    if(length(makeFounderIdx))
      message("Converting to founder: ", toString(labs[makeFounderIdx]))
  }

  # The actual reduction. Using the as.matrix trick anticipating marker data a.s.o.
  xmatr = as.matrix(x)
  new = xmatr[-remov, , drop = FALSE]

  if(nrow(new) == 0) {
    if(verbose) message("Remaining pedigree is empty!")
    return(NULL)
  }

  if(length(makeFounderIdx))
    new[match(makeFounderIdx, new[,1]), 2:3] = 0

  attrs = attributes(xmatr)

  # Remove labels
  attrs$LABELS = attrs$LABELS[-remov]

  # Remove founder inbreeding
  finb = attrs$FOUNDER_INBREEDING
  if (!is.null(finb) && length(intersect(remov, FOU))) {
    aut = finb$autosomal
    xchr = finb$x
    attrs$FOUNDER_INBREEDING = list(autosomal = aut[!names(aut) %in% labs[remov]],
                                            x = xchr[!names(xchr) %in% labs[remov]])
  }

  restorePed(new, attrs = attrs)
}

#' @rdname ped_modify
#' @export
branch = function(x, id) {
  if(!is.ped(x)) stop2("Input is not a `ped` object")
  if(length(id) == 0)
    stop2("`id` cannot be empty")
  if(length(id) > 1)
    stop2("`id` must contain a single ID label: ", id)
  desc = descendants(x, id)
  spous = unlist(lapply(c(id, desc), spouses, x = x))
  ids = unique.default(c(id, desc, spous))

  # sort (since subset() does not sort)
  ids = ids[order(internalID(x, ids))]

  subset(x, subset = ids)
}


#' @param subset A character vector (or coercible to such) with ID labels
#'   forming a connected sub-pedigree.
#' @param ... Not used.
#'
#' @rdname ped_modify
#' @export
subset.ped = function(x, subset, ...) {
  if(!is.null(x$LOOP_BREAKERS))
    stop2("`subset()` is not yet implemented for pedigrees with broken loops")

  sub_idx = internalID(x, subset)

  if(anyDuplicated.default(subset))
    stop2("Duplicated ID label: ", unique(subset[duplicated(subset)]))

  pedm = as.matrix(x)
  subped = pedm[sub_idx, , drop = FALSE]

  # set FID = 0 if father is not in subset
  subped[!(subped[, 2] %in% sub_idx), 2] = 0L

  # set MID = 0 if mother is not in subset
  subped[!(subped[, 3] %in% sub_idx), 3] = 0L

  # Fix labels
  attrs = attributes(pedm)
  attrs$LABELS = attrs$LABELS[sub_idx]

  restorePed(subped, attrs, ...)
}

