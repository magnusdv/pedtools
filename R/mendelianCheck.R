#' Check for Mendelian errors
#'
#' Check marker data for Mendelian inconsistencies
#'
#' @param x a [ped()] object
#' @param remove a logical. If FALSE, the function returns the indices of
#' markers found to incorrect.  If TRUE, a new `ped` object is
#' returned, where the incorrect markers have been deleted.
#' @param verbose a logical. If TRUE, details of the markers failing the tests
#' are shown.
#'
#' @return
#' A numeric containing the indices of the markers
#' that did not pass all tests, or (if `remove = TRUE`) a new `ped`
#' object where the failing markers are removed.
#' @author Magnus Dehli Vigeland
#'
#' @examples
#'
#' x = nuclearPed()
#'
#' # Add a SNP with Mendelian error
#' m = marker(x, '1' = "1/1", '2' = "1/1", '3' = "1/2")
#' x = setMarkers(x, m)
#'
#' mendelianCheck(x)
#'
#' @export
mendelianCheck = function(x, remove = FALSE, verbose = !remove) {
  if (is.singleton(x) || !hasMarkers(x))
    return(if(remove) x else numeric(0))

  nucs = subnucs(x)

  chromX = which(isXmarker(x))
  chromAUT = .mysetdiff(seq_markers(x), chromX)

  errorlist = vector(length = pedsize(x), mode = "list")
  names(errorlist) = labels(x)
  nuc_errors = numeric()  # container for allele count errors...belongs to the whole subnuc.

  ### AUTOSOMAL
  if (length(chromAUT) > 0) {
    allelematr = do.call(cbind, x$MARKERS[chromAUT])
    for (sub in nucs) {
      fa = allelematr[sub$father, ]
      mo = allelematr[sub$mother, ]
      offs = sub$children
      labs = attr(sub, 'labels')

      # Check trios
      for (of in offs) {
        triopass = trioCheckFast(fa, mo, allelematr[of, ])
        if (all(triopass))
          next

        new_errs = chromAUT[!triopass]
        errorlist[[of]] = c(errorlist[[of]], new_errs)
        if (verbose)
          message(sprintf("Individual `%s` incompatible with parents at markers: %s",
                          labs[of], toString(new_errs)))
      }

      # Check sibshibs
      if (length(offs) == 1)
        next

      sibpass = sibshipCheck(allelematr[offs, ])
      if(all(sibpass))
        next

      new_errs = chromAUT[!sibpass]
      nuc_errors = c(nuc_errors, new_errs)
      if (verbose)
        message(sprintf("Sibship with parents `%s` and `%s` have too many alleles at markers: %s",
                        labs[sub$father], labs[sub$mother], toString(new_errs)))
    }
  }

  ### X
  if (length(chromX) > 0) {
    sex = x$SEX
    allelematr = do.call(cbind, x$MARKERS[chromX])

    # Identify & report male heterozygosity
    even = 2 * seq_along(chromX)
    odd = even - 1
    maleXhet = allelematr[sex == 1, odd, drop = FALSE] != allelematr[sex == 1, even, drop = FALSE]
    if(any(maleXhet)) {
      maleXhet_errors = which(maleXhet, arr.ind = TRUE)
      error_males_int = which(sex == 1)[maleXhet_errors[, 1]] # modify first col from index *among males*, to *among all*
      error_markers = maleXhet_errors[, 2]
      for(i in unique.default(error_males_int)) {
        new_errs = chromX[error_markers[error_males_int == i]]
        errorlist[[i]] = c(errorlist[[i]], new_errs)
        if (verbose)
          message(sprintf("Male `%s` is heterozygous at markers: %s", labels(x)[i], toString(new_errs)))
      }
    }

    for (sub in nucs) {
      fa = allelematr[sub$father, ]
      mo = allelematr[sub$mother, ]
      offs = sub$children
      labs = attr(sub, 'labels')

      # Check trios - X
      for (of in offs) {
        ofdat = allelematr[of, ]
        if (sex[of] == 1) {
          maleXpass = maleXCheck(mo, ofdat)
          if (all(maleXpass))
            next
          new_errs = chromX[!maleXpass]
          errorlist[[of]] = c(errorlist[[of]], new_errs)
          if (verbose)
            message(sprintf("Male `%s` incompatible with mother at markers: %s",
                            labs[of], toString(new_errs)))
        }
        else {
          triopass = trioCheckFast(fa, mo, ofdat)
          if (all(triopass))
            next

          new_errs = chromX[!triopass]
          errorlist[[of]] = c(errorlist[[of]], new_errs)
          if (verbose)
            message(sprintf("Female `%s` incompatible with parents at markers: %s",
                            labs[of], toString(new_errs)))
        }
      }

      # Sibships - X
      # TODO: skip this for now - needs fixing!
      next

      if (length(offs) <= 2)
        next

      sibpass = sibshipCheck(allelematr[offs, ])
      if(all(sibpass))
        next

      new_errs = chromX[!sibpass]
      nuc_errors = c(nuc_errors, new_errs)
      if (verbose)
        message(sprintf("Sibship with parents `%s` and `%s` have too many alleles at markers: %s",
                        labs[sub$father], labs[sub$mother], toString(new_errs)))
    }
  }
  err_index = sort.int(unique.default(c(unlist(errorlist), nuc_errors)))

  if (remove)
      return(removeMarkers(x, err_index)) else return(err_index)
}


trioCheckFast = function(fa, mo, of) {
  even = 2 * seq_len(length(fa)/2)
  odd = even - 1
  fa_odd = fa[odd]
  fa_even = fa[even]
  mo_odd = mo[odd]
  mo_even = mo[even]
  of_odd = of[odd]
  of_even = of[even]
  fa0 = (fa_odd == 0 | fa_even == 0)
  mo0 = (mo_odd == 0 | mo_even == 0)
  of_odd0 = (of_odd == 0)
  of_even0 = (of_even == 0)
  ff1 = (fa0 | of_odd0 | of_odd == fa_odd | of_odd == fa_even)
  ff2 = (fa0 | of_even0 | of_even == fa_odd | of_even == fa_even)
  mm1 = (mo0 | of_odd0 | of_odd == mo_odd | of_odd == mo_even)
  mm2 = (mo0 | of_even0 | of_even == mo_odd | of_even == mo_even)
  (ff1 & mm2) | (ff2 & mm1)
}

maleXHomoz = function(of) {
  even = 2 * seq_len(length(of)/2)
  odd = even - 1
  of[odd] == of[even]
}

maleXCheck = function(mo, of) {
  even = 2 * seq_len(length(of)/2)
  odd = even - 1
  mo_odd = mo[odd]
  mo_even = mo[even]
  of = of[odd]
  mo0 = (mo_odd == 0 | mo_even == 0)
  of == 0 | mo0 | of == mo_odd | of == mo_even
}

sibshipCheck = function(offs) {
  # offs = matrix with 2*N columns
  even = 2 * seq_len(ncol(offs)/2)

  # loop through markers
  unlist(lapply(even, function(i) {
    ###offs_als = unique.default(offs[, (i - 1):i])
    genos = offs[, (i - 1):i]

    # number of (different) alleles occuring in homozygous state
    homoz_alleles = genos[genos[,1] == genos[,2], 1]
    n_homoz = length(.mysetdiff(homoz_alleles, 0))

    # number of different alleles in total
    n_alleles = length(.mysetdiff(genos, 0))

    # if no homoz: consistent if number of alleles <= 4.
    # for each new allele observed homoz, "4" is reduced by 1.
    n_alleles <= (4 - n_homoz)
  }))
}
