#' Set a mutation model
#'
#' This function offers a convenient way to set or modify mutation models to
#' markers attached to a pedigree. It wraps [pedmut::mutationModel()], which
#' does the main work of creating the models, but relieves the user from having
#' to loop through the markers in order to supply the correct alleles and
#' frequencies for each marker. (This function supersedes
#' `pedprobr::setMutationModel()`.)
#'
#' Currently, the following models are implemented in the `pedmut` package:
#'
#' * `equal` :  All mutations equally likely; probability \eqn{1-rate} of no
#' mutation
#'
#' * `proportional` : Mutation probabilities are proportional to the target
#' allele frequencies
#'
#' * `onestep`: A mutation model for microsatellite markers, allowing mutations
#' only to the nearest neighbours in the allelic ladder. For example, '10' may
#' mutate to either '9' or '11', unless '10' is the lowest allele, in which case
#' '11' is the only option. This model is not applicable to loci with
#' non-integral microvariants.
#'
#' * `stepwise`: A common model in forensic genetics, allowing different
#' mutation rates between integer alleles (like '16') and non-integer
#' "microvariants" like '9.3'). Mutations also depend on the size of the
#' mutation if the parameter 'range' differs from 1.
#'
#' * `custom` : Allows any mutation matrix to be provided by the user, in the
#' `matrix` parameter
#'
#' * `random` : This produces a matrix of random numbers, where each row is
#' normalised so that it sums to 1
#'
#' * `trivial` : The identity matrix; i.e. no mutations are possible.
#'
#' @param x A `ped` object or a list of such.
#' @param markers A vector of names or indices referring to markers attached to
#'   `x`. (Default: All markers.)
#' @param ... Arguments forwarded to [pedmut::mutationModel()], e.g., `model`,
#'   `rate`, etc.
#' @param update A logical. If TRUE, existing mutation models (if present) are
#'   updated with the parameters specified in `...`. If FALSE (default), any
#'   previous models are ignored, and new mutation models are created from the
#'   parameters in `...`.
#'
#' @return An object similar to `x`.
#'
#' @examples
#'
#' ### Example requires the pedmut package ###
#'
#' if (requireNamespace("pedmut", quietly = TRUE)){
#'
#' # A pedigree with 1 empty marker; attach 'equal' mutation model
#' x = nuclearPed(1) |>
#'   addMarker() |>
#'   setMutmod(model = "equal", rate = 0.01)
#'
#' mutmod(x, 1)
#'
#' # Update rate (but still "equal" model)
#' y = setMutmod(x, rate = 0.05, update = TRUE)
#' mutmod(y, 1)
#'
#' # Change to stepwise model
#' z = setMutmod(x, model = "stepwise",
#'               rate = list(female = 0.01, male = 0.02),
#'               range = 0.1, rate2 = 1e-6)
#' mutmod(z, 1)
#'
#' # Remove mutation model
#' w = setMutmod(x, model = NULL)
#' mutmod(w, 1)
#'
#' }
#'
#' @importFrom pedmut getParams mutationModel
#' @export
setMutmod = function(x, markers = NULL, ..., update = FALSE) {
  if (!requireNamespace("pedmut", quietly = TRUE))
    stop2("Package `pedmut` must be installed in order to include mutation models")

  opts = list(...)
  markers = markers %||% seq_len(nMarkers(x))
  mIdx = whichMarkers(x, markers)

  # Remove all models?
  if("model" %in% names(opts) && is.null(opts$model)) {
    for(i in mIdx)
      mutmod(x, i) = NULL
    return(x)
  }

  for(i in mIdx) {

    useopts = opts

    if(update && !is.null(oldmut <- mutmod(x, i))) {
      oldpar = getParams(oldmut)
      oldparList = lapply(oldpar, function(v) list(female = v[1], male = v[2]))
      useopts = modifyList(oldparList, opts)
    }

    fr = afreq(x, i)
    args = c(list(alleles = names(fr), afreq = fr), useopts)

    mutmod(x, i) = do.call(pedmut::mutationModel, args)
  }

  x
}
