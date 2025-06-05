#' Set a mutation model
#'
#' This function offers a convenient way to set or modify mutation models to
#' markers attached to a pedigree. It wraps [pedmut::mutationModel()], which
#' does the main work of creating the models, but relieves the user from having
#' to loop through the markers in order to supply the correct alleles and
#' frequencies for each marker.
#'
#' Currently, the following models are supported:
#'
#' * `equal`:  All mutations equally likely; probability `1 - rate` of no
#' mutation
#'
#' * `proportional`: Mutation probabilities are proportional to the target
#' allele frequencies
#'
#' * `onestep`: A simple model for microsatellite markers, in which mutations
#' are only allowed to the nearest neighbours in the allelic ladder. For
#' example, '10' may mutate to either '9' or '11' (unless '10' is the lowest
#' allele, in which case '11' is the only option). Not applicable to loci with
#' non-integral microvariants.
#'
#' * `stepwise`: A common model for microsatellite markers. Mutation rates
#' depend on the step size in the allelic ladder, and also the allelic classes:
#' integral repeats like '16', versus non-integer microvariants like '16.3'.
#'
#' * `custom`: Allows any mutation matrix to be provided by the user, in the
#' `matrix` parameter
#'
#' * `random`: This produces a matrix of random numbers, where each row is
#' normalised so that it sums to 1
#'
#' * `trivial`: The identity matrix; no mutations are possible
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

  # Remove all models?
  if("model" %in% names(opts) && is.null(opts$model)) {
    mutmod(x, markers) = NULL
    return(x)
  }

  markers = markers %||% seq_len(nMarkers(x))
  mIdx = whichMarkers(x, markers)

  for(i in mIdx) {

    useopts = opts

    if(update && !is.null(oldmut <- mutmod(x, i))) {
      oldpar = pedmut::getParams(oldmut)
      oldparList = lapply(oldpar, function(v) list(female = v[1], male = v[2]))
      useopts = modifyList(oldparList, opts)
    }

    fr = afreq(x, i)
    args = c(list(alleles = names(fr), afreq = fr), useopts)

    mutmod(x, i) = do.call(pedmut::mutationModel, args)
  }

  x
}
