# Development version

## Breaking change

* The `plot.ped()` argument `id.labels` is now deprecated in favor of the new `labs`. This works *almost* as before, with some exceptions documented here. The `labs` argument should be thought of as *who should be labelled* rather than *what are the labels*. For example, in previous versions, `plot(singleton(1), id.labels = "2")` would rename the singleton to "2". In contrast, `plot(singleton(1), labs = "2")` will not show any label. In general `intersect(labs, labels(x))` determines who gets a label.

Another change is that if `labs` is a function, it is now applied to the pedigree `x`, not to `labels(x)`. This makes it very easy to apply standard pedigree functions like `females()`, `nonfounders()` and `typedMembers()`, since they can be referred to simply by name: `plot(x, labs = females)`.

# pedtools 0.9.3

## New features

* New function `readFrequencyDatabase()` reads databases. Both list formats and
allelic ladders are supported.

* Marker attributes "chrom" and "name" are now easier to get/set in ped lists.

* The `relabel()` function now also works for ped lists.

## Bug fixes

* `relabel()` now works correctly in pedigrees with broken loops

* `mendelianCheck()` didn't always print as intended


# pedtools 0.9.2

## New features

* The `labels()` function now also works for ped lists (returning a list of vectors).

## Bug fixes

* The previous version of `getSex()` was buggy; this has been rewritten and made more efficient.


# pedtools 0.9.1

## New features

* New functions for extracting marker properties: `emptyMarkers()` and `nTyped()`. 
These are generic, with methods for `marker`, `ped` and `list`.

* The functions `allowsMutations()`, `isXmarker()` and `nAlleles()` are now generic, 
with methods for `marker`, `ped` and `list`.

* `plot.ped()` now accepts functional forms of the arguments `id.labels`, `shaded` 
and `starred`. This simplifies certain plotting tasks, allowing calls like 
`plot(cousinPed(1), shaded = founders, starred = leaves)`.

* `mutmod<-()` now allows to set the same mutation model for multiple markers 
in one call.

* Many utility functions now operate not only on single pedigrees but also on 
lists of pedigrees. These include `chrom()`, `name()`, `selectMarkers()`, 
`setMarkers()`, `typedMembers()` and `untypedMembers()`, 

* `selectMarkers()` and friends now accepts boolean marker selection, meaning that 
the `markers` argument may be a logical vector (of length equal to the number of 
attached markers).

## Bug fixes

* `readPed()` is now more careful regarding marker names. In particular, it should
now preserve all names exactly as given, and raise an error if encountering duplicated 
names.

# pedtools 0.9.0

* Initial CRAN release.
