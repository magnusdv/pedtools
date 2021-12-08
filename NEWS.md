# pedtools 1.1.0

The main theme of this version is to make `pedtools` more adapted to piping, e.g., allowing chains of commands like `nuclearPed() |> addSon(1) |> addMarker(alleles = 1:2)`.

## New features

* New functions `setAfreq()`, `setChrom()`, `setGenotype()`, `setMarkername()`, `setPosition()` for modifying marker attributes. These are alternatives to the previous in-place modifiers `afreq<-()` a.s.o..

* New function `addMarker()` which simplifies that common task of creating and attaching a single marker. The command `addMarker(x, ...)` is equivalent to `addMarkers(x, marker(x, ...))`.

* The new `addMarker()` accepts ped lists, so that one can write e.g. `list(singleton(1), singleton(2)) |> addMarker("1" = "a/b", alleles = c("a", "b"))`

* `readPed()` gains the argument `colSep`, which fixes the previous inability to handle names with spaces.

* New function `descentPaths()`, mostly intended for use in other ped suite packages.

* `relabel(x, new = "generations")` now gives automatic, generation-aware labelling: I-1, I-2, II-1, ...

* `generations()` gains argument `maxOnly`, by default TRUE. If FALSE, the function returns the generation number of each individual.


# pedtools 1.0.1

## New features

* New function `generations()` for counting generations in pedigrees.

* New function `newMarker()` (mostly for internal use).

* `plot.ped()` gains a new parameter `twins`.

* `father()` and `mother()` now accepts ped lists as input.

* Added info and links to **ped suite** in README.

## Bug fixes

* Fixed bug in `getGenotypes()` affecting pedigrees with numerical labels.

* Fixed bug in `doubleCousins()`.


# pedtools 0.9.7

## Breaking changes

* The rarely-used function `cousins()` (not to be confused with `cousinPed()`) is temporarily retracted, since it did not work as intended.

## New features

* New constructor `newPed()` (mainly for internal use).

* New function `foundersFirst()`, moved from the **ribd** package.

* In `addChildren()`, unspecified `nch` is now allowed, and defaults to `length(ids)` or `length(sex)`.

* `transferMarkers()` has a new argument `checkSex()`, and has been made more efficient by skipping redundant validation steps. 

* The functions `swapSex()`, `alleles()` and `internalID()` now work for lists of pedigrees.

* `getComponent()` gained a new argument `errorIfUnknown()`.


## Bug fixes

* `unrelated()` and `siblings()` have been improved and cleaned of bugs.

* Fixed an obscure bug in `plot.singleton()`.


# pedtools 0.9.6

## Breaking changes

* `getMap(na.action = 1)` is re-implemented and now behaves slightly differently. (This was necessary to improve the handling of linked markers in `pedprobr::merlin()`.)

* The order of individuals in `linearPed()` now always follows the "asPlot" pattern, as for the other basic pedigrees. (Missed this in the previous version.)

## New features

* `plot.ped()` gains arguments `textInside`, `textAbove` and `carrier`.

* `transferMarkers()` has new arguments `fromIds` and `toIds` enabling transfer between differently-named individuals.

* In `setMarkers()` and friends, the shortcut `locusAttributes = "snp-12"` may be used to indicate that all supplied markers are SNPs with alleles 1 and 2. Further shortcuts are "snp-ab" and "snp-AB".

* `setMap()` is extended to ped lists.


# pedtools 0.9.5

## Breaking changes
* Built-in pedigree structures are now labelled according to default plotting order. In particular, this means that pedigrees made by `halfSibPed()`, `cousinPed()` and `halfCousinPed()` are ordered differently than before.

* In `plot.ped()`,  the parameter `skipEmptyGenotypes` is replaced by `showEmpty`, with default value `FALSE`.

* Function `xxxFrequencyDatabase()` have been renamed to `xxxFreqDatabase()`

* The marker attribute `posCm` has been removed, to avoid confusion with the physical position.

* `marker()` now checks for duplicated allele labels.

* `setMarkers()` now checks for duplicated marker names (and allele labels, through `marker()`; see previous point).

## New features
* `readPed()` and friends now automatically recognises allele separator "/" when genotypes are written like "a/b". Other separators must be indicated with `sep` as before, e.g., `readPed(..., sep = ",")`.

* New function `getGenotypes()`, which is similar to `getAlleles()`, but returns a matrix of genotypes written as "a/b".

* More flexible conversion of pedigrees to data frames, with new arguments `sep` and `missing` in `as.data.frame.ped()`.

* New function `setMap()`, facilitating setting chromosome and physical position attributes.

* `marker()` has a new argument `geno`, allowing commands like `marker(nuclearPed(1), geno = c("a/a", NA, "a/b"))`.

* `print.marker()` has been overhauled and gives a more coherent output.

* `halfSibPed()` has a new argument `type`, either "paternal" (default) or "maternal".

* `reorderPed()` by default orders by numerical value, if all labels are numeric.

* `plot.ped()` has a new argument `hint`, which is forwarded to `kinship2::plot.pedigree()`. This is necessary in some cases where the automatic plotting fails to give a nice pedigree. An example is given in `?plot.ped`.

* `plot.ped()` gains argument `hatched`, which will eventually replace `shaded`.

* Added default values allows executing `singleton()` and `nuclearPed()` with no input.

* Parts of `plotPedList()` have been restructured. In particular, the new argument `groups` makes it easier to control grouping and frames. Previous argument `frametitles` has been renamed to `titles`, because it also works without frames.



# pedtools 0.9.4

## Breaking changes

* The `plot.ped()` argument `id.labels` is now deprecated in favour of the new `labs`. This works *almost* as before, with some exceptions documented here. The `labs` argument should be thought of as *who should be labelled* rather than *what are the labels*. For example, with `x = singleton(1)`, the previous `plot(x, id.labels = "2")` would rename the singleton to "2". In contrast, `plot(x, labs = "2")` will not show any label (since `x` doesn't have a member named "2"). In general `intersect(labs, labels(x))` determines who gets a label.  

  Another change is that if `labs` is a function, it is now applied to the pedigree `x`, not to `labels(x)`. This makes it very easy to apply standard pedigree functions like `females()`, `nonfounders()` and `typedMembers()`, since they can be referred to simply by name: `plot(x, labs = females)`.

* The implementation of `doubleCousins()` is improved, and some edge cases smoothed out, but the final ordering of individuals may be different in some cases now.

* `writePed()` has been partially rewritten, to make it more similar to `readPed()`. By default, only the "ped" file is written. New logical arguments "famid" and "header" provide further control of this file. 

  Writing files in merlin format (indicated by `merlin = TRUE`) is internally now done in a separate function. This option is rarely needed by end users, but is called by e.g. `pedprobr::likelihoodMerlin()`.

## New features

* Genotype assignment in `marker()` is more user-friendly now, allowing inputs like `marker(singleton("s"), s = "A/B")`. Previously, heterozygous genotypes had to be provided allele-wise, e.g., `marker(singleton("s"), s = c("A", "B"))`. The character "/" must be used as allele separator and will always be interpreted as such.  
Given the simplicity of the new syntax I recommend that homozygous genotypes are also written out fully, e.g. `s = "B/B"` instead of the previous (but still functional) `s = "B"`.

* New functions `commonAncestors()` and `commonDescendants()` for finding common ancestors/descendants of members in a pedigree.

* The functions `ancestors()` and `descendants()` have a new logical argument, `inclusive`, indicating if the person itself should be included.

* New function `setSex()`. This is inverse to `getSex()` in the sense that `setSex(x, sex = getSex(x, named = T))` is identical to `x`, whether `x` is a single `ped` object or a list of such (with unique ID labels).  
The old `swapSex()` is often more convenient in practise, since it automatically deals with spouses. One situation where `setSex()` is the only option, is when one wants to assign unknown sex (`sex = 0`) to someone.

* New function `setMap()`, which can be used for assigning chromosome and position attributes to marker objects.


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
