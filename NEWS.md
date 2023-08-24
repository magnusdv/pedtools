# pedtools 2.3.1

## New features

*	New function `singletons()` for creating a list of singletons.

*	New S3 `summary()` method handling lists of (possibly disconnected) pedigrees.

*	`setGenotype()` can set the genotype of multiple individuals *or* markers in the same call. Also, individuals can be specified with selector functions like `leaves()` and `founders()`, as in: `nuclearPed() |> addMarker() |> setGenotype(ids = founders, geno = "1/2")`.

*	`plot.ped()` gains argument `trimLabs`, defaulting to TRUE, which removes line breaks at the start of ID labels. For example: `nuclearPed() |> addMarker(geno = c("a/a", "b/b", NA)) |> plot(marker = 1, labs = "1")`. (Compare with `plot(..., trimLabs = F)`.)

*	The `showEmpty` argument of `plot.ped()` is more user friendly, allowing for instance: `nuclearPed() |> addMarker() |> plot(marker = 1, showEmpty = "1")`.

*	`readFreqDatabase()` gains an optional argument `df`, a data frame of allele frequencies in either "list" or "allelic ladder" format. This is useful in cases where the raw data must be read or modified manually for some reason.

*	Add CITATION file.

## Bug fixes

*	Forgotten argument `strict` in `setAfreq()` for pedlists.
*	Fix glitch affecting `textAbove` and `textInside`.
*	Fix `title` in (still experimental) `plot.list()`.
* Fixed interpretation of marker names in `setMarkers(x, alleleMatrix)`.


# pedtools 2.2.0

## New features

* The **pedprobr** function `setMutationModel()` has been moved to **pedtools** and renamed to `setMutmod()`. Unlike its predecessor, this has a new argument `update`, allowing to update existing models (i.e., leaving unspecified parameters unchanged) instead of creating new models from scratch.

* In `swapSex()` and `setSex()`, the `ids` parameter may now be the name of a selector function operating on the input pedigree. This is convenient when piping; for example, `x |> setSex(ids = leaves, sex = 0)` sets unknown sex for all leaves of `x`.

* `print.ped()` now returns the `ped` object (not the data frame, as before) invisibly.

* `plot.ped()` gains new arguments `fill`, `lty`, `lwd` and `hatchDensity`.

* `generations()` have been rewritten, with new parameter `what`.

* `addMarker()` gains argument `locusAttr`.
 
## Bug fixes

* Fixed a bug in `randomPed()`, which caused the function to run out of mating options sometimes.

* The pedigree plot alignment fails in some cases (see https://github.com/mayoverse/kinship2/issues/13). When this happens the plot method reverts to DAG mode and gives a warning.


# pedtools 2.1.1

## New features

* The `plot.ped()` method has been internally refactored into 5 functions. Three of these calculate various parameters: `.pedAlignment()`, `.pedScaling()` and `.pedAnnotation()`. The remaining two, `drawPed()` and `.annotatePed()` actually draw stuff on the graphics device. As indicated by the dot prefixes, these functions are primarily intended for internal use. Nevertheless, they are documented and exported, to make them available for other packages requiring special plot methods. (For example, the latest version of `ibdsim2::haploDraw()` use this to compute automatic margins.)

* The function `randomPed()` has been completely rewritten, ensuring that the output is always a connected pedigree. The new version takes as input the total pedigree size (n) and the number of founders (f).

* `removeIndividuals()` gains an argument `remove`, taking as value either "ancestors" or "descendants" (possibly abbreviated. The default value ("descendants") behaves as the previous version. A typical application of `remove = "ancestors"` is to remove founders, as in `linearPed(2) |> removeIndividuals(1, remove = "anc")`. 

* Both `relabel()` and `removeIndividuals()` gain an argument `returnLabs`. If TRUE, the functions return a vector of the pedigree members about to be modified/removed, instead of actually performing any changes.

* Pedigree construction is generally faster, due to various code improvements. Also, `ped()` gains an argument `detectLoops`, which if set to FALSE may cut runtime significantly in some cases. 

* New function `setFounderInbreeding()`, which is more flexible (and pipe friendly!) than the previous `founderInbreeding<-()`. The latter will continue to exist. 

* More informative `summary()` for pedigrees.


# pedtools 2.0.0

This version introduces a number of changes in the pedigree plotting. The alignment of individuals is still done with `kinship2`, but the calculation of scaling, margins and symbol sizes are now done in pedtools. As a result, pedigrees plotted with old code may look slightly different.

## Breaking changes

* Pedigree symbols should always have the same height/width (e.g., perfect squares for males). Previously, symbols were squished in many cases, sometimes heavily so.

* Pedigrees now always span the entire plot region, which was often not the case before.

* Plotting of singletons has been completely rewritten, and is now done in `plot.ped()`. The previous method `plot.singleton()` has been removed. As a result, singleton plots are much more consistent and appear centred in the plot region.

* Some efforts are done to prevent unneeded duplication of founders, for instance in the case of 3/4-siblings: `nuclearPed(2) |> addSon(3) |> addSon(4:5) |> plot()`

* The default plot margins have been set to 1.

* `plotPedList()` is better at guessing relative widths.

* Several minor tweaks in `plotPedList()` in response to changes in `plot.ped()`.

* The function `relabel()` now has "asPlot" as default, i.e., using the numerical plot order.

## New features

* Pedigrees can now be plotted as directed acyclic graphs (DAGs), by adding `arrows = TRUE` in `plot()`.

* Plotting pedigrees with selfing is now supported, by automatically switching to DAG mode.

* The argument `margins` in `plot()` now accepts a single number, making it more user friendly. The default, `margins = 1` corresponds to `par(mar = c(1,1,1,1))`.

## Bug fixes

* Reset graphical parameters after `plotPedList()`

*

# pedtools 1.3.0

## New features

* New function `avuncularPed()` for creating aunt/uncle - nephew/niece pedigrees.

* New function `addAllele()` for extending the allele set of a marker.

* `addSon()` and `addDaughter()` are now more flexible. The previous argument `parent` has been renamed to `parents` and accepts one or two parents in any order.

* `mergePed()` has been overhauled. In particular the new argument `by` makes it much more user friendly.

* `setAfreq()` gains argument `strict`.


## Other

* Minor improvements of README and vignette.

* Fixed bug in `setGenotype()` when setting multiple markers.

* Fixed bug ignoring alleles in `distributeMarkers()`.


# pedtools 1.2.0

## New features

* **pedtools** now depends on R 4.1 (or later) because of the pipe operator `|>`.

* New function `setSNPs()` for creating and attaching a set of SNP markers with given positions and allele frequencies.

* New function `distributeMarkers()` for creating and attaching equal markers evenly across a set of chromosomes (by default, the human autosomes).

* New function `halfSibTriangle()` implementing an interesting breeding pattern.

## Bug fixes

* `transferMarkers()` now ignores members of unknown sex when checking compatibility.

* Fixed bug in `addMarker()` when input is a list of pedigrees.

* Fixed glitches in `setMap()`.

## Other

* Various improvements in code and docs.

* Added many tests.

* Rewrite README example to show piping.


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
