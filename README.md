
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pedtools <img src="man/figures/logo.png" align="right" height=140/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/pedtools)](https://CRAN.R-project.org/package=pedtools)
[![](https://cranlogs.r-pkg.org/badges/grand-total/pedtools?color=yellow)](https://cran.r-project.org/package=pedtools)
[![](https://cranlogs.r-pkg.org/badges/last-month/pedtools?color=yellow)](https://cran.r-project.org/package=pedtools)
<!-- badges: end -->

## Introduction

The goal of **pedtools** is to provide a lightweight, but comprehensive
tool set for creating, manipulating and visualizing pedigrees with or
without marker data. Common pedigree structures are quickly produced
with tailor-made functions, while a range of utilities enable
modifications like adding or removing individuals, extracting subsets,
loop breaking, and merging pedigrees. The plotting functionality is
imported from the
[kinship2](https://CRAN.R-project.org/package=kinship2) package.

**pedtools** is a continuation of the
[paramlink](https://CRAN.R-project.org/package=paramlink) package, which
is no longer actively developed.

## Installation

To get **pedtools**, install from CRAN as follows:

``` r
install.packages("pedtools")
```

Alternatively, you can obtain the latest development version from
GitHub:

``` r
# install.packages("devtools") # install devtools if needed
devtools::install_github("magnusdv/pedtools")
```

## Example

The following example illustrates how pedigrees and markers may be built
from scratch.

``` r
library(pedtools)

# Create pedigree
x = cousinPed(degree = 0, removal = 2)
x = addChildren(x, father = 3, nch = 2, sex = 2)

# Relabel according to plot order
x = relabel(x, "asPlot")

# Create marker and attach to pedigree
m = marker(x, "7" = "a/b", "11" = "b/b")

# Plot pedigree with genotypes
plot(x, marker = m, hatched = leaves(x))
```

<img src="man/figures/README-example-1.png" width="40%" />

For details about what **pedtools** can do, and many other examples, the
vignette is the recommended place to start.
