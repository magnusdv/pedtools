
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
loop breaking, and merging pedigrees. The plotting feature imports
machinery from the
[kinship2](https://CRAN.R-project.org/package=kinship2) package.

**pedtools** is the hub of the **ped suite**, a collection of R packages
for pedigree analysis, including applications in forensic and medical
genetics. The **ped suite** has its own [GitHub
repository](https://github.com/magnusdv/pedsuite) and a dedicated
[website](https://magnusdv.github.io/pedsuite/) offering more
information.

<span style="color:red;"> Try the **QuickPed** app for building and
analysing pedigrees here: <https://magnusdv.shinyapps.io/quickped/>
</span>

## Installation

To get **pedtools**, install from CRAN as follows:

``` r
install.packages("pedtools")
```

Alternatively, fetch the latest development version from GitHub:

``` r
devtools::install_github("magnusdv/pedtools")
```

## Example

The following example illustrates a step-by-step creation of a pedigree
with a marker object.

``` r
library(pedtools)

# Start with two half brothers
x = halfSibPed(type = "paternal")

# Make 5 female
x = swapSex(x, 5)

# Add a sister to 5 (parents are 2 and 3)
x = addDaughter(x, parents = 2:3)

# Add inbred child
x = addSon(x, parents = 4:5)

# Create marker
x = addMarker(x, "7" = "a/b")

# Plot pedigree with genotypes
plot(x, marker = 1, hatched = 7)
```

<img src="man/figures/README-example-1.png" width="40%" style="display: block; margin: auto;" />

The process of building pedigrees is perfectly suited for the pipe
operator `|>` recently introduced in R. For example, the above pedigree
could have been created as follows:

``` r
x = halfSibPed(type = "paternal") |>
  swapSex(5) |>
  addDaughter(parents = 2:3) |>
  addSon(parents = 4:5) |>
  addMarker("7" = "a/b")
```

For details about what **pedtools** can do, and many other examples,
[the
vignette](https://cran.r-project.org/package=pedtools/vignettes/pedtools.html)
is a good place to start.
