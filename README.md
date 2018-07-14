<!-- README.md is generated from README.Rmd. Please edit that file -->
pedtools
========

Introduction
------------

The goal of pedtools is to provide a lightweight, but comprehensive tool set for creating, manipulating and visualizing pedigrees with or without marker data. Common pedigree structures are quickly produced with tailor-made functions, while a range of utilities enable modifications like adding or removing individuals, extracting subsets, loop breaking, and merging pedigrees. The plotting functionality is imported from the [kinship2](https://CRAN.R-project.org/package=kinship2) package.

pedtools is a rewritten and improved version of one part of the many-faceted package [paramlink](https://CRAN.R-project.org/package=paramlink), which is now heading towards retirement. The functionality of paramlink dealing with basic manipulation of pedigrees and markers is moved to pedtools, in most cases using the same function names. Under the hood, however, many things are changed, including a new S3 class `ped` for basic pedigree objects.

Installation
------------

pedtools is under development and not on CRAN yet. However, you can install the latest version from GitHub as follows:

``` r
 # First install devtools if needed
if(!require(devtools)) install.packages("devtools")

# Install pedtools from github
devtools::install_github("magnusdv/pedtools")
```

Example
-------

The following shows how to create a pedigree with a consanguineous mating between first cousins. We attach a SNP marker for which the child is homozygous for the 'a' allele.

``` r
library(pedtools)

x = cousinsPed(degree=1, child=TRUE)
m = marker(x, '9'='a', alleles=c('a','b'))
plot(x, m)
```

![](man/figures/README-example-1.png)
