## Resubmission
This is a resubmission. In this version I have:

* Removed "Tools for" from the package title

* Reset the par() settings in the example where this was changed.

* Added `on.exit(par(opar))` in functions where par() settings are changed

In one of these functions, `plot.ped()` it is sometimes useful to keep the graphical parameters in order to enable further modifications. To achieve this I have added an argument `keep.par` (default = FALSE) to the function, and inside the function I use 
```
if(!keep.par)
    on.exit(par(opar))
```

## Resubmission
This is a resubmission. In this version I have:

As requested:
* Removed "in R" from the package title.

* Modified the code of `getAlleles()`, so that it always returns
  a character matrix, as documented (instead of NULL in some cases).

And also:
* Added a URL entry to the DESCRIPTION file.

* Made some very minor style changes.

* Changed the name of a single variable throughout (`markerdata` --> `MARKERS`).
  
I have rerun all checks to ensure no errors/warnings/notes were introduced.
  

## Test environments
* local Windows 10 install, R 3.6.1
* devtools::check_win_devel()
* rhub::check_for_cran()

## R CMD check results

0 errors | 0 warnings | 1 note

## Comments

NOTE: This is a new release. 

