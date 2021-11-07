# NegativeDatasets

Reviewed sequences and associated information were downloaded from UniProtKB release 2020_06.

Link to Dropbox directory with data files: https://www.dropbox.com/sh/n7hcu1byp1izuwv/AAB6irXnv8S5dE-LEW4QkM-ya?dl=0

## Getting started

This repository uses [renv](https://CRAN.R-project.org/package=renv) and [targets](https://CRAN.R-project.org/package=targets) packages to control the workflow and assure the reproducibility. To reproduce the results clone the repo and:

``` r
renv::restore()
targets::tar_make()
```

## Content

**publication_results**: all files and tables for the final publication.
