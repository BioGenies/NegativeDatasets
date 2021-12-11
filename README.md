# NegativeDatasets

This repository contains the data and code necessary to reproduce the results from the preprint: Katarzyna Sidorczuk, Przemysław Gagat, Filip Pietluch, Jakub Kała, Dominik Rafacz, Laura Bąkała, Jadwiga Słowik, Rafał Kolenda, Stefan Rödiger, Legana C H W Fingerhut, Ira R Cooke, Paweł Mackiewicz, Michał Burdukiewicz *The impact of negative data sampling on antimicrobial peptide prediction*. 


## Getting started

This repository uses [renv](https://CRAN.R-project.org/package=renv) and [targets](https://CRAN.R-project.org/package=targets) packages to control the workflow and assure the reproducibility. 

Some of the data files are too large to store them on GitHub but they can be downloaded using the links below:

- [UniProt data](https://www.dropbox.com/sh/n7hcu1byp1izuwv/AAB6irXnv8S5dE-LEW4QkM-ya?dl=0) - Data directory with reviewed sequences and their annotation downloaded from UniProtKB release 2020_06. These sequences and their annotations were used to create negative data sets in our study. 

- [Prediction results for architectures](https://www.dropbox.com/sh/iuytufcl92kd61a/AAArrO0P9XhZavDxfTpqjIhua?dl=0) - Results directory with prediction results for all 660 models trained and tested in our study. These files are necessary for calculation of models' performance and generation of plots and tables from the paper.

To reproduce the results clone the repo, set your path to the directories with data files and:

``` r
renv::restore()
targets::tar_make()
```


## Content

**\_targets.R** - reproducible pipeline for generation of all data sets and results processing, 

**data** - data files used during the study, e.g. for creation of the positive dataset,

**drafts** - draft codes used for initial exploratory analyses,

**functions** - all functions used for running the pipeline and obtaining results,

**presentations** - presentation files for this project,

**renv** - renv package files,

**reports** - reports with initial analyses,

**third-party** - third-party executables used in the pipeline.
