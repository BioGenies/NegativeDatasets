
# NegativeDatasets

# Citation

Katarzyna Sidorczuk, Przemysław Gagat, Filip Pietluch, Jakub Kała,
Dominik Rafacz, Laura Bąkała, Jadwiga Słowik, Rafał Kolenda, Stefan
Rödiger, Legana C H W Fingerhut, Ira R Cooke, Paweł Mackiewicz, Michał
Burdukiewicz [*Benchmarks in antimicrobial peptide prediction are biased
due to the selection of negative
data.*](https://doi.org/10.1101/2022.05.30.493946)

## Getting started

This repository contains the data and code necessary to reproduce the
results from the paper *Benchmarks in antimicrobial peptide prediction
are biased due to the selection of negative data*. It uses
[renv](https://CRAN.R-project.org/package=renv) and
[targets](https://CRAN.R-project.org/package=targets) packages to
control the workflow and assure the reproducibility.

Some of the data files are too large to store them on GitHub but they
can be downloaded using the links below:

-   [UniProt
    data](https://www.dropbox.com/sh/n7hcu1byp1izuwv/AAB6irXnv8S5dE-LEW4QkM-ya?dl=0) -
    Data directory with reviewed sequences and their annotation
    downloaded from UniProtKB release 2020_06. These sequences and their
    annotations were used to create negative data sets in our study.

-   [Prediction results for
    architectures](https://www.dropbox.com/sh/iuytufcl92kd61a/AAArrO0P9XhZavDxfTpqjIhua?dl=0) -
    Results directory with prediction results for all 660 models trained
    and tested in our study. These files are necessary for calculation
    of models’ performance and generation of plots and tables from the
    paper.

To reproduce the results clone the repo, set your path to the
directories with data files and:

``` r
renv::restore()
targets::tar_make()
```

## Content

**\_targets.R** - reproducible pipeline for generation of all data sets
and results processing,

**data** - data files used during the study, e.g. for creation of the
positive dataset,

**drafts** - draft codes used for initial exploratory analyses,

**functions** - all functions used for running the pipeline and
obtaining results,

**presentations** - presentation files for this project,

**renv** - renv package files,

**reports** - reports with initial analyses,

**third-party** - third-party executables used in the pipeline.

# Important links

-   <https://github.com/BioGenies/NegativeDatasets>: the repository
    containing the code necessary to reproduce results of our analysis.
-   <https://github.com/BioGenies/NegativeDatasetsArchitectures>: the
    repository containing all architectures considered in our analysis.
-   <https://github.com/BioGenies/AMPBenchmark>: the source code of
    AMPBenchmark.

# Contact

If you have any questions, suggestions or comments, contact [Michal
Burdukiewicz](mailto:michalburdukiewicz@gmail.com).
