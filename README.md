
[![R build status](https://github.com/federicomarini/ideal/workflows/R-CMD-check/badge.svg)](https://github.com/federicomarini/ideal/actions)
[![codecov.io](https://codecov.io/github/federicomarini/ideal/coverage.svg?branch=master)](https://codecov.io/github/federicomarini/ideal?branch=master)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)

<img src="man/figures/ideal.png" align="right" alt="" width="120" />

# `ideal` - Interactive Differential Expression AnaLysis in RNA-seq data 

<a href="https://doi.org/10.1186/s12859-020-03819-5"><img src="https://img.shields.io/badge/doi-ideal-blue.svg"><a>
<a href="https://doi.org/10.1002/cpz1.411"><img src="https://img.shields.io/badge/doi-ideal_protocol-blue.svg"><a>

`ideal` is a Bioconductor package containing a Shiny application for analyzing RNA-Seq data in the context of differential expression. 
This enables an interactive and at the same time reproducible analysis, keeping the functionality accessible, and yet providing a comprehensive selection of graphs and tables to mine the dataset at hand.

`ideal` is an R package which fully leverages the infrastructure of the Bioconductor project in order to deliver an interactive yet reproducible analysis for the detection of differentially expressed genes in RNA-Seq datasets. 
Graphs, tables, and interactive HTML reports can be readily exported and shared across collaborators. 
The dynamic user interface displays a broad level of content and information, subdivided by thematic tasks. 
All in all, it aims to enforce a proper analysis, by reaching out both life scientists and experienced bioinformaticians, and also fosters the communication between the two sides, offering robust statistical methods and high standard of accessible documentation.

It is structured in a similar way to the `pcaExplorer`, also designed  as an interactive companion tool for RNA-seq analysis focused rather on the exploratory data analysis e.g. using principal components analysis as a main tool.

The interactive/reactive design of the app, with a dynamically generated user interface makes it easy and immediate to apply the gold standard methods in a way that is information-rich and accessible also to the bench biologist, while also providing additional insight also for the experienced data analyst. 
Reproducibility is supported via state saving and automated report generation.

## Installation

`ideal` can be easily installed using `BiocManager::install()`:

``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ideal")
```

Note that this should be the preferred way to install the latest stable release version.

To also install the packages listed in the `Suggests:` field, you can run

``` r
BiocManager::install("ideal", dependencies = TRUE)
```

to make sure to have for example the required demo dataset (`airway`) when running the app - or if you want to follow through the vignette entirely.

Optionally, if you want to install the development version from GitHub, you can use:

``` r
BiocManager::install("federicomarini/ideal", dependencies = TRUE)
# or alternatively...
devtools::install_github("federicomarini/ideal", dependencies = TRUE)
```

Setting `dependencies = TRUE` should ensure that all packages, including the ones in the `Suggests:` field of the `DESCRIPTION`, are installed - this can be essential if you want to reproduce the code in the vignette, for example.

### Installation troubleshooting

If using `devtools` or `remotes` to install packages, you could run into the warning

``` r
# ... after launching the install_github command
Error: (converted from warning) package ´IRanges´ was built under R version 3.6.2
Execution halted
ERROR: lazy loading failed for package ´ideal´
*removing ´Library/Frameworks/R.framework/Versions/3.6/Resources/library/ideal´
Error: Failed to install 'ideal' from GitHub:
  (converted from warning) installation of package ´....../ideal_1.11.2.tar.gz´ had non zero exit status
```

In this case, you can follow the instructions found at https://remotes.r-lib.org/index.html#environment-variables, which specifically suggest to set `R_REMOTES_NO_ERRORS_FROM_WARNINGS` to `true`. You can do so directly in R before installing the package by entering

``` r
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
# and then again
devtools::install_github("federicomarini/ideal", dependencies = TRUE)
```

### Which version should I use?

If you are a **regular user**, you should **install the latest stable release version**. 
This can be done at best by using `BiocManager::install("ideal")`, as recommended in https://www.bioconductor.org/install/#troubleshoot-bioconductor-packages.
Please follow the general instructions in https://www.bioconductor.org/install to make sure you are using the correct version, matched to the version of the R software in use.

If you are a **software developer** and want to have access to the latest features that are currently in the **devel branch of Bioconductor** (i.e. experimental functionality, and more), you can do so by calling first `BiocManager::install(version = "devel")` as specified in https://bioconductor.org/developers/how-to/useDevel/, then followed by `BiocManager::install("ideal")`.
Keep in mind that according to the release cycle you might need to install the devel version of R itself.

If you just want to use the **bleeding edge version**, which is the one you can find on GitHub, you can install that by calling `BiocManager::install("federicomarini/ideal")` (which is basically a wrapper around `remotes::install("federicomarini/ideal")`).
This approach might be recommended for **experienced users** - based on which Bioconductor version you might be using, you might encounter mismatches in the dependencies if you mix up versions from release and devel branches.

## Quick start

This command loads the `ideal` package

``` r
library("ideal")
```

The main parameters for `ideal` are

- `dds_obj` - a `DESeqDataSet` object. If not provided, then a `countmatrix` and a 
`expdesign` need to be provided. If none of the above is provided, it is possible
to upload the data during the execution of the Shiny App
- `res_obj` -  a `DESeqResults` object. If not provided, it can be computed during
the execution of the application
- `annotation_obj` - a `data.frame` object, with row.names as gene identifiers 
(e.g. ENSEMBL ids) and a column, `gene_name`, containing e.g. HGNC-based gene
symbols. If not provided, it can be constructed during the execution via the 
`org.eg.XX.db` packages
- `countmatrix` - a count matrix, with genes as rows and samples as columns.
If not provided, it is possible to upload the data during the execution of
the Shiny App
- `expdesign` -a `data.frame` containing the info on the experimental covariates
of each sample. If not provided, it is possible to upload the data during the
execution of the Shiny App

The `ideal` app can be launched in different modes:

- `ideal(dds_obj = dds, res_obj = res, annotation_obj = anno)`, where the objects 
are precomputed in the current session and provided as parameters
- `ideal(dds_obj = dds)`, as in the command above, but where the result object is
assembled at runtime 
- `ideal(countmatrix = countmatrix, expdesign = expdesign)`, where instead of 
passing the defined `DESeqDataSet` object, its components are given, namely the 
count matrix (e.g. generated after a run of featureCounts or HTSeq-count) and a 
data frame with the experimental covariates. The design formula can be constructed
interactively at runtime
- `ideal()`, where the count matrix and experimental design can simply be uploaded
at runtime, where all the derived objects can be extracted and computed live. These 
files have to be formatted as tabular text files, and a function in the package 
tries to guess the separator, based on heuristics of occurrencies per line of 
commonly used characters

## Accessing the public instance of `ideal` 

To use `ideal` without installing any additional software, you can 
access the public instance of the Shiny Server made available at the Institute of 
Medical Biostatistics, Epidemiology and Informatics (IMBEI) in Mainz.

This resource is accessible at this address: 

http://shiny.imbei.uni-mainz.de:3838/ideal

## Deploying to a Shiny Server

A deployment-oriented version of the package is available at 
https://github.com/federicomarini/ideal_serveredition. This repository contains also
detailed instruction to setup the running instance of a Shiny Server, where `ideal` 
can be run without further installation for the end-users.

## Code of Conduct

Please note that the ideal project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

## Contact

For additional details regarding the functions of **ideal**, please consult the documentation or 
write an email to marinif@uni-mainz.de. 

### Bug reports/Issues/New features

Please use https://github.com/federicomarini/ideal/issues for reporting bugs, issues or for 
suggesting new features to be implemented.
