# 🧬 MicroBioMeta <img src="man/figures/logo.jpg" alt="Logo" width="150" align="right"/>

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/Steph0522/MicroBioMeta)](https://github.com/Steph0522/MicroBioMeta/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/Steph0522/MicroBioMeta)](https://github.com/Steph0522/MicroBioMeta/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

The goal of `MicroBioMeta` is to support beginners in microbiome data
analysis, whether working with metagenomic or metabarcoding datasets. It
enables users with limited experience in R to explore, visualize, and
analyze results with minimal coding and customization. This helps save
time on scripting, while still producing statistically robust results
and publication-ready figures in an intuitive and accessible way

## 👩‍💻 Contributors

- Stephanie Hereira - lead developer.
- Nina Montoya - co-developer and functions.
- Karla Zarco - co-developer and functions.

## :computer: Installation instructions

```r
install.packages("devtools")
library(devtools)
install_github("Steph0522/MicroBioMeta")
library(MicroBioMeta)

or with BiocManager:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("Steph0522/MicroBioMeta")
```

## :green_book: Examples

For examples and uses go to the [web page](https://steph0522.github.io/MicroBioMeta/)


## :raised_hands: Code of Conduct

Please note that the `MicroBioMeta` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## :hammer: Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.18/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://Steph0522.github.io/MicroBioMeta)
  is automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.18/biocthis)*.
