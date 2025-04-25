
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MicroBioMeta

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/Steph0522/MicroBioMeta)](https://github.com/Steph0522/MicroBioMeta/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/Steph0522/MicroBioMeta)](https://github.com/Steph0522/MicroBioMeta/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/MicroBioMeta.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MicroBioMeta)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/MicroBioMeta.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/MicroBioMeta)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/MicroBioMeta.svg)](http://bioconductor.org/packages/stats/bioc/MicroBioMeta/)
[![Bioc
support](https://bioconductor.org/shields/posts/MicroBioMeta.svg)](https://support.bioconductor.org/tag/MicroBioMeta)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/MicroBioMeta.svg)](https://bioconductor.org/packages/release/bioc/html/MicroBioMeta.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/MicroBioMeta.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/MicroBioMeta/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/MicroBioMeta.svg)](https://bioconductor.org/packages/release/bioc/html/MicroBioMeta.html#since)
[![R-CMD-check-bioc](https://github.com/Steph0522/MicroBioMeta/actions/workflows/R-CMD-check-bioc.yaml/badge.svg)](https://github.com/Steph0522/MicroBioMeta/actions/workflows/R-CMD-check-bioc.yaml)
<!-- badges: end -->

The goal of `MicroBioMeta` is to support beginners in microbiome data
analysis, whether working with metagenomic or metabarcoding datasets. It
enables users with limited experience in R to explore, visualize, and
analyze results with minimal coding and customization. This helps save
time on scripting, while still producing statistically robust results
and publication-ready figures in an intuitive and accessible way

## Contributors

- Stephanie Hereira - lead developer.
- Nina Montoya - co-developer and functions.
- Karla Zarco - co-developer and functions.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `MicroBioMeta` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("MicroBioMeta")
```

And the development version from
[GitHub](https://github.com/Steph0522/MicroBioMeta) with:

``` r
BiocManager::install("Steph0522/MicroBioMeta")
```

## Examples

This is a basic example which shows you how to solve a common problem:

``` r
library("MicroBioMeta")
## basic example code
```

## Citation

Below is the citation output from using `citation('MicroBioMeta')` in R.
Please run this yourself to check for any updates on how to cite
**MicroBioMeta**.

``` r
print(citation('MicroBioMeta'), bibtex = TRUE)
#> To cite package 'MicroBioMeta' in publications use:
#> 
#>   Steph0522 (2025). _MicrobioMeta_. doi:10.18129/B9.bioc.MicroBioMeta
#>   <https://doi.org/10.18129/B9.bioc.MicroBioMeta>,
#>   https://github.com/Steph0522/MicroBioMeta/MicroBioMeta - R package
#>   version 0.99.0, <http://www.bioconductor.org/packages/MicroBioMeta>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {MicrobioMeta},
#>     author = {{Steph0522}},
#>     year = {2025},
#>     url = {http://www.bioconductor.org/packages/MicroBioMeta},
#>     note = {https://github.com/Steph0522/MicroBioMeta/MicroBioMeta - R package version 0.99.0},
#>     doi = {10.18129/B9.bioc.MicroBioMeta},
#>   }
#> 
#>   Steph0522 (2025). "MicrobioMeta." _bioRxiv_. doi:10.1101/TODO
#>   <https://doi.org/10.1101/TODO>,
#>   <https://www.biorxiv.org/content/10.1101/TODO>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {MicrobioMeta},
#>     author = {{Steph0522}},
#>     year = {2025},
#>     journal = {bioRxiv},
#>     doi = {10.1101/TODO},
#>     url = {https://www.biorxiv.org/content/10.1101/TODO},
#>   }
```

Please note that the `MicroBioMeta` was only made possible thanks to
many other R and bioinformatics software authors, which are cited either
in the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `MicroBioMeta` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

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
