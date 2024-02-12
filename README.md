
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MicroBioMeta

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/Steph0522/MicroBioMeta)](https://github.com/Steph0522/MicroBioMeta/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/Steph0522/MicroBioMeta)](https://github.com/Steph0522/MicroBioMeta/pulls)
<!-- badges: end -->

The goal of `MicroBioMeta` is to analyze, visualize, and explore data from microbiome studies of amplicon and metagenomic sequencing at taxonomic, functional, and phylogenetic levels.

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


## Citation

Below is the citation output from using `citation('MicroBioMeta')` in R.
Please run this yourself to check for any updates on how to cite
**MicroBioMeta**.

``` r
print(citation('MicroBioMeta'), bibtex = TRUE)
#> To cite package 'MicroBioMeta' in publications use:
#> 
#>   Hereira S (2024). _MicroBioMeta: What the Package Does (One Line,
#>   Title Case)_. R package version 0.99.0,
#>   <https://github.com/Steph0522/MicroBioMeta>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {MicroBioMeta: What the Package Does (One Line, Title Case)},
#>     author = {Stephanie Hereira},
#>     year = {2024},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/Steph0522/MicroBioMeta},
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
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
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
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
