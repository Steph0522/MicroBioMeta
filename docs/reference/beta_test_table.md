# Table of Permanova or Betadisper

This function create a table with the results of permanova or betadisper

## Usage

``` r
beta_test_table(
  table,
  metadata,
  formula_str,
  method = "euclidean",
  test = c("permanova", "betadisper"),
  permutations = 999,
  strata_var = NULL,
  decimales = 3
)
```

## Arguments

- metadata:

  Data frame of characteristics or important information of the samples

- formula_str:

  Model formula

- method:

  Method for calculate pairwise distances of a matrix

- test:

  Permanova or betadisper

- permutations:

  Number of permutations required

- strata_var:

  Group or variable within which permutations are restricted

- decimales:

  Number of decimales required

- matriz:

  Distance matrix or data frame with taxonomy, where, the columns are
  the samples and rows are ASV's or taxa.

## Value

A table with the results of R2, F and p value

## Examples
