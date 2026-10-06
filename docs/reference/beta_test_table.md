# Table of Permanova or Betadisper

This function create a table with the results of permanova or betadisper

## Usage

``` r
beta_test_table(
  table,
  metadata,
  formula_str,
  distance = "euclidean",
  test = c("permanova", "betadisper"),
  permutations = 999,
  mc_samples = 1,
  strata_var = NULL,
  digits = 3,
  save_table = FALSE,
  table_filename = "beta_test_results.txt"
)
```

## Arguments

- table:

  A precomputed distance matrix (`matrix` or `dist` object, e.g. from
  [`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)),
  or a data frame with taxonomy, where the columns are the samples and
  rows are ASVs or taxa.

- metadata:

  Data frame of characteristics or important information of the samples

- formula_str:

  Model formula

- distance:

  Method for calculating pairwise distances. Same convention as
  `beta_div_plot`'s `distance` argument: `"compositional"` runs ALDEx2's
  CLR transform
  ([`ALDEx2::aldex.clr`](https://rdrr.io/pkg/ALDEx2/man/aldex.clr.function.html))
  on the raw counts and then a Euclidean distance on the CLR values
  (requires `table` to be a raw abundance data frame with a taxonomy
  column, not a precomputed distance matrix);
  `"aitchison"`/`"robust.aitchison"` use vegan's built-in Aitchison
  distance
  ([`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  with a pseudocount); any other value (e.g. `"euclidean"`, `"bray"`) is
  passed straight to
  [`vegan::vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  on the raw values (`"euclidean"` default; case-insensitive).

- test:

  Statistical test to run: one of `"permanova"` (default) or
  `"betadisper"`. Case-insensitive.

- permutations:

  Number of permutations required

- mc_samples:

  Number of ALDEx2 Monte Carlo instances used when
  `distance = "compositional"`. With `1` (default) the clr values of one
  random instance are used: fast, but the result changes a little
  between runs (use [`set.seed()`](https://rdrr.io/r/base/Random.html)).
  With more, the clr values are averaged across instances, which gives
  an almost identical result in every run; `128` (ALDEx2's default) is
  suggested for final analyses, and takes longer. Ignored for other
  distances.

- strata_var:

  Group or variable within which permutations are restricted

- digits:

  Number of decimal places for the numeric columns (except the p-value).

- save_table:

  Logical. If `TRUE`, saves the results table to disk. Default `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"beta_test_results.txt"`.

## Value

A data frame (class `mbm_test_table`) with the test results: one row per
term and the columns returned by
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
(`Df`, `SumOfSqs`, `R2`, `F`, `Pr(>F)`) or
[`vegan::permutest()`](https://vegandevs.github.io/vegan/reference/anova.cca.html),
plus `Term`. Printing it (e.g. typing its name) draws the formatted
table figure;
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
returns that figure as a `ggplot` object, e.g. to combine it with other
plots
([`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html),
`patchwork`) or save it with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

The first column of `metadata` must hold the sample IDs; metadata rows
are matched to the samples by ID, so their order doesn't matter. A
precomputed distance is used as-is (`distance` is ignored). With
`distance = "compositional"` and `mc_samples = 1`, the clr values come
from one random Monte Carlo instance of
[`ALDEx2::aldex.clr()`](https://rdrr.io/pkg/ALDEx2/man/aldex.clr.function.html);
call [`set.seed()`](https://rdrr.io/r/base/Random.html) before the
function to make the result reproducible, or use `mc_samples = 128`.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

# Example using a data frame
beta_test_table(
  table       = table,
  metadata    = metadata,
  formula_str = "Location*Treatment",
  distance    = "bray",
  test        = "permanova",
  permutations = 999,
  strata_var  = "Plot"
)


# Example using a distance matrix
dist_matrix <- vegan::vegdist(
  t(table[, setdiff(colnames(table), "taxonomy")]),
  method = "bray"
)
beta_test_table(
  table       = dist_matrix,
  metadata    = metadata,
  formula_str = "Location",
  test        = "betadisper"
)


# Compositional PERMANOVA (CLR via ALDEx2, then Euclidean)
beta_test_table(
  table       = table,
  metadata    = metadata,
  formula_str = "Location",
  distance    = "compositional",
  test        = "permanova",
  permutations = 999
)
#> no conditions provided: forcing denom = 'all'
#> no conditions provided: forcing conds = 'NA'
#> conditions vector supplied
#> operating in serial mode
#> computing center with all features
```
