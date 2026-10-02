# ANCOMBC2 differential abundance bar/heatmap plot

Runs ANCOMBC2 on a counts table and metadata, then visualizes
differentially abundant taxa. For two-group comparisons a bar plot of
log-fold changes is returned; for three or more groups a heatmap is
returned; for a continuous `group_col` (e.g. `"dist_km"`) a bar plot of
the effect size per unit increase is returned instead, colored by the
direction of the effect.

## Usage

``` r
ancombc_plot(
  table,
  metadata,
  group_col,
  level = "genus",
  min_prevalence = 0.1,
  p_adjust_method = "holm",
  formula = NULL,
  rand_formula = NULL,
  ref_level = NULL,
  diverging_palette = "BuOr",
  x_axis_title = NULL,
  y_axis_title = NULL,
  bar_colors = c("#56B4E9", "#E69F00"),
  save_table = FALSE,
  table_filename = "ancombc_results.txt",
  ...
)
```

## Arguments

- table:

  Data frame with taxa as rows and samples as columns. The last column
  must contain taxonomy strings (named "taxonomy", "Taxonomy", "taxon",
  "taxa", "Taxa", or "Taxon").

- metadata:

  Data frame with samples as rows. The first column must contain sample
  IDs that match the column names of `table`.

- group_col:

  Character. Name of the column in `metadata` that defines the grouping
  variable. If this column is numeric (a continuous variable), it's
  treated as a covariate instead of a group: ANCOMBC2's
  `group`/structural-zero machinery (which requires discrete groups) is
  disabled, and the resulting plot shows the effect size per unit
  increase rather than a group-vs-group comparison.

- level:

  Character or `NULL`. Taxonomic level to agglomerate to before running
  `ancombc2` (e.g. `"genus"`, `"family"`; case-insensitive). Default
  `"genus"`. Pass `NULL` to skip agglomeration and run ANCOMBC2 directly
  on the ASV/OTU-level table (rows of `table`, as-is).

- min_prevalence:

  Numeric. Prevalence cut-off passed to `ancombc2` (default `0.1`).
  Lower values retain more taxa.

- p_adjust_method:

  Character. Multiple-testing correction method passed to `ancombc2`
  (default `"holm"`). Use `"BH"` for a less strict correction when
  sample sizes are small.

- formula:

  Character. Right-hand side of the fixed-effects formula passed to
  `ancombc2` (e.g. `"group + age"`). If `NULL` (default) `group_col` is
  used alone.

- rand_formula:

  Character. Random-effects formula passed to `ancombc2` for mixed
  models (e.g. `"~ 1 | subject_id"`). Default `NULL`.

- ref_level:

  Character. Reference level for `group_col`. If `NULL` (default) the
  first factor level is used as reference. Use this to change which
  group appears as the baseline in comparisons (e.g. `ref_level = "P2"`
  to compare all other groups against P2).

- diverging_palette:

  Character. Name of the colorblind-friendly diverging palette used for
  the 3+-group heatmap's log-fold-change fill scale. One of `"BuOr"`
  (blue-orange, default), `"BuVm"` (blue-vermillion), `"BuPk"`
  (blue-pink), or `"GnPk"` (green-pink). Ignored for the 2-group /
  continuous bar plot, which uses `bar_colors` instead.

- x_axis_title, y_axis_title:

  Titles for the x- and y-axis. Default `NULL`: the log fold change of
  the comparison on x of the bar plots, and no title otherwise.

- bar_colors:

  Character vector of (at least) 2 colors used for the bar plot (2-group
  or continuous `group_col`). First color is the "positive" direction
  (the non-reference group / increases with the variable); second color
  is the "negative" direction (the reference group / decreases with the
  variable). Default `c("#56B4E9", "#E69F00")` (the same
  colorblind-friendly blue/orange pairing used as the 2-group default
  throughout the package). Ignored for the 3+-group heatmap, which uses
  `diverging_palette` instead.

- save_table:

  Logical. If `TRUE`, saves the full ANCOMBC2 results table to disk.
  Default `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"ancombc_results.txt"`.

- ...:

  Old names of renamed arguments (`col_cond`, `tax_level`, `prv_cut`,
  `p_adj_method`), still accepted with a warning. Any other extra
  argument is an error.

## Value

A `ggplot2` object: a bar plot (2 groups) or a heatmap (\\\geq\\3
groups).

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

# Not run automatically because ANCOMBC2's internal bias-correction step
# can intermittently error on some random bootstrap draws (a known
# ANCOMBC2 edge case, not specific to this dataset), which would make an
# always-run example a flaky check. p_adjust_method = "BH" is less strict
# than the "holm" default; min_prevalence is raised above the 0.1 default to
# filter out rare/sparse taxa before testing.
# \donttest{
ancombc_plot(
  table        = table,
  metadata     = metadata,
  group_col     = "Location",
  min_prevalence      = 0.3,
  p_adjust_method = "BH"
)
#> Registered S3 method overwritten by 'lme4':
#>   method           from
#>   na.action.merMod car 
#> Checking the input data type ...
#> The input data is of type: phyloseq
#> PASS
#> Checking the sample metadata ...
#> The specified variables in the formula: Location
#> The available variables in the sample metadata: BarcodeSequence, LinkerPrimerSequence, Loc, FW, Location, Root_L, Stem_L, Root_DW, Treatment, TOC, Root_FW, ReversePrimer, Leaves_No, SMC, pH, Silt, WHC, TN, Plot, EC, DW, Sand, Clay, Arbus_per, Cod_No, Code, Month, Description
#> PASS
#> Checking other arguments ...
#> The number of groups of interest is: 2
#> Warning: The group variable has < 3 categories 
#> The multi-group comparisons (global/pairwise/dunnet/trend) will be deactivated
#> The sample size per group is: Rhizosphere = 24, Roots = 22
#> PASS
#> Obtaining initial estimates ...
#> Estimating sample-specific biases ...
#> Loading required package: foreach
#> Loading required package: rngtools
#> ANCOM-BC2 primary results ...
#> Conducting sensitivity analysis for pseudo-count addition to 0s ...
#> For taxa that are significant but do not pass the sensitivity analysis,
#> they are marked in the 'passed_ss' column and will be treated as non-significant in the 'diff_robust' column.
#> For detailed instructions on performing sensitivity analysis, please refer to the package vignette.
#> Term 'Location': 1 significant taxa

# }
```
