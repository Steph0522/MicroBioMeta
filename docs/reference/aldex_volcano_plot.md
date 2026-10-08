# ALDEx2 volcano or effect size plot

Runs ALDEx2 between two groups and plots the results as a volcano plot
or an effect size plot.

## Usage

``` r
aldex_volcano_plot(
  table,
  metadata,
  group_col,
  type = "volcano",
  col_inf = "#56B4E9",
  col_sup = "#E69F00",
  threshold_lower = -1.5,
  threshold_upper = 1.5,
  cond = NULL,
  pval_threshold = 0.05,
  p_adjust_method = "BH",
  show_labels = TRUE,
  taxa = NULL,
  x_axis_title = NULL,
  y_axis_title = NULL,
  label_size = 3.5,
  filter_uncultured = FALSE,
  save_table = FALSE,
  table_filename = "aldex_pval_effect.txt"
)
```

## Arguments

- table:

  A data frame with taxa in rows and samples in columns. The last column
  must be named `taxonomy`, containing full taxonomic strings.

- metadata:

  A data frame containing sample metadata. Its first column must hold
  the sample IDs (the column names of `table`).

- group_col:

  Character. Name of the column in `metadata` that defines the groups.
  It must have exactly two groups.

- type:

  Type of plot to generate: `"volcano"` (default, volcano plot) or
  `"effect"` (effect-size plot). Case-insensitive.

- col_inf:

  Character. Color of the taxa lower in `cond`. Default `"#56B4E9"`
  (Okabe-Ito blue).

- col_sup:

  Character. Color of the taxa higher in `cond`. Default `"#E69F00"`
  (Okabe-Ito orange).

- threshold_lower:

  Numeric. Lower cutoff on the x-axis (effect size or difference between
  groups), drawn as a dashed line. Default `-1.5`.

- threshold_upper:

  Numeric. Upper cutoff on the x-axis (effect size or difference between
  groups), drawn as a dashed line. Default `1.5`.

- cond:

  Name of the reference condition for the plot labels: positive values
  on the x-axis mean higher in `cond`. Default `NULL`, the condition of
  the first sample in `table`.

- pval_threshold:

  Numeric. P-value cutoff (after `p_adjust_method`), drawn as the dashed
  horizontal line. Default `0.05`.

- p_adjust_method:

  Character. `"BH"` (default) or `"none"`: use the Benjamini-Hochberg
  adjusted p-values of ALDEx2 (`wi.eBH`) or the raw ones (`wi.ep`) for
  the y-axis and `pval_threshold`. ALDEx2 only computes the BH
  correction.

- show_labels:

  Logical. If `TRUE` (default), shows the "Higher in"/"Lower in" `cond`
  labels (only for `type = "effect"`).

- taxa:

  Optional data frame with columns `Feature.ID` and `Taxon`, used to
  label the taxa. If `NULL` (default), the `taxonomy` column of `table`
  is used.

- x_axis_title, y_axis_title:

  Titles for the x- and y-axis. Default `NULL`: `"Effect size"`
  (`"effect"`) or log2 fold change (`"volcano"`) on x, and the
  (adjusted) p-value on y.

- label_size:

  Numeric. Font size of the taxon labels drawn on significant points in
  the "volcano" plot. Default `3.5`.

- filter_uncultured:

  Logical. If `TRUE`, taxa whose name contains "uncultured"/"unculture"
  are still plotted as points but not labelled, keeping the volcano
  plot's text annotations to named taxa. Default `FALSE`.

- save_table:

  Logical. If `TRUE`, saves the ALDEx2 results table as a tab-delimited
  file. Default `FALSE`. `effect` and `diff.btw` are saved as ALDEx2
  returns them (second group in alphabetical order minus the first); in
  the plot, positive means higher in `cond`.

- table_filename:

  Character. Name or path of the saved file (used when
  `save_table = TRUE`). Default `"aldex_pval_effect.txt"`.

## Value

A ggplot object with the selected plot.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)

# group_col must have exactly two groups; Location has two
# (Rhizosphere and Roots) in the bundled example data
aldex_volcano_plot(
  table           = table,
  metadata        = metadata,
  group_col        = "Location",
  type            = "effect",
  col_inf         = "#56B4E9",
  col_sup         = "#E69F00",
  threshold_lower = -0.5,
  threshold_upper = 0.5,
  cond            = "Rhizosphere",
  show_labels     = TRUE
)
#> aldex.clr: generating Monte-Carlo instances and clr values
#> conditions vector supplied
#> operating in serial mode
#> computing center with all features
#> aldex.ttest: doing t-test
#> aldex.effect: calculating effect sizes
```
