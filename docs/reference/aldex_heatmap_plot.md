# ALDEx2 differential abundance heatmap

Runs ALDEx2 between two groups and shows the differentially abundant
taxa as a heatmap (with ComplexHeatmap) of their CLR values, annotated
with the effect size and the p-value.

## Usage

``` r
aldex_heatmap_plot(
  table,
  metadata,
  group_col,
  effect_threshold = 0,
  pval_threshold = 0.05,
  p_adjust_method = "BH",
  mc_samples = 128,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  heatmap_colors = NULL,
  effect_colors = c("#009E73", "white", "#CC79A7"),
  pvalue_colors = c(`<0.001` = "#000000", `<0.01` = "#D55E00", `<0.05` = "#F0E442",
    `>0.05` = "grey85"),
  group_colors = NULL,
  save_table = FALSE,
  table_filename = "aldex_pval_effect.txt",
  draw = TRUE
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

- effect_threshold:

  Numeric. Minimum absolute effect size to retain. Default `0` (no
  effect-size filter), so by default the heatmap shows the same taxa
  that `aldex_volcano_plot` colors as significant.

- pval_threshold:

  Numeric or `NULL`. P-value cutoff (after `p_adjust_method`) to keep
  only significant taxa. Default `0.05`. If `NULL`, only
  `effect_threshold` is applied.

- p_adjust_method:

  Character. `"BH"` (default) or `"none"`: use the Benjamini-Hochberg
  adjusted p-values of ALDEx2 (`wi.eBH`) or the raw ones (`wi.ep`) for
  `pval_threshold` and the p-value annotation. ALDEx2 only computes the
  BH correction.

- mc_samples:

  Integer. Number of ALDEx2 Monte Carlo instances. Default `128`
  (ALDEx2's default), recommended for final analyses. Fewer instances
  are faster but the p-values and effect sizes are less stable (e.g.
  `16` for a quick look); call
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the function
  to make the result reproducible.

- cluster_rows:

  Logical. Cluster heatmap rows (default `FALSE`).

- cluster_columns:

  Logical. Cluster heatmap columns (default `FALSE`).

- heatmap_colors:

  Controls the color scale of the main heatmap body (CLR values). Three
  options: `NULL` (default, same as `"viridis"`) uses a sequential
  colorblind-friendly viridis scale (the same family used in
  `abundance_heatmap_plot`), with range computed automatically from the
  data; a preset name string: `"viridis"` or one of the diverging
  presets `"BuOr"` (blue-orange), `"BuVm"` (blue-vermillion), `"BuPk"`
  (blue-pink), `"GnPk"` (green-pink); or a
  [`circlize::colorRamp2`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
  function for full manual control.

- effect_colors:

  Three colors for the effect size annotation strip: negative, zero and
  positive effect (mapped to -1.5, 0 and 1.5). Default
  `c("#009E73", "white", "#CC79A7")` (green-white-pink,
  colorblind-friendly, distinct from the heatmap body and p-value
  defaults). A
  [`circlize::colorRamp2()`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
  function is also accepted.

- pvalue_colors:

  Named vector of colors for the p-value annotation strip, with names
  `"<0.001"`, `"<0.01"`, `"<0.05"` and `">0.05"`. Default
  black/vermillion/yellow/grey (distinct from the heatmap body and
  effect size defaults). A list with one such vector named `"p-value"`
  is also accepted.

- group_colors:

  Optional character vector of colors for the group annotation, one per
  group, named after the groups or in their order. If `NULL` (default),
  the colorblind-friendly Okabe-Ito palette is used (orange/blue for two
  groups).

- save_table:

  Logical. If `TRUE`, saves the ALDEx2 results of the plotted taxa
  (effect size, `diff.btw` and p-value category) as a tab-delimited
  file. Default `FALSE`. `effect` and `diff.btw` are saved as ALDEx2
  returns them (second group in alphabetical order minus the first); the
  `higher_in` column says in which group each taxon is higher.

- table_filename:

  Character. Name or path of the saved file (used when
  `save_table = TRUE`). Default `"aldex_pval_effect.txt"`.

- draw:

  Logical. If `TRUE` (default), the heatmap is drawn on the current
  device. Use `FALSE` to only build the returned grob without drawing
  it.

## Value

Invisibly, a `gTree` (grid grob) with the heatmap. Printing it (e.g.
typing its name) draws the heatmap; it can also be combined with other
plots (e.g.
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html),
[`patchwork::wrap_elements()`](https://patchwork.data-imaginist.com/reference/wrap_elements.html)).

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)

# group_col must have exactly two groups; Location has two
# (Rhizosphere and Roots) in the bundled example data. Here taxa are
# filtered by effect size alone (pval_threshold = NULL), since this small
# (46-sample) dataset has few taxa with significant p-values.
aldex_heatmap_plot(
    table            = table,
    metadata         = metadata,
    group_col        = "Location",
    mc_samples       = 16,
    effect_threshold = 0.5,
    pval_threshold   = NULL
)
#> Warning: values are unreliable when estimated with so few MC smps


# All the colors have colorblind-friendly defaults, but each can be set by
# hand: group_colors (difference bars, named after the groups),
# effect_colors (negative, zero and positive effect size), pvalue_colors
# (one color per p-value class) and heatmap_colors (median clr values)
# \donttest{
aldex_heatmap_plot(
    table = table,
    metadata = metadata,
    group_col = "Location",
    mc_samples = 16,
    effect_threshold = 0.5,
    pval_threshold = NULL,
    group_colors = c(Rhizosphere = "#56B4E9", Roots = "#009E73"),
    effect_colors = c("#0072B2", "white", "#E69F00"),
    pvalue_colors = c(
        "<0.001" = "black", "<0.01" = "grey30",
        "<0.05" = "grey60", ">0.05" = "grey90"
    ),
    heatmap_colors = "BuOr"
)
#> Warning: values are unreliable when estimated with so few MC smps

# }
```
