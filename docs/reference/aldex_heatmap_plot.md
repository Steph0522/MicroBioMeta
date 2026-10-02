# ALDEx2 differential abundance heatmap

Runs ALDEx2 on a counts table and metadata and returns a ComplexHeatmap
showing differentially abundant taxa, their effect size, p-value, and
difference between groups.

## Usage

``` r
aldex_heatmap_plot(
  table,
  metadata,
  group_col,
  effect_threshold = 0.8,
  pval_threshold = NULL,
  p_adjust_method = "BH",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  heatmap_colors = NULL,
  effect_colors = circlize::colorRamp2(c(-1.5, 0, 1.5), c("#009E73", "white", "#CC79A7")),
  pvalue_colors = list(`p-value` = c(`<0.001` = "#000000", `<0.01` = "#D55E00", `<0.05` =
    "#F0E442", `>0.05` = "grey85")),
  group_colors = NULL,
  save_table = FALSE,
  table_filename = "aldex_pval_effect.txt",
  draw = TRUE,
  ...
)
```

## Arguments

- table:

  Data frame with taxa as rows and samples as columns. Must contain
  exactly one taxonomy column (named "taxonomy", "Taxonomy", "taxon",
  "taxa", "Taxa", or "Taxon").

- metadata:

  Data frame with one row per sample. Must contain the column specified
  in `group_col`.

- group_col:

  Character. Name of the column in `metadata` that defines the two
  groups to compare. Exactly two unique values are required.

- effect_threshold:

  Numeric. Minimum absolute effect size to retain (default `0.8`).

- pval_threshold:

  Numeric or NULL. Maximum p-value to retain (adjusted or not, depending
  on `p_adjust_method`). If NULL (default) only `effect_threshold` is
  applied.

- p_adjust_method:

  `"BH"` (default) or `"none"`: whether `pval_threshold` and the p-value
  annotation use ALDEx2's Benjamini-Hochberg adjusted p-values
  (`wi.eBH`) or the raw ones (`wi.ep`). ALDEx2 only computes the BH
  correction.

- cluster_rows:

  Logical. Cluster heatmap rows (default `FALSE`).

- cluster_columns:

  Logical. Cluster heatmap columns (default `FALSE`).

- heatmap_colors:

  Controls the color scale of the main heatmap body (CLR values). Three
  options: `NULL` (default, same as `"viridis"`) uses a sequential
  colorblind-friendly viridis scale (the same family used in
  `abundance_heatmap_plot`), with range computed automatically from the
  data; a preset name string — `"viridis"` or one of the diverging
  presets `"BuOr"` (blue-orange), `"BuVm"` (blue-vermillion), `"BuPk"`
  (blue-pink), `"GnPk"` (green-pink); or a
  [`circlize::colorRamp2`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
  function for full manual control.

- effect_colors:

  Color function for the effect size annotation strip (default:
  green-white-pink colorblind-friendly scale, distinct from the heatmap
  body and p-value defaults).

- pvalue_colors:

  Named list of colors for the p-value annotation strip (default:
  black/vermillion/yellow/grey categorical scale, distinct from the
  heatmap body and effect size defaults).

- group_colors:

  Optional character vector of colors for the difference barplot
  annotation, either named after the conditions (e.g.
  `c(Rhizosphere = "#56B4E9", Roots = "#009E73")`) or one color per
  condition in the order they appear in `metadata[[group_col]]`. If
  `NULL` (default) an orange/blue colorblind-friendly palette is used
  (matching the col_sup/col_inf convention used elsewhere, e.g.
  `aldex_volcano_plot`), cycling through the rest of the Okabe-Ito
  palette as needed for more than two groups.

- save_table:

  Logical. If `TRUE`, saves the underlying ALDEx2 results table
  (filtered taxa, effect size, diff.btw, p-value category) to disk.
  Default `FALSE`. `effect` and `diff.btw` are saved as ALDEx2 returns
  them (alphabetically second group minus the first); the `seccion`
  column says in which group each taxon is higher.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"aldex_pval_effect.txt"`.

- draw:

  Logical. If `TRUE` (default), the heatmap is drawn on the current
  device. Use `FALSE` to only build the returned grob without drawing
  it.

- ...:

  Old names of renamed arguments (`col_cond`, `pvalue_BH`), still
  accepted with a warning. Any other extra argument is an error.

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
colnames(metadata)[1] <- "SampleID"

# group_col must have exactly two groups; Location has two
# (Rhizosphere and Roots) in the bundled example data. effect_threshold
# alone (no pval_threshold) is used here since this small (46-sample) dataset
# rarely has taxa that pass both an effect-size and a significance
# threshold at once - combine both for a stricter, real analysis.
aldex_heatmap_plot(
  table            = table,
  metadata         = metadata,
  group_col         = "Location",
  effect_threshold = 0.5
)
#> aldex.clr: generating Monte-Carlo instances and clr values
#> conditions vector supplied
#> operating in serial mode
#> computing center with all features
#> aldex.ttest: doing t-test
#> aldex.effect: calculating effect sizes

```
