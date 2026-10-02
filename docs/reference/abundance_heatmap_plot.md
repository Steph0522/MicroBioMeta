# Heatmap of relative abundance

This function create a heatmap to visualize the relative abundance of
ASV's, features o bacterial groups.

## Usage

``` r
abundance_heatmap_plot(
  table,
  metadata,
  condition1 = NULL,
  condition2 = NULL,
  condition3 = NULL,
  colors_condition1 = NULL,
  colors_condition2 = NULL,
  colors_condition3 = NULL,
  name_legend_condition1 = NULL,
  name_legend_condition2 = NULL,
  name_legend_condition3 = NULL,
  top_n,
  exclude_unclassified = TRUE,
  cluster = TRUE,
  show_column_names = TRUE,
  cell_size = NULL,
  annotation_height = NULL,
  save_table = FALSE,
  table_filename = "abundance_heatmap_table.txt",
  feature_prefix = "",
  max_label_length = 35,
  composite_names = TRUE,
  draw = TRUE
)
```

## Arguments

- table:

  Data frame with taxonomy, where, the columns are the samples and rows
  are ASV's or taxa.

- metadata:

  Data frame of characteristics or important information of the samples

- condition1:

  Variable of the first horizontal annotation

- condition2:

  Variable of the second horizontal annotation

- condition3:

  Variable of the third horizontal annotation

- colors_condition1:

  Color vector for condition 1. If named, colors are matched to the
  condition values by name (e.g.
  `c(Roots = "#009E73", Rhizosphere = "#56B4E9")`); otherwise they are
  assigned in alphabetical order of the values. Defaults to the
  colorblind-friendly Okabe-Ito palette.

- colors_condition2:

  Color vector for condition 2 (same rules as `colors_condition1`).
  Defaults to colors from the colorblind-friendly "Safe" palette (rose,
  indigo, olive...), which contrast with the Okabe-Ito colors of
  condition 1.

- colors_condition3:

  Color vector for condition 3 (same rules as `colors_condition1`).

- name_legend_condition1:

  Title assigned to legend of condition 1

- name_legend_condition2:

  Title assigned to legend of condition 2

- name_legend_condition3:

  Title assigned to legend of condition 3

- top_n:

  Number of features to plot.

- exclude_unclassified:

  Logical. If `TRUE` (default), taxa with no recognizable classification
  at any level (labeled "Unclassified") are dropped *before* selecting
  the `top_n` most abundant features, so `top_n` always returns
  identified taxa. Set to `FALSE` to keep the previous behavior and
  allow "Unclassified" rows into the plot.

- cluster:

  Logical indicating whether to cluster rows (TRUE) or order by
  abundance (FALSE)

- show_column_names:

  Logical indicating whether to show column names (TRUE) or not (FALSE)

- cell_size:

  Numeric or `NULL`. Side, in millimeters, of each (square) heatmap
  cell. The Phylum and condition annotations use the same size, so every
  tile in the figure has the same shape. If `NULL` (default), the
  largest size that fits the current plotting device is used; when there
  are too many samples for square cells to leave room for readable row
  labels, cells keep the available width but are made just tall enough
  for the labels. Increase the figure width (e.g. `fig.width` in R
  Markdown) to get squarer cells with many samples.

- annotation_height:

  Numeric or `NULL`. Height, in millimeters, of each column annotation
  bar (condition1/condition2/condition3). If `NULL` (default), it
  matches the cell height, so annotation tiles have exactly the same
  size and shape as the heatmap cells.

- save_table:

  Logical. If `TRUE`, saves the underlying abundance table to disk.
  Default `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"abundance_heatmap_table.txt"`.

- feature_prefix:

  Character. Prefix used to label each row in the heatmap, immediately
  before the row number (e.g. `feature_prefix = "ASV"` labels rows
  `"ASV1"`, `"ASV2"`...). Default `""` (rows are labeled just `"1"`,
  `"2"`...) since the right label depends on how `table`'s features were
  generated - set it to whatever fits (`"ASV"`, `"OTU"`, `"Taxon"`,
  `"Species"`...).

- max_label_length:

  Integer or `NULL`. Taxon names longer than this many characters are
  cut with an ellipsis in the row labels, so a single long name doesn't
  squeeze the heatmap. Use `NULL` to never cut. Default `35`.

- composite_names:

  Logical. If `TRUE` (default), SILVA composite names of three or more
  genera are labeled with the last genus plus "group" (e.g.
  "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" becomes
  "Rhizobium group"); two-genus names (e.g. "Escherichia-Shigella") are
  kept whole. The full names are kept in the saved table.

- draw:

  Logical. If `TRUE` (default), the heatmap is drawn on the current
  device. Use `FALSE` to only build the returned grob without drawing
  it, e.g. to combine it with other plots; then, unless `cell_size` is
  given, cells stretch to fill the panel they are placed in instead of
  having a fixed size.

## Value

Invisibly, a `gTree` (grid grob) with the heatmap of the `top_n` most
abundant features. Printing it (e.g. typing its name) draws the heatmap;
it can also be combined with other plots (e.g.
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html),
[`patchwork::wrap_elements()`](https://patchwork.data-imaginist.com/reference/wrap_elements.html)).

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

abundance_heatmap_plot(
  table                  = table,
  metadata               = metadata,
  condition1             = "Location",
  condition2             = "Treatment",
  condition3             = "Plot",
  top_n                  = 50,
  cluster                = TRUE,
  show_column_names      = FALSE,
  name_legend_condition1 = "Location",
  name_legend_condition2 = "Treatment",
  name_legend_condition3 = "Plot"
)
#> Warning: Note: Some bacterial phylum names have been updated to match NCBI's revised taxonomy:
#> 
#> Reference: https://ncbiinsights.ncbi.nlm.nih.gov/2021/12/10/ncbi-taxonomy-prokaryote-phyla-added/
```
