# Plot Correlation Between Environmental Variables and Taxonomic Groups

This function calculates the relative abundances of taxa at a specified
taxonomic level (phylum, genus, or species) from a count table, computes
correlations between these abundances and environmental variables, and
visualizes the results as either a heatmap (tile) or a bubble plot
(circle).

## Usage

``` r
corr_env_abund_plot(
  table,
  env_data,
  metadata = NULL,
  env_vars = NULL,
  method = "spearman",
  hc.order = TRUE,
  geom = c("tile", "circle"),
  show_labels = TRUE,
  col_palette = NULL,
  diverging_palette = "BuOr",
  invert_axes = TRUE,
  taxonomy_db = "silva",
  level = "genus",
  pval_threshold = NULL,
  p_adjust_method = "BH",
  x_label_angle = 45,
  save_table = FALSE,
  table_filename = "corr.txt"
)
```

## Arguments

- table:

  A data frame containing a taxonomic abundance table with a taxonomy
  column and sample columns with numeric counts.

- env_data:

  A data frame or matrix of environmental variables, with samples as row
  names.

- metadata:

  A data frame containing sample metadata. The first column must
  correspond to sample identifiers.

- env_vars:

  Character vector of environmental variables to include in the
  correlation analysis. If NULL, all available variables are used.

- method:

  Character. Correlation method passed to stats::cor and
  stats::cor.test. Supported options include ("spearman", "pearson",
  "kendall").

- hc.order:

  Logical. If TRUE, applies hierarchical clustering to reorder taxa and
  environmental variables in the plot.

- geom:

  Character. Type of visualization to generate: `"tile"` (heatmap,
  default) or `"circle"` (bubble plot). Case-insensitive.

- show_labels:

  Logical. If TRUE, displays correlation values on the plot.

- col_palette:

  Character vector defining the color palette for correlation values. If
  NULL, the palette is chosen via `diverging_palette`.

- diverging_palette:

  Character. Name of a built-in colorblind-friendly diverging palette to
  use when `col_palette` is NULL. One of `"BuOr"` (default;
  blue-white-orange, the same Okabe-Ito blue/ orange pairing used for
  the two-group color convention elsewhere in the package, e.g.
  `aldex_volcano_plot`'s col_inf/col_sup and the Rhizosphere/Roots
  colors in the bundled examples), `"BuVm"` (blue-vermillion), `"BuPk"`
  (blue-pink), `"GnPk"` (green-pink), `"PuYl"` (purple-white-yellow,
  viridis endpoints); or `"viridis"` for the plain sequential
  purple-to-yellow scale used in `aldex_heatmap_plot` (no neutral
  midpoint - not recommended for correlations, where 0 should look
  distinct from either extreme). All presets have a true white midpoint
  at 0 except `"viridis"`.

- invert_axes:

  Logical. If TRUE, swaps x and y axes in the plot.

- taxonomy_db:

  Character. Taxonomic database whose prefix style is used for
  parsing/annotation. One of `"silva"` (default), `"gg2"` (also accepts
  `"gg"` / `"greengenes2"`), `"unite"`, or `"Kraken2"` (also accepts
  `"kraken"`). Case-insensitive.

- level:

  Character. Taxonomic level to collapse taxa to. One of `"kingdom"`,
  `"phylum"`, `"class"`, `"order"`, `"family"`, `"genus"` (default), or
  `"species"`. Case-insensitive.

- pval_threshold:

  Numeric. Optional p-value threshold to retain only taxa showing
  significant correlations with at least one environmental variable. If
  `NULL`, no significance filtering is applied.

- p_adjust_method:

  Multiple-comparison correction applied to the p-values of all taxon x
  variable correlations before filtering with `pval_threshold`; any
  method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"BH"` (false discovery rate); `"none"` uses the raw p-values.

- x_label_angle:

  Numeric. Rotation (in degrees) of the x-axis labels, so long variable
  or taxon names don't overlap. Default `45`; `0` for horizontal labels.

- save_table:

  Logical. If TRUE, saves the plotted correlations (one row per taxon
  and variable, with the taxon name, raw p-value and p-value adjusted
  with `p_adjust_method`) as a tab-delimited text file.

- table_filename:

  Character. Name of the output file used when save_table = TRUE.

## Value

A ggplot2 object.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

# env_data must have rownames matching the sample names in `table`
env_data <- metadata
rownames(env_data) <- env_data$SampleID

# Uses the default colorblind-friendly diverging palette (purple-white-yellow)
corr_env_abund_plot(
  table          = table,
  env_data      = env_data,
  metadata       = metadata,
  env_vars      = c("pH", "TOC", "FW", "Root_FW", "DW", "Root_L", "Stem_L"),
  method         = "pearson",
  geom           = "tile",
  hc.order       = FALSE,
  invert_axes    = TRUE,
  show_labels    = FALSE,
  level          = "phylum",
  taxonomy_db    = "silva"
)

# pval_threshold = 0.05 would keep only taxa with a significant correlation
# after the p_adjust_method correction (BH by default).
```
