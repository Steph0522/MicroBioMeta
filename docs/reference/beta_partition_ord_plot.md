# Beta diversity partition (Jaccard/Sorensen)

This function computes beta diversity partition (Jaccard or Sorensen)
using the betapart package, and plots the ordination (PCoA/NMDS) with
gg_ordiplot.

## Usage

``` r
beta_partition_ord_plot(
  table,
  metadata,
  family = "jaccard",
  group_col = NULL,
  shape_col = NULL,
  legend_title = NULL,
  point_size = 3,
  group_colors = NULL,
  panel_label_case = "upper",
  panel_labels = NULL,
  panel_label_bold = TRUE,
  save_table = FALSE,
  table_filename = "SAMPLE1"
)
```

## Arguments

- table:

  Abundance matrix or data frame with taxa/features as rows and samples
  as columns (same orientation as the rest of the package). If a
  taxonomy column is present it is detected and removed automatically.

- metadata:

  Data frame with sample metadata. First column must be SampleID.

- family:

  Dissimilarity family for the partition: `"jaccard"` (default) or
  `"sorensen"`. Case-insensitive.

- group_col:

  Column in metadata to use as color grouping.

- shape_col:

  Optional column in metadata for point shapes.

- legend_title:

  Optional custom legend title.

- point_size:

  Numeric. Size of points in ordination plots. Default `3`.

- group_colors:

  Optional named vector of colors for groups.

- panel_label_case:

  Character. Case of the auto-generated A/B/C panel tags. One of
  `"upper"` (default, "A", "B", "C") or `"lower"` ("a", "b", "c").
  Ignored if `panel_labels` is supplied.

- panel_labels:

  Optional character vector of 3 custom panel tags (one per
  jaccard/turnover/nestedness panel), used as-is (e.g.
  `c("(a)", "(b)", "(c)")` or `c("a.", "b.", "c.")`) — for journal
  styles that `panel_label_case` alone can't produce. Overrides
  `panel_label_case` when provided.

- panel_label_bold:

  Logical. If `TRUE` (default), panel tags are bold. Set to `FALSE` for
  journals that require plain (non-bold) panel tags.

- save_table:

  Logical. If `TRUE`, saves the dissimilarity table to disk. Default
  `FALSE`.

- table_filename:

  Character. Base name for the saved table file. Default `"SAMPLE1"`.

## Value

A `patchwork` object with the three partition ordinations
(Jaccard/Sorensen, turnover, nestedness). It can still be modified:
`p & theme(...)` changes every panel, `p[[2]] + labs(...)` one.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

# Always tagged A/B/C (see panel_label_case and panel_labels to customize)
beta_partition_ord_plot(
  table      = table,
  metadata   = metadata,
  group_col  = "Location",
  point_size = 4
)
```
