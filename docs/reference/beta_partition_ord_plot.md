# Beta diversity partition ordination

Partitions the beta diversity (Jaccard or Sorensen) into its turnover
and nestedness components with betapart, and plots a PCoA of each
component (vegan::betadisper()) with the samples linked to their group
centroid.

## Usage

``` r
beta_partition_ord_plot(
  table,
  metadata,
  family = "jaccard",
  group_col,
  shape_col = NULL,
  legend_title = NULL,
  point_size = 3,
  group_colors = NULL,
  panel_label_case = "upper",
  panel_labels = NULL,
  panel_label_bold = TRUE,
  save_table = FALSE,
  table_filename = "beta_partition"
)
```

## Arguments

- table:

  A data frame with taxa in rows and samples in columns. The last column
  must be named `taxonomy`, containing full taxonomic strings.

- metadata:

  A data frame containing sample metadata. Its first column must hold
  the sample IDs (the column names of `table`).

- family:

  Dissimilarity family for the partition: `"jaccard"` (default) or
  `"sorensen"`. Case-insensitive.

- group_col:

  Character. Name of the column in `metadata` that defines the groups.
  Required: the dispersion and the lines to the centroid are computed
  per group.

- shape_col:

  Character. Name of the column in `metadata` used for the point shapes.
  Optional; `NULL` (default) for one shape.

- legend_title:

  Character. Title of the legend. If `NULL` (default), the name of
  `group_col` is used.

- point_size:

  Numeric. Size of the points. Default `3`.

- group_colors:

  Optional character vector of colors, one per group, named after the
  groups or in their order. If `NULL` (default), the colorblind-friendly
  Okabe-Ito palette is used.

- panel_label_case:

  Character. Case of the panel tags: `"upper"` (default; A, B, C) or
  `"lower"` (a, b, c). Ignored if `panel_labels` is given.

- panel_labels:

  Optional character vector of custom panel tags, one per panel
  (Jaccard/Sorensen, turnover, nestedness), used as-is (e.g.
  `c("(a)", "(b)", "(c)")`). Overrides `panel_label_case`.

- panel_label_bold:

  Logical. If `TRUE` (default), panel tags are bold. Set to `FALSE` for
  journals that require plain (non-bold) panel tags.

- save_table:

  Logical. If `TRUE`, saves the PCoA coordinates of the three
  ordinations as tab-delimited files. Default `FALSE`.

- table_filename:

  Character. Base name or path of the saved files (used when
  `save_table = TRUE`): three files, `<name>_jacs.txt`,
  `<name>_jtus.txt` and `<name>_jnes.txt`. Default `"beta_partition"`.

## Value

A `patchwork` object with the three ordinations (Jaccard/Sorensen,
turnover and nestedness).

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)

beta_partition_ord_plot(
    table      = table,
    metadata   = metadata,
    group_col  = "Location",
    point_size = 4
)
```
