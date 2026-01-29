# Beta diversity partition (Jaccard/Sorensen)

This function computes beta diversity partition (Jaccard or Sorensen)
using the betapart package, and plots the ordination (PCoA/NMDS) with
gg_ordiplot.

## Usage

``` r
beta_partition_plot(
  table,
  metadata,
  index = "jaccard",
  group_col = NULL,
  shape_col = NULL,
  legend_title = NULL,
  point_size = 3,
  colors = NULL,
  save_table = TRUE,
  table_filename = "SAMPLE1"
)
```

## Arguments

- table:

  Abundance matrix (samples in rows, species/features in columns).

- metadata:

  Data frame with sample metadata. First column must be SampleID.

- index:

  Family of dissimilarity: "jaccard" (default) or "sorensen".

- group_col:

  Column in metadata to use as color grouping.

- shape_col:

  Optional column in metadata for point shapes.

- legend_title:

  Optional custom legend title.

- colors:

  Optional named vector of colors for groups.

## Value

A combined cowplot panel of beta diversity partition plots.
