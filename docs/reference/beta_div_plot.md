# Beta diversity plot with multiple distance and ordination methods

This function computes beta diversity using several distance metrics and
ordination methods (PCA, PCoA, NMDS). It requires an abundance table
with taxonomy, metadata, and allows customization of color and shape
aesthetics. It also supports compositional transformation via ALDEx2.

## Usage

``` r
beta_div_plot(
  table,
  metadata,
  distance = "compositional",
  ordination = "PCA",
  group_col = NULL,
  group_colors = NULL,
  shape_col = NULL,
  legend_title = NULL,
  arrows_size = 10,
  n_taxa = 5
)
```

## Arguments

- table:

  A data frame with abundances. The last column must contain taxonomy
  information.

- metadata:

  A data frame with sample metadata. The first column must contain the
  sample IDs.

- distance:

  Distance method: one of "euclidean", "bray", "jaccard", "sorensen",
  "compositional" (default), "aitchison", or "robust.aitchison".

- ordination:

  Ordination method: one of "PCA" (default), "PCoA", or "NMDS".

- group_col:

  Column in `metadata` to fill color points.

- group_colors:

  Optional named vector of colors.

- shape_col:

  Optional column in `metadata` to shape points.

- legend_title:

  Optional legend title.

- n_taxa:

  Number of top contributing taxa to display as arrows in PCA.

## Value

A `ggplot2` object.
