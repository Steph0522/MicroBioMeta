# Alpha diversity plot

This function generates a boxplot or barplot to visualize alpha
diversity Hill numbers (q = 0, 1, 2) for a given dataset, faceted by one
or two categorical variables (e.g., sample type or treatment). It
supports palette customization, faceting, and statistical comparison.

## Usage

``` r
alpha_diversity_plot(
  table,
  metadata,
  type = "boxplot",
  stat = NULL,
  x_col,
  fill_col,
  facet_by = NULL,
  facet_by2 = NULL,
  facet_orientation = "horizontal",
  fill_palette = "colorb",
  custom_palette = NULL,
  n_cols = NULL,
  n_rows = NULL,
  strip_color = "grey",
  show_legend = TRUE,
  figure_title = NULL,
  legend_title = NULL,
  legend_position = "bottom",
  axis_x_title = NULL,
  axis_y_title = "Diversity measure",
  free_y = FALSE,
  rarefy_depth = NULL,
  save_table = TRUE,
  table_filename = "diversity.txt"
)
```

## Arguments

- table:

  A data frame or matrix with samples as columns and taxa as rows. The
  first column must contain the OTUID, ASV, or species name.

- metadata:

  A data frame with metadata. The first column must match sample names
  in `table`.

- type:

  Type of plot: either "boxplot" or "barplot". Default is "boxplot".

- stat:

  Optional. A string indicating the test used for comparing means (e.g.,
  "wilcox.test").

- x_col:

  Column in `metadata` to be used on the x-axis.

- fill_col:

  Column in `metadata` to define fill color.

- facet_by:

  Optional. A metadata column to facet (e.g., Treatment, Site).

- facet_by2:

  Optional. A metadata column to double facet (e.g., Treatment, Site).

- facet_orientation:

  Whether `facet_by` appears in columns ("horizontal", default) or rows
  ("vertical").

- fill_palette:

  Color palette to use: "colorb", "grey", "viridis", or "brewer".
  Default: "colorb".

- custom_palette:

  A vector of custom colors. Overrides `fill_palette` if provided.

- n_cols:

  Number of columns in facet wrap (optional).

- n_rows:

  Number of rows in facet wrap (optional).

- strip_color:

  Background color of facet strips. Default: "grey".

- show_legend:

  Logical. Show legend? Default: TRUE.

- figure_title:

  Title for the entire plot.

- legend_title:

  Title for the legend.

- legend_position:

  Position of the legend: "bottom", "top", "right", or "left". Default
  is "bottom".

- axis_x_title:

  Title for the x-axis.

- axis_y_title:

  Title for the y-axis.

- free_y:

  Logical. Wether if scales in y are free or not.

## Value

A ggplot object showing alpha diversity with Hill numbers.

## Examples

``` r
library(vegan)
#> Loading required package: permute
data(dune)
data(dune.env)
alpha_hill_plot(
    table = t(dune),
    metadata = dune.env %>% tibble::rownames_to_column("SampleID"),
    x_col = "Management",
    fill_col = "Management",
    facet_by = "Use",
    facet_orientation = "vertical"
)
#> Error in dune.env %>% tibble::rownames_to_column("SampleID"): could not find function "%>%"
```
