# Generate a Lollipop plot from random forest results

Generate a Lollipop plot from random forest results

## Usage

``` r
randomf_lollipop_plot(
  table,
  metadata,
  top_n = 15,
  size = 8,
  variable_to_predict,
  col_palette = NULL,
  legend_figure = NULL
)
```

## Arguments

- table:

  Data frame where columns are samples and rows are ASVs or taxa.

- metadata:

  Data frame containing sample metadata.

- top_n:

  Number of top features to plot (default = 15).

- size:

  the size of the point of the lollipop.

- variable_to_predict:

  Variable to predict from the metadata.

- col_palette:

  Custom color palette (optional).

- legend_figure:

  Main title for the figure.

## Value

A lollipop plot showing top important features from random forest
analysis.

## Examples

``` r
randomf_lollipop_plot(table = otu_table, 
                     metadata = sample_metadata,
                     variable_to_predict = "season",
                     legend_figure = "Top important features (Random Forest)")
#> Error: object 'otu_table' not found
```
