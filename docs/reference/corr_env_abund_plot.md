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
  env_table,
  metadata,
  cond_vect = NULL,
  method = "spearman",
  hc.order = TRUE,
  geom = c("tile", "circle"),
  show_labels = TRUE,
  col_palette = NULL,
  invert_axes = TRUE,
  taxonomy_db = "silva",
  level = "genus",
  pval_threshold = NULL,
  save_table = TRUE,
  table_filename = "corr.txt"
)
```

## Arguments

- table:
- env_table:
- metadata:
- cond_vect:
- method:
- hc.order:
- geom:
- show_labels:
- col_palette:
- invert_axes:
- taxonomy_db:
- level:
- pval_threshold:

## Value

A ggplot2 object.

## Examples

``` r
colores<- c("pink","white","purple")
```
