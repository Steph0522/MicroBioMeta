# Heatmap ALDEx

This function generates a heatmap to visualize alpha diversity Hill
numbers (q = 0, 1, 2) for a given dataset, faceted by one or two
categorical variables (e.g., sample type or treatment). It supports
palette customization, faceting, and statistical comparison.

## Usage

``` r
aldex_heatmap_plot(
  table,
  metadata,
  col_cond,
  effect_threshold = 0.8,
  pvalue_BH = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  heatmap_colors = circlize::colorRamp2(c(0, 0.5, 1), c("#cbdbcd", "#759c8a", "#00544d")),
  effect_colors = circlize::colorRamp2(c(-1.5, 0, 1.5), c("lightsalmon4", "white",
    "lightseagreen")),
  pvalue_colors = list(`p-value` = c(`<0.001` = "#C70039", `<0.01` = "#FF5733", `<0.05` =
    "#FFC300", `>0.05` = "#F9E79F")),
  treatment_colors = c(Higher = "#808000", Lower = "#1B5E20")
)
```

## Arguments

- table, :

  Data frame with taxonomy, where, the columns are the samples and rows
  are ASV's or taxa.

- pvalue_BH, :

  Value p-ajusted (optional)

- conditions, :

  Vector that defines the categories or classes to compare.

- effect, :

  Effect size (default: \>= 0.8)

## Value

A plot with the deferentially abundant taxonomic groups between two
categories or groups of samples.

## Examples

``` r
aldex_plot(table = table,
                     conditions = conditions,
                     effect = 0.8,
                     pvalue_BH = NULL)
#> Error in aldex_plot(table = table, conditions = conditions, effect = 0.8,     pvalue_BH = NULL): could not find function "aldex_plot"
```
