# barplot ancombc

This function generates a heatmap to visualize alpha diversity Hill
numbers (q = 0, 1, 2) for a given dataset, faceted by one or two
categorical variables (e.g., sample type or treatment). It supports
palette customization, faceting, and statistical comparison.

## Usage

``` r
ancombc_plot(
  table,
  metadata,
  col_cond,
  prv_cut = 0.1,
  formula = NULL,
  rand_formula = NULL
)
```

## Arguments

- table, :

  Data frame with taxonomy, where, the columns are the samples and rows
  are ASV's or taxa.

- formula, :

  Formula of the model

- conditions, :

  Vector that defines the categories or classes to compare.

- effect, :

  Effect size (default: \>= 0.8)

- pvalue_BH, :

  Value p-ajusted (optional)

## Value

A plot with the deferentially abundant taxonomic groups between two
categories or groups of samples.

## Examples

``` r
ancomb_plot(table = table,
                     conditions = conditions,
                     effect = 0.8,
                     formula = formula, 
                     pvalue_BH = NULL)
#> Error in ancomb_plot(table = table, conditions = conditions, effect = 0.8,     formula = formula, pvalue_BH = NULL): could not find function "ancomb_plot"
```
