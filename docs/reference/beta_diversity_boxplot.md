# Beta Diversity Boxplot

This function calculates beta diversity (shared species, turnover, or
nestedness) and generates a ggplot2 boxplot with facets and custom
coloring. The user only needs to specify the metadata column(s) to be
used for comparisons; the function constructs the comparison pairs
internally and removes duplicates (A_vs_B = B_vs_A).

## Usage

``` r
beta_diversity_boxplot(
  table,
  metadata,
  comparison_condition1 = NULL,
  condition1_col,
  condition2_col = NULL,
  color_facets_x = NULL,
  color_axis_x = NULL,
  title_axis_x = "Condition",
  partition = c("shared", "turnover", "nestedness"),
  family = c("sorensen", "jaccard"),
  save_table = TRUE,
  table_filename = "betadiv_table.txt"
)
```

## Arguments

- table:

  Abundance matrix (samples in columns, species/features in rows).

- metadata:

  Data frame with sample metadata. First column must match sample names
  in table.

- comparison_condition1:

  Optional vector of comparison labels for the first condition.

- condition1_col:

  Column name in metadata for the first condition.

- condition2_col:

  Optional column name in metadata for the second condition (used as
  facet).

- color_facets_x:

  Optional vector of colors for facet strips. Defaults to RColorBrewer
  "Set2".

- color_axis_x:

  Optional named vector of colors for x-axis groups. Defaults to
  RColorBrewer "Set2".

- title_axis_x:

  Title for the x-axis.

- partition:

  Type of beta diversity to compute: "shared", "turnover", or
  "nestedness".

- family:

  Family for turnover/nestedness calculation: "sorensen" or "jaccard".

## Value

A ggplot2 figure object.
