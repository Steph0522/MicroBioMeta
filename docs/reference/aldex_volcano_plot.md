# This function generates either an effect size plot or a volcano plot based on ALDEx2 results.

This function generates either an effect size plot or a volcano plot
based on ALDEx2 results.

## Usage

``` r
aldex_volcano_plot(
  table,
  metadata,
  col_cond,
  type = "volcano",
  col_inf = "blue",
  col_sup = "red",
  threshold_lower = -1.5,
  threshold_upper = 1.5,
  cond = NULL,
  cutoff.pval = 0.05,
  show_labels = TRUE,
  taxa = NULL,
  save_table = TRUE,
  table_filename = "aldex_pval_effect.txt"
)
```

## Arguments

- table:

  Data frame with count data; columns represent samples, rows represent
  features.

- metadata:

  Data frame containing metadata for the samples.

- col_cond:

  Name of the column in `metadata` that contains the experimental
  conditions.

- type:

  Type of plot to generate: "effect" for effect size plot or "volcano"
  for volcano plot.

- col_inf:

  Color for points with effect size/difference lower than threshold
  (default for "effect" plot).

- col_sup:

  Color for points with effect size/difference higher than threshold
  (default for "effect" plot).

- threshold_lower:

  Lower threshold for effect size/difference (x-axis).

- threshold_upper:

  Upper threshold for effect size/difference (x-axis).

- cond:

  Name of the condition that appears first in `table` (used in plot
  labels).

- cutoff.pval:

  p-value cutoff for significance (default = 0.05).

- show_labels:

  Logical. Whether to display "Higher/Lower in cond" labels (for
  "effect" plot only, default is TRUE).

- taxa:

  Data frame with taxonomic information (required for "volcano" plot
  only).

## Value

A `ggplot` object with the selected plot.
