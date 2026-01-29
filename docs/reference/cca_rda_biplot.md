# CCA/RDA Biplot with ggplot2

Performs Canonical Correspondence Analysis (CCA) or Redundancy Analysis
(RDA) based on a species abundance table and selected environmental
variables, returning a biplot with ggplot2 that visualizes sample scores
and environmental vectors.

## Usage

``` r
cca_rda_biplot(
  table,
  env_data,
  env_vars,
  method = "hell",
  metadata = NULL,
  group_col = NULL,
  group_colors = NULL,
  legend_title = NULL,
  scale_env = TRUE,
  pval_threshold = 0.05,
  show_all_env_vectors = FALSE,
  analysis = "CCA",
  seed = 126,
  scale_arrows = 1,
  title = NULL
)
```

## Arguments

- table:

  A data frame or matrix of species abundances (samples as rows, species
  as columns).

- env_data:

  A data frame of environmental variables (rows must match `table`).

- env_vars:

  A character vector with the names of environmental variables to
  include in the analysis.

- method:

  Transformation method passed to `decostand` (default is `"hell"` for
  Hellinger).

- metadata:

  Optional data frame with sample metadata for grouping in the plot.

- group_col:

  Optional name of the column in `metadata` used to define sample
  groups.

- group_colors:

  Optional named vector of colors to use for each group.

- legend_title:

  Optional custom title for the group legend.

- scale_env:

  Logical; whether to scale environmental variables (default is `TRUE`).

- pval_threshold:

  P-value threshold for selecting significant environmental variables
  (default is `0.05`).

- show_all_env_vectors:

  Logical; if TRUE, plot all environmental vectors regardless of
  significance.

- analysis:

  Either `"CCA"` or `"RDA"` (default is `"CCA"`).

- seed:

  Random seed for reproducibility (default is `126`).

- scale_arrows:

  Numeric value to scale environmental vectors in the plot.

- title:

  Optional plot title.

## Value

A `ggplot` object displaying the biplot with sample scores and
environmental vectors.
