# Compare taxon abundance ratios between two conditions

Computes relative abundances from a taxonomic abundance table and
compares two experimental conditions by calculating a directional
abundance ratio for each taxon. Taxa are ranked by mean abundance and
visualized as a bubble plot, where bubble size represents mean relative
abundance and color indicates the dominant condition.

## Usage

``` r
ratios_bubble_plot(
  table,
  metadata,
  group_col,
  condition_A,
  condition_B,
  taxonomy_db = "silva",
  top_n = 30,
  level = "genus",
  x_axis_title = "Taxon",
  legend_title = NULL,
  group_colors = NULL,
  save_table = FALSE,
  table_filename = "ratios_bubble_table.txt",
  ...
)
```

## Arguments

- table:

  A data frame containing a taxonomic abundance table with a taxonomy
  column and sample columns with numeric counts.

- metadata:

  A data frame containing sample metadata. The first column must
  correspond to sample IDs and include a column defining the
  experimental conditions.

- group_col:

  Character. Name of the metadata column defining the experimental
  condition.

- condition_A:

  Character. Name of the first condition to compare.

- condition_B:

  Character. Name of the second condition to compare.

- taxonomy_db:

  Character. Taxonomic database used for annotation.
  ("silva","Kraken2").

- top_n:

  Integer. Number of taxa with the highest mean abundance to display.

- level:

  Character. Taxonomic level to use for comparison (e.g. "phylum",
  "genus", "species").

- x_axis_title:

  Character. Label for the x-axis (taxon names).

- legend_title:

  Character. Title for the fill legend and the axis showing the dominant
  condition. Defaults to `group_col` when `NULL` (default).

- group_colors:

  Character vector of colors used to represent the dominant condition.

- save_table:

  Logical. If `TRUE`, saves the underlying ratio table to disk. Default
  `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"ratios_bubble_table.txt"`.

- ...:

  Old names of renamed arguments (`condition_col`), still accepted with
  a warning. Any other extra argument is an error.

## Value

A ggplot2 object showing abundance ratios between the two conditions.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

ratios_bubble_plot(
  table         = table,
  metadata      = metadata,
  group_col = "Location",
  condition_A   = "Rhizosphere",
  condition_B   = "Roots",
  taxonomy_db   = "silva",
  level         = "genus",
  top_n         = 20
)
```
