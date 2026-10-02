# Generate a Venn diagram of taxa shared between sample groups

This function creates a Venn diagram using either the `ggVennDiagram` or
`ggenn` package based on the user's preference. It includes additional
customization options like minimum prevalence filtering and custom group
color scales (manual).

## Usage

``` r
venn_plot(
  table,
  metadata,
  merge_by = NULL,
  selected_samples = NULL,
  min_prevalence = 0,
  title = NULL,
  method = "ggvenn",
  group_colors = NULL,
  save_table = FALSE,
  table_filename = "venn_taxa_sets.txt"
)
```

## Arguments

- table:

  A data frame containing taxonomic abundance data with a column named
  `taxonomy` and subsequent columns as sample IDs.

- metadata:

  A data frame containing metadata with a column named `SAMPLEID` that
  matches the sample columns in `table`.

- merge_by:

  A character string specifying the metadata column by which to group
  and merge samples.

- selected_samples:

  Optional character vector specifying a subset of sample IDs to include
  in the analysis.

- min_prevalence:

  Optional numeric value (0-1) to filter taxa based on minimum
  prevalence across groups.

- title:

  Optional character string for the title of the plot.

- method:

  Character. Package used to draw the Venn diagram: `"ggvenn"` (default)
  or `"ggVennDiagram"`. Case-insensitive.

- group_colors:

  Optional vector of colors for the groups, either named after the
  groups or in the order of the groups. If NULL, the colorblind-friendly
  Okabe-Ito palette is used.

- save_table:

  Logical. If `TRUE`, saves a long-format table of taxa membership per
  group to disk. Default `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"venn_taxa_sets.txt"`.

## Value

A ggplot object or other plot depending on the method.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SAMPLEID"

venn_plot(
  table          = table,
  metadata       = metadata,
  merge_by       = "Location",
  min_prevalence = 0
)


## Filtering taxa by prevalence
venn_plot(
  table          = table,
  metadata       = metadata,
  merge_by       = "Location",
  min_prevalence = 0.2
)


## Custom colors
venn_plot(
  table          = table,
  metadata       = metadata,
  merge_by       = "Location",
  min_prevalence = 0,
  group_colors   = c("#1B9E77", "#D95F02")
)
```
