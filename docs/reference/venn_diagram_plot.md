# Generate a Venn diagram of taxa shared between sample groups

This function creates a Venn diagram using either the `ggVennDiagram` or
`ggenn` package based on the user's preference. It includes additional
customization options like minimum prevalence filtering and custom group
color scales (manual).

## Usage

``` r
venn_diagram_plot(
  table,
  metadata,
  merge_by = NULL,
  selected_samples = NULL,
  min_prevalence = 0,
  title = NULL,
  method = "ggvenn",
  group_colors = NULL
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
  and merge samples. Default is `Tratamiento`.

- selected_samples:

  Optional character vector specifying a subset of sample IDs to include
  in the analysis.

- min_prevalence:

  Optional numeric value (0-1) to filter taxa based on minimum
  prevalence across groups.

- title:

  Optional character string for the title of the plot.

- method:

  Character string: `ggVennDiagram` (default) or `ggvenn`, specifying
  the package to use for Venn diagram generation.

- group_colors:

  Optional vector of colors for the groups. If NULL, a default
  `distiller` scale with `Set3` palette will be used.

## Value

A ggplot object or other plot depending on the method.
