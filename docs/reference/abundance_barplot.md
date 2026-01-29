# Relative abundance barplot

Generates a bar plot of relative abundance (%) for the most abundant
taxa groups across samples or sample groups.

## Usage

``` r
abundance_barplot(
  table,
  metadata,
  taxonomy_db = "silva",
  level = "genus",
  x_col,
  facet_col = NULL,
  width_equal = FALSE,
  label = "taxonomy",
  top_n_groups = 15,
  x_axis_title = "Samples",
  add_remained = FALSE,
  save_table = TRUE,
  table_filename = "relative_abundance.txt"
)
```

## Arguments

- table:

  A data frame with taxa in rows and samples in columns. The first
  column must be named `taxonomy`, containing full taxonomic strings.

- metadata:

  A data frame containing sample metadata. Must include a `SAMPLEID`
  column matching sample names in `table`.

- taxonomy_db:

  Character. Reference taxonomy database: `"silva"` (default). Affects
  how taxonomic names are simplified in the plot.

- level:

  Character. Taxonomic level to collapse: `"genus"` (default) or
  `"phylum"`.

- x_col:

  Character. Column name in `metadata` to use for the x-axis (e.g.,
  environment, condition).

- facet_col:

  Optional. Character. Column name in `metadata` to facet the plot by
  (e.g., treatment group). Default is `NULL`.

- label:

  Character. Legend title for the taxa groups. Default is `"taxonomy"`.

- top_n_groups:

  Integer. Number of most abundant taxa groups to display. Default is
  `15`.

- x_axis_title:

  Character. The tittle that should be in the x-axis (deault =
  "Samples")

- add_remained:

  Logical indicating whether to include an "Other" category to sum
  remaining groups; default is FALSE.

## Value

A `ggplot2` object showing a stacked barplot of relative abundances.

## Details

- Relative abundances are calculated per sample (%).

- Taxa names are collapsed to the specified taxonomic `level` ("genus"
  or "phylum").

- Only the top `top_n_groups` taxa are shown; others are filtered out.

- Samples are grouped and ordered according to `x_col`.

- Optional faceting by `facet_col` if provided.

- Taxonomic strings matching `"d__Bacteria;__;__;__;__;__"` are
  automatically removed.

## Examples

``` r
# relative_abundance_plot(table = your_table, metadata = your_metadata, ...)
```
