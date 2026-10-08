# Sankey diagram of relative abundances

Creates an interactive Sankey diagram (with networkD3) of how the
relative abundance flows across taxonomic ranks.

## Usage

``` r
abundance_sankey_plot(
  table,
  output_file = NULL,
  maxn = 25,
  taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S"),
  taxonomy_db = "silva",
  width = NULL,
  height = NULL,
  save_table = FALSE,
  table_filename = "sankey_nodes_links.txt"
)
```

## Arguments

- table:

  A data frame with taxa in rows and samples in columns. The last column
  must be named `taxonomy`, containing full taxonomic strings.

- output_file:

  Character or `NULL`. Path of an HTML file to save the interactive
  diagram to. If `NULL` (default), nothing is written to disk; the
  diagram is only returned.

- maxn:

  Integer. Maximum number of taxa per rank. Default `25`.

- taxRanks:

  Character vector of the ranks to show, as one-letter codes: `"D"`
  (domain), `"K"` (kingdom), `"P"`, `"C"`, `"O"`, `"F"`, `"G"` and
  `"S"`. Default all of them.

- taxonomy_db:

  Character. Database the taxonomy strings come from: `"silva"`
  (default), `"gg2"` (Greengenes2, also `"gg"`), `"unite"` or
  `"Kraken2"` (also `"kraken"`). Case-insensitive.

- width, height:

  Numeric or `NULL`. Width and height of the diagram in pixels. If
  `NULL` (default), the widget's default size is used (it fills the
  available width).

- save_table:

  Logical. If `TRUE`, saves the Sankey nodes and links (one table, with
  a `table_type` column: `"node"` or `"link"`) as a tab-delimited file.
  Default `FALSE`.

- table_filename:

  Character. Name or path of the saved file (used when
  `save_table = TRUE`). Default `"sankey_nodes_links.txt"`.

## Value

Invisibly returns the Sankey diagram object (an htmlwidget that is shown
when printed). If `output_file` is given, it is also saved as an HTML
file.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

abundance_sankey_plot(
  table        = table,
  output_file  = file.path(tempdir(), "sankey_output.html"),
  maxn         = 10,
  taxRanks     = c("P", "C", "G", "S"),
  taxonomy_db  = "silva"
)
#> Sankey diagram saved to: C:\Users\HP\AppData\Local\Temp\RtmpyaHdP0/sankey_output.html
```
