# Sankey diagram of relative abundances

Generate a Sankey diagram from a table with taxonomy

## Usage

``` r
abundance_sankey_plot(
  table,
  output_file = NULL,
  maxn = 25,
  taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S"),
  taxonomy_db = "gg",
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

  Maximum number of taxa per level to include in the diagram (default:
  25).

- taxRanks:

  Taxonomic levels to display (default:
  c("D","K","P","C","O","F","G","S")).

- taxonomy_db:

  Reference taxonomy database whose prefix style the taxonomy strings
  follow. One of `"gg"` (default; Greengenes, also accepts `"gg2"` /
  `"greengenes2"`), `"silva"`, `"unite"`, or `"kraken2"`.
  Case-insensitive.

- width, height:

  Numeric or `NULL`. Width and height of the diagram in pixels. If
  `NULL` (default), the widget's default size is used (it fills the
  available width).

- save_table:

  Logical. If `TRUE`, saves a combined table of the Sankey nodes and
  links to disk, distinguished by a `table_type` column (`"node"` or
  `"link"`). Default `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
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
#> Sankey diagram saved to: C:\Users\HP\AppData\Local\Temp\Rtmpkxiywq/sankey_output.html
```
