# Relative abundance barplot

Generates a bar plot of relative abundance (%) for the most abundant
taxa groups across samples or sample groups.

## Usage

``` r
abundance_bar_plot(
  table,
  metadata = NULL,
  taxonomy_db = "silva",
  level = "genus",
  x_col = NULL,
  facet_by = NULL,
  width_equal = FALSE,
  label = "taxonomy",
  top_n = 15,
  x_axis_title = "Samples",
  y_axis_title = "Relative abundance (%)",
  x_label_angle = 0,
  strip_text_bold = FALSE,
  strip_color = "grey",
  aspect_ratio = NULL,
  add_remained = FALSE,
  save_table = FALSE,
  table_filename = "relative_abundance.txt"
)
```

## Arguments

- table:

  A data frame with taxa in rows and samples in columns. The last column
  must be named `taxonomy`, containing full taxonomic strings.

- metadata:

  A data frame containing sample metadata. Its first column must hold
  the sample IDs (the column names of `table`). Optional: if `NULL`
  (default), one bar is drawn per sample.

- taxonomy_db:

  Character. Reference taxonomy database. One of `"silva"` (default),
  `"gg2"` (Greengenes2; also accepts `"gg"` / `"greengenes2"`),
  `"unite"` (fungal ITS), or `"Kraken2"` (also accepts `"kraken"`).
  Case-insensitive.

- level:

  Character. Taxonomic level to collapse to. One of `"kingdom"`,
  `"phylum"`, `"class"`, `"order"`, `"family"`, `"genus"` (default), or
  `"species"`. Case-insensitive.

- x_col:

  Character. Column name in `metadata` to use for the x-axis (e.g.,
  environment, condition). If `NULL` (default), one bar per sample.

- facet_by:

  Optional. Character. Column name in `metadata` to facet the plot by
  (e.g., treatment). Default is `NULL`.

- width_equal:

  Logical. If `TRUE`, all bars have equal width regardless of sample
  count per group. Default `FALSE`.

- label:

  Character. Legend title for the taxa groups. Default is `"taxonomy"`.

- top_n:

  Integer. Number of most abundant taxa groups to display. Default is
  `15`.

- x_axis_title:

  Character. The title for the x-axis (default = "Samples")

- y_axis_title:

  Character. The title for the y-axis (default = "Relative abundance
  (%)").

- x_label_angle:

  Numeric. Rotation (in degrees) of the x-axis tick labels. Default `0`
  (horizontal); set to `45` or `90` when sample names are long enough to
  overlap.

- strip_text_bold:

  Logical. If `TRUE`, facet strip labels are bold. Default `FALSE`
  (plain).

- strip_color:

  Background color of facet strips (only used when `facet_by` is set).
  Default: `"grey"`.

- aspect_ratio:

  Numeric. Aspect ratio (height/width) of the panel. Default `NULL`
  (automatic).

- add_remained:

  Logical indicating whether to include an "Other" category to sum
  remaining groups; default is FALSE.

- save_table:

  Logical. If `TRUE`, saves the relative-abundance table to disk.
  Default `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"relative_abundance.txt"`.

## Value

A `ggplot2` object showing a stacked barplot of relative abundances.

## Details

- Relative abundances are calculated per sample (%).

- Taxa names are collapsed to the specified taxonomic `level` ("genus"
  or "phylum").

- Only the top `top_n` taxa are shown; others are filtered out.

- Samples are grouped and ordered according to `x_col`.

- Optional faceting by `facet_by` if provided.

- Taxonomic strings matching `"d__Bacteria;__;__;__;__;__"` are
  automatically removed.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

abundance_bar_plot(
  table        = table,
  metadata     = metadata,
  taxonomy_db  = "silva",
  level        = "genus",
  x_col        = "Location",
  label        = "Genus",
  facet_by    = "Treatment",
  top_n        = 30,
  add_remained = TRUE
)
```
