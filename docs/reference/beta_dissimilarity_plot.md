# Beta Diversity Boxplot

This function calculates beta diversity (shared species, turnover, or
nestedness) and generates a ggplot2 boxplot with facets and custom
coloring. The user only needs to specify the metadata column(s) to be
used for comparisons; the function constructs the comparison pairs
internally and removes duplicates (A_vs_B = B_vs_A).

## Usage

``` r
beta_dissimilarity_plot(
  table,
  metadata,
  comparison_condition1 = NULL,
  condition1_col,
  condition2_col = NULL,
  facet_colors = NULL,
  group_colors = NULL,
  x_axis_title = "Condition",
  y_axis_title = NULL,
  partition = c("shared", "turnover", "nestedness"),
  family = c("sorensen", "jaccard"),
  stat = NULL,
  p_adjust_method = "holm",
  show_x_labels = TRUE,
  x_label_angle = 0,
  strip_text_bold = FALSE,
  strip_text_color = "black",
  aspect_ratio = NULL,
  save_table = FALSE,
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

  Optional column name in metadata for the second condition, used as
  facet. Each facet keeps only pairs of samples that share that value
  (e.g. two samples of treatment TC), so comparisons are split by the
  condition they come from instead of mixing them.

- facet_colors:

  Optional vector of colors for facet strips. Defaults to a neutral
  `"grey85"` background.

- group_colors:

  Optional named vector of colors for x-axis groups. Defaults to the
  package's colorblind-friendly Okabe-Ito palette (`.mbm_colors`,
  orange/blue first).

- x_axis_title:

  Title for the x-axis.

- y_axis_title:

  Title for the y-axis. Default `NULL`: built from `partition` (e.g.
  `"Beta diversity (shared features)"`).

- partition:

  Type of beta diversity to compute: `"shared"` (default), `"turnover"`,
  or `"nestedness"`. Case-insensitive.

- family:

  Dissimilarity family for the turnover/nestedness partition:
  `"sorensen"` (default) or `"jaccard"`. Case-insensitive.

- stat:

  Character or `NULL`. Statistical test to compare the boxes within each
  panel/facet. `"wilcox.test"` or `"t.test"` compare every pair of
  boxes, each with its own bracket and p-value
  ([`ggpubr::stat_pwc()`](https://rpkgs.datanovia.com/ggpubr/reference/geom_pwc.html));
  `"kruskal.test"` or `"anova"` give one global p-value per panel
  ([`ggpubr::stat_compare_means()`](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html)).
  Default `NULL` (no test shown).

- p_adjust_method:

  Multiple-comparison correction for the pairwise tests
  (`stat = "wilcox.test"` or `"t.test"`), applied within each panel; any
  method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"holm"`; `"none"` shows the raw p-values.

- show_x_labels:

  Logical. If `TRUE` (default), x-axis tick labels are shown. The
  comparison groups are also in the legend, so set to `FALSE` to hide
  the (often long) tick labels when that's redundant.

- x_label_angle:

  Numeric. Rotation (in degrees) of the x-axis tick labels when
  `show_x_labels = TRUE`. Default `0` (horizontal); use e.g. `45` or
  `90` when comparison names are long enough to overlap.

- strip_text_bold:

  Logical. If `TRUE`, facet strip labels are bold. Default `FALSE`.

- strip_text_color:

  Color of the facet strip labels, which sit on the `facet_colors`
  backgrounds. Default `"black"` (the default `facet_colors` is a light
  `"grey85"`; pass `"white"` if you supply darker `facet_colors`).

- aspect_ratio:

  Numeric. Sets the aspect ratio (height/width) of each panel. Default
  `NULL` (automatic).

- save_table:

  Logical. If `TRUE`, saves the beta diversity table to disk. Default
  `FALSE`.

- table_filename:

  Character. File path/name for the saved table. Default
  `"betadiv_table.txt"`.

## Value

A ggplot2 figure object.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

beta_dissimilarity_plot(
  table                = table,
  metadata             = metadata,
  comparison_condition1 = c("Rhizosphere_vs_Roots"),
  condition1_col       = "Location",
  condition2_col       = "Treatment",
  x_axis_title         = "Samples",
  partition            = "shared",
  family               = "sorensen"
)
```
