# Box plot of beta diversity

Box plot of beta diversity

## Usage

``` r
beta_turnover_plot(
  table,
  metadata,
  comparison_condition1,
  comparison_condition2 = NULL,
  condition1_col,
  condition2_col = NULL,
  facet_colors = NULL,
  group_colors = NULL,
  x_axis_title = "Section",
  y_axis_title = "Proportion of feature turnover",
  show_x_labels = FALSE,
  x_label_angle = 0,
  strip_text_bold = FALSE,
  strip_text_color = "white",
  aspect_ratio = NULL,
  stat = NULL,
  p_adjust_method = "holm",
  save_table = FALSE,
  table_filename = "betadiv_turnover.txt"
)
```

## Arguments

- table:

  table Data frame where columns are samples and rows are ASVs or taxa.

- metadata:

  A data frame with sample metadata. The first column must match sample
  names in "table".

- comparison_condition1:

  Vector of `"A_vs_B"` group-pair labels to keep (matched against the
  values of `condition1_col`). Matching ignores order - listing
  `"A_vs_B"` also matches pairs the pairwise self-join happened to
  record as `"B_vs_A"`, so each pair only needs to be listed once. Each
  matched pair becomes one x-axis/fill group, labelled exactly as
  written (e.g. `"Rhizosphere_vs_Roots"`), same as `condition1_group` in
  [`beta_dissimilarity_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_dissimilarity_plot.md).

- comparison_condition2:

  Optional. Same as `comparison_condition1`, for the values of
  `condition2_col` (e.g. `"TC_vs_TC"` to keep only within-group pairs of
  a given treatment). condition1 stays the main comparison
  (x-axis/legend); condition2 is the secondary one - supplying it both
  filters to the named pairs *and* facets the plot by them (one panel
  per pair, e.g. `"TC_vs_TC"`, `"TD_vs_TD"`), the same role
  `condition2_col` plays in
  [`beta_dissimilarity_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_dissimilarity_plot.md).
  Default `NULL`: every pair matching `comparison_condition1` is kept
  and the plot only facets by q.

- condition1_col:

  Metadata column (e.g. `"Type_of_soil"`) used to build
  `comparison_condition1`'s group pairs, same as `condition1_col` in
  [`beta_dissimilarity_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_dissimilarity_plot.md).

- condition2_col:

  Optional. Metadata column used to build `comparison_condition2`'s
  group pairs, and to facet the plot. Required when
  `comparison_condition2` is set. If given without
  `comparison_condition2`, the plot is faceted by it keeping only pairs
  of samples with the same value (e.g. `"TC_vs_TC"`, `"TD_vs_TD"`).

- facet_colors:

  Optional color vector for facet strips. Defaults to a neutral
  `"grey85"` background for each facet, like
  [`beta_dissimilarity_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_dissimilarity_plot.md)'s
  `facet_colors`.

- group_colors:

  Optional named color vector for x-axis groups. Defaults to the
  package's colorblind-friendly Okabe-Ito palette (`.mbm_colors`,
  orange/blue first), like
  [`beta_dissimilarity_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_dissimilarity_plot.md)'s
  `group_colors`.

- x_axis_title:

  Title for the x-axis. Default `"Section"`.

- y_axis_title:

  Title for the y-axis. Default `"Proportion of feature turnover"`.

- show_x_labels:

  Logical. If `TRUE`, x-axis tick labels are shown. Default `FALSE`,
  since the same groups are already named in the legend.

- x_label_angle:

  Numeric. Rotation (in degrees) of the x-axis tick labels when
  `show_x_labels = TRUE`. Default `0` (horizontal).

- strip_text_bold:

  Logical. If `TRUE`, facet strip labels are bold. Default `FALSE`
  (plain).

- strip_text_color:

  Color of the top (x) facet strip labels, which sit on the
  `facet_colors` backgrounds. Default `"white"`, unless `facet_colors`
  is left at its own light `"grey85"` default, in which case this
  defaults to `"black"` instead so it stays legible.

- aspect_ratio:

  Numeric. Aspect ratio (height/width) of each panel. Default `NULL`
  (automatic).

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

- save_table:

  Logical. If `TRUE`, saves the underlying turnover table to disk.
  Default `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"betadiv_turnover.txt"`.

## Value

A ggplot2 figure with beta diversity partitions across conditions.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "OTUID"

# comparison_condition1/2 match group pairs regardless of order, so listing
# "Rhizosphere_vs_Roots" also catches pairs the self-join recorded the
# other way around ("Roots_vs_Rhizosphere") - no need to list both.
beta_turnover_plot(
  table                 = table,
  metadata              = metadata,
  comparison_condition1 = c("Rhizosphere_vs_Roots", "Rhizosphere_vs_Rhizosphere"),
  comparison_condition2 = c("Control_vs_Control", "Moderate_drought_vs_Moderate_drought"),
  condition1_col        = "Location",
  condition2_col        = "Treatment",
  facet_colors        = c("#5D478B", "#8B668B"),
  group_colors          = c("Rhizosphere_vs_Roots"       = "#56B4E9",
                            "Rhizosphere_vs_Rhizosphere" = "#E69F00")
)
```
