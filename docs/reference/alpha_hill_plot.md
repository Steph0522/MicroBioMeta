# Alpha diversity plot

This function generates a boxplot or barplot to visualize alpha
diversity Hill numbers (q = 0, 1, 2) for a given dataset, faceted by one
or two categorical variables (e.g., sample type or treatment). It
supports palette customization, faceting, and statistical comparison.

## Usage

``` r
alpha_hill_plot(
  table,
  metadata,
  type = "boxplot",
  stat = NULL,
  p_adjust_method = "holm",
  x_col,
  fill_col,
  facet_by = NULL,
  facet_by2 = NULL,
  facet_orientation = "horizontal",
  palette = "colorb",
  group_colors = NULL,
  n_cols = NULL,
  n_rows = NULL,
  strip_color = "grey",
  show_legend = TRUE,
  title = NULL,
  legend_title = NULL,
  legend_position = "bottom",
  x_axis_title = NULL,
  y_axis_title = "Effective number of features",
  free_y = FALSE,
  panel_label_case = "upper",
  panel_labels = NULL,
  panel_label_bold = TRUE,
  x_label_angle = 0,
  strip_text_bold = FALSE,
  aspect_ratio = NULL,
  save_table = FALSE,
  table_filename = "hill.txt"
)
```

## Arguments

- table:

  A data frame or matrix with samples as columns and taxa as rows. The
  first column must contain the OTUID, ASV, or species name.

- metadata:

  A data frame with metadata. The first column must match sample names
  in `table`.

- type:

  Type of plot: `"boxplot"` (default) or `"barplot"`. Case-insensitive.

- stat:

  Optional. Statistical test to compare the groups within each panel.
  `"wilcox.test"` or `"t.test"` compare every pair of groups, each with
  its own bracket and p-value
  ([`ggpubr::stat_pwc()`](https://rpkgs.datanovia.com/ggpubr/reference/geom_pwc.html));
  `"kruskal.test"` or `"anova"` give one global p-value per panel.

- p_adjust_method:

  Multiple-comparison correction for the pairwise tests
  (`stat = "wilcox.test"` or `"t.test"`), applied within each panel; any
  method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"holm"`; `"none"` shows the raw p-values. Panel tags (A, B, C...) are
  added regardless of whether `stat` is set; `stat` only adds the
  p-value annotations on top of them.

- x_col:

  Column in `metadata` to be used on the x-axis.

- fill_col:

  Column in `metadata` to define fill color.

- facet_by:

  Optional. A metadata column to facet (e.g., Treatment, Site).

- facet_by2:

  Optional. A metadata column to double facet (e.g., Treatment, Site).

- facet_orientation:

  Whether `facet_by` appears in columns (`"horizontal"`, default) or
  rows (`"vertical"`). Case-insensitive.

- palette:

  Color palette to use: `"colorb"` (default, colorblind-friendly
  Okabe-Ito), `"grey"`, `"viridis"`, or `"brewer"`. Case-insensitive.

- group_colors:

  A vector of custom colors. Overrides `palette` if provided.

- n_cols:

  Number of columns in facet wrap (optional).

- n_rows:

  Number of rows in facet wrap (optional).

- strip_color:

  Background color of facet strips. Default: "grey".

- show_legend:

  Logical. Show legend? Default: TRUE.

- title:

  Title for the entire plot.

- legend_title:

  Title for the legend.

- legend_position:

  Position of the legend: "bottom", "top", "right", or "left". Default
  is "bottom".

- x_axis_title:

  Title for the x-axis.

- y_axis_title:

  Title for the y-axis.

- free_y:

  Logical. Whether y-axis scales are free across facets. Default
  `FALSE`.

- panel_label_case:

  Character. Case of the auto-generated panel tags (A, B, C... added to
  every panel by default). One of `"upper"` (default, "A", "B", "C") or
  `"lower"` ("a", "b", "c"). Ignored if `panel_labels` is supplied.

- panel_labels:

  Optional character vector of custom panel tags, one per panel, used
  as-is (e.g. `c("(a)", "(b)", "(c)")` or `c("a.", "b.", "c.")`) — for
  journal styles that `panel_label_case` alone can't produce. Overrides
  `panel_label_case` when provided.

- panel_label_bold:

  Logical. If `TRUE` (default), panel tags are bold. Set to `FALSE` for
  journals that require plain (non-bold) panel tags.

- x_label_angle:

  Numeric. Rotation (in degrees) of the x-axis tick labels. Default `0`
  (horizontal); use e.g. `45` or `90` when group names are long enough
  to overlap.

- strip_text_bold:

  Logical. If `TRUE`, facet strip labels are bold. Default `FALSE`
  (plain).

- aspect_ratio:

  Numeric. Aspect ratio (height/width) of each panel. Default `NULL`,
  which lets the panels fill the available space.

- save_table:

  Logical. If `TRUE`, saves the diversity table to disk. Default
  `FALSE`.

- table_filename:

  Character. File path/name for the saved table. Default `"hill.txt"`.

## Value

A ggplot object showing alpha diversity with Hill numbers.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

# Panel tags (A/B/C) are always added (see panel_label_case and
# panel_labels to customize their case/format); stat additionally
# overlays p-value annotations on each panel
alpha_hill_plot(
  table            = table,
  metadata         = metadata,
  type             = "boxplot",
  x_col            = "Location",
  fill_col         = "Location",
  facet_by         = "Treatment",
  legend_position  = "top",
  stat             = "kruskal.test",
  panel_label_case = "upper"
)
```
