# MicroBioMeta 0.99.0

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.
* Bundled a real (subset) 16S example dataset in `inst/extdata`
  (`tabla_bacteria.txt` / `metadata_bacteria.txt`), used throughout the
  documentation examples.

SIGNIFICANT USER-VISIBLE CHANGES

* Standardized the optional table-export interface across all plotting and
  processing functions: every function now uses `save_table`/`table_filename`
  (default `FALSE`) instead of the previous, inconsistent argument names
  (e.g. `export_txt`/`file_name`).
* Renamed the example metadata column `Type_of_soil` to `Location`, and
  recoded `Treatment` from numeric codes (1/2/3) to descriptive labels
  (`Control`, `Moderate_drought`, `Severe_drought`).
* `beta_ord_plot()` now displays the NMDS stress value on the plot when
  `ordination = "NMDS"`.
* `beta_test_table()` now labels `betadisper` output with the actual grouping
  variable name instead of vegan's generic `"Groups"` label.
* Added `panel_label_case`/`panel_labels`/`panel_label_bold` arguments to
  multi-panel plotting functions for journal-style figure labeling, and
  removed hardcoded letter prefixes (a./b./c.) from lollipop-style plots.
* Ratio and lollipop plots now share a consistent visual style, and legend
  titles are derived dynamically from the chosen condition column instead of
  being hardcoded.
* Genus/species-level taxon labels are now italicized across plots.

BUG FIXES

* Fixed a crash in `abundance_heatmap_plot()` when the phylum annotation
  contained `NA` values, when more than 8 groups were requested from
  `RColorBrewer`, and when annotation columns were numeric.
* Fixed a bug in `corr_env_abund_plot()` where a `for (taxon in ...)` loop
  variable shadowed the `taxon` object used downstream, causing a
  `auto_copy(): x and y must share the same src` error.
* Fixed taxonomy columns being dropped unintentionally in `beta_turnover_plot()`.
* Fixed italic formatting of Hill-number (q0/q1/q2) facet labels in
  `beta_turnover_plot()`.
* Fixed `beta_partition_ord_plot()`'s shape legend showing the raw
  `if (!is.null(shape_col)) ...` R expression as its title instead of the
  column name, and its color/shape legends overlapping the panel titles
  (`legend.box` is now horizontal instead of stacking them vertically).
* Fixed `aldex_volcano_plot()`'s "Lower/Higher in `<condition>`" corner labels
  having backwards `hjust` values, which pushed the text past the panel edge
  and clipped it; they're now anchored to stay fully inside the panel.
* Fixed `beta_ord_plot()`'s taxon-loading labels rendering as one long
  unbroken word (e.g. `uncultured_Acidobacteriaceae`) instead of wrapping
  onto multiple lines, and added `"metagenome"` to the placeholder taxonomy
  strings excluded from vector labels.
* Fixed `aldex_heatmap_plot()` row (taxon) names never rendering despite the
  underlying data being correct - `ComplexHeatmap`'s built-in row-name
  mechanism wasn't drawing anything in this composite annotation layout, so
  row labels are now drawn via an explicit `anno_text()` row annotation
  instead; genus-level names are italicized, higher-rank fallback labels
  (e.g. `"other Chytridiomycetes"`) are not.
* Fixed `beta_turnover_plot()` printing `hillR`'s raw pairwise-comparison
  progress bar straight to the console/output.
* Standardized the package's 2-group default color pair
  (`.mbm_colors_2group`, used by `alpha_hill_plot()`, `alpha_diversity_plot()`,
  `beta_dissimilarity_plot()`, `beta_ord_plot()`, `venn_plot()`, and others)
  to use the same light blue as the qualitative 3+-group palette, instead of
  a separate darker navy - a group's color no longer changes shade
  depending on whether it's plotted alongside 1 or 3+ other groups.
* Fixed a duplicate, empty `Alpha diversity` section in `_pkgdown.yml`'s
  reference index.
* Fixed stale references to pre-rename function names
  (`alpha_hill_corrplot` -> `alpha_hill_corr_plot`, `beta_partition_plot` ->
  `beta_partition_ord_plot`) left over in code comments and one `@param` doc.

OTHER

* Removed all `install.packages()` calls from within functions; missing
  Suggested/Imported packages now fail with an informative `stop()` message
  instead of installing silently, per Bioconductor guidelines.
* Removed `assign(..., envir = .GlobalEnv)` side effects; functions return
  their results directly.
* Replaced `T`/`F` literals with `TRUE`/`FALSE`, and `1:length(x)`-style
  indexing (which breaks on zero-length input) with `seq_along()`/`seq_len()`.
* Added a `Depends: R (>= 4.1)` floor to `DESCRIPTION`, reflecting the
  package's use of the native lambda (`\(x) ...`) syntax.
* Removed `shared_plot()` and `beta_plot_flexible()`, two unexported and
  unused legacy functions superseded by `beta_dissimilarity_plot()`.
* Added `x_label_angle`, `strip_text_bold`, and `aspect_ratio` arguments to
  `alpha_hill_plot()`, `alpha_diversity_plot()`, `beta_dissimilarity_plot()`,
  `beta_turnover_plot()`, and `abundance_bar_plot()` for consistent control
  over tick-label rotation, facet strip weight, and panel proportions across
  the package's bar/boxplot functions.
* Added a `stat` argument to `beta_dissimilarity_plot()` and
  `beta_turnover_plot()` for an optional `ggpubr::stat_compare_means()`
  statistical comparison, matching the alpha diversity plotting functions.
* Added `label_size` and `filter_uncultured` arguments to
  `aldex_volcano_plot()`.
