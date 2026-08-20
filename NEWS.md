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
* Fixed taxonomy columns being dropped unintentionally in `beta_turnover_plot()`,
  `shared_plot()`, and `beta_plot_flexible()`.
* Fixed italic formatting of Hill-number (q0/q1/q2) facet labels in
  `beta_turnover_plot()`.

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
