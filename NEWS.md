# MicroBioMeta 0.99.0

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.
* Bundled a real (subset) 16S example dataset in `inst/extdata`
  (`tabla_bacteria.txt` / `metadata_bacteria.txt`), used throughout the
  documentation examples.

SIGNIFICANT USER-VISIBLE CHANGES

* `aldex_volcano_plot()` and `aldex_heatmap_plot()` gain `mc_samples`, the
  number of ALDEx2 Monte Carlo instances. The default, `128` (ALDEx2's
  default), gives the same results as before; fewer instances are faster
  (e.g. `16` for a quick look) but less stable.
* `ancombc_plot()` is about 2.5 times faster: it no longer runs the ANCOMBC2
  sensitivity analysis (`pseudo_sens`), which the plot did not use. The
  plotted taxa are the same; the saved table no longer has the `passed_ss_*`
  and `diff_robust_*` columns.
* `abundance_heatmap_plot()`: `top_n` now defaults to `15`, as in
  `abundance_bar_plot()` (before it had no default).
* `beta_partition_ord_plot()`: `group_col` is required (the dispersion and the
  lines to the centroid are computed per group; before, its `NULL` default
  failed with an unclear error), and `table_filename` defaults to
  `"beta_partition"` (files `beta_partition_jacs.txt`, `_jtus.txt`, `_jnes.txt`)
  instead of `"SAMPLE1"`.
* `venn_plot()`: `merge_by` is required (its `NULL` default failed with an
  unclear error).
* `taxonomy_db` and `level` work the same way in every function that has
  them: the same accepted values and aliases (`"silva"`, `"gg2"`/`"gg"`,
  `"unite"`, `"Kraken2"`/`"kraken"`; kingdom to species, case-insensitive)
  and a clear error for anything else. `ratios_bubble_plot()` now names taxa
  like `abundance_bar_plot()` at every level and for every database (before,
  `level = "species"` only worked for Kraken2, and class, order and family
  kept the full taxonomy string). `abundance_sankey_plot()` defaults to
  `taxonomy_db = "silva"`, like the other functions.
* `alpha_hill_corr_plot()`: `title = "auto"` (default) shows the default
  title and `NULL` shows none, as in `beta_ord_plot()` and `cca_rda_biplot()`
  (before, `"default"` and `"none"`).
* `collapse_table()` no longer takes `metadata` (it was only used to order
  the sample columns): it keeps every sample column of the table, in its
  order. Call it as `collapse_table(table, level = ...)`.
* `cca_rda_biplot()` and `corr_env_abund_plot()` take the environmental
  variables directly from `metadata`: `env_data` is now optional (default
  `NULL`) and `env_vars` names the metadata columns to use, so there is no
  need to build a separate table with row names and the same sample order.
  A separate `env_data` table still works as before.
* `abundance_bar_plot()` and `abundance_heatmap_plot()` no longer need
  `metadata`: without it (or without `x_col`), they show one bar or column
  per sample. Passing something that is not a table (e.g. R's `table()`
  function by mistake) now gives a clear error.
* `beta_ord_plot()` and `beta_test_table()` gain `mc_samples` for
  `distance = "compositional"`: the number of ALDEx2 Monte Carlo instances.
  `1` (default) uses one random instance, as before, but only that one is
  drawn (about 15x faster; with the same seed the values differ from the
  previous version, which drew 128 and used the first). Values above 1
  average the clr values over the instances (e.g. `128`, ALDEx2's default),
  which gives stable results between runs.
* Faster, with identical results: `corr_env_abund_plot()` computes the
  correlation p-values for all pairs at once (same values as `cor.test()`),
  `abundance_heatmap_plot()` reads each taxonomy string only once,
  `beta_turnover_plot()` computes the pairwise Hill partition with matrices
  (same values as `hillR::hill_taxa_parti_pairwise()`; 134 s -> 1 s with 53
  samples), and the package theme is built once per session.
* SILVA composite genus names of three or more genera are now shown as
  "<last genus> group" (e.g. "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"
  becomes "Rhizobium group") in `abundance_heatmap_plot()`,
  `abundance_bar_plot()`, `abundance_sankey_plot()`, `aldex_heatmap_plot()`
  and `aldex_volcano_plot()`; two-genus names (e.g. "Escherichia-Shigella")
  are kept whole. `abundance_heatmap_plot()` also cuts other labels longer
  than `max_label_length` and gains `composite_names` to turn the rule off.
* `abundance_heatmap_plot()`: annotation tiles have the same size as the
  heatmap cells; named `colors_condition*` are matched by name; condition 2
  defaults to the colorblind-friendly "Safe" palette; new `draw` argument
  (the heatmap is drawn only once; `draw = FALSE` builds a grob that fills
  the panel it is placed in, e.g. with `cowplot::plot_grid()`). Printing
  the returned object draws the heatmap.
* `aldex_heatmap_plot()` now filters taxa like `aldex_volcano_plot()` by
  default: `pval_threshold = 0.05` and `effect_threshold = 0` (before:
  `|effect| >= 0.8` and no p-value filter), so both show the same taxa.
* Old argument names (e.g. `col_cond`, `env_table`, `index`) are no longer
  accepted; functions no longer take `...`, so a misspelled or old argument
  name gives R's usual "unused argument" error.
* `aldex_heatmap_plot()` now returns a grob (instead of a
  `ComplexHeatmap::HeatmapList`) that can be combined with other plots
  (e.g. `cowplot::plot_grid()`) and is drawn when printed; new `draw`
  argument; `effect_colors` takes 3 plain colors and
  `pvalue_colors` a named vector (the previous forms still work).
* `abundance_sankey_plot()` no longer writes an HTML file by default
  (`output_file = NULL`) and gains `width`/`height`.
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

* `beta_ord_plot()` works without `group_col` (all points in one color);
  before it failed although `group_col` defaults to `NULL`.
* Fixed the `ancombc_plot()` heatmap (3+ groups) centering its color scale
  on the middle of the LFC range instead of 0, which made small negative
  log fold changes look enriched. White is now LFC = 0, with symmetric limits.
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

* Fewer dependencies: `reshape2`, `purrr`, `tidyselect`, `RColorBrewer` and
  `viridis` were removed from Imports (their few uses are now done with
  `tidyr`, base R and `scales`, with identical results).
* Code restyled with 4-space indentation (styler).
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
