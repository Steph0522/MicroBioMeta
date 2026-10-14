# MicroBioMeta 0.99.0

NEW FEATURES

* Initial submission to Bioconductor.
* Functions for metabarcoding (16S/18S/ITS) and shotgun metagenomic
  (Kraken2/Bracken) tables with a taxonomy column, plus a metadata table:
  - Data preparation: `merge_feature_taxonomy()`, `collapse_table()`,
    `from_phyloseq()` and `from_tse()`.
  - Composition: `abundance_bar_plot()`, `abundance_heatmap_plot()` and
    `abundance_sankey_plot()`.
  - Alpha diversity: `alpha_hill_corr_plot()`, `alpha_hill_plot()`,
    `alpha_diversity_plot()`, `alpha_decay_plot()` and `venn_plot()`.
  - Beta diversity: `beta_ord_plot()`, `beta_test_table()`,
    `beta_partition_ord_plot()`, `beta_dissimilarity_plot()`,
    `beta_turnover_plot()` and `beta_decay_plot()`.
  - Differential abundance and important taxa: `aldex_volcano_plot()`,
    `aldex_heatmap_plot()`, `ancombc_plot()`, `random_forest_lollipop_plot()`
    and `ratios_bubble_plot()`.
  - Environmental variables: `corr_env_abund_plot()` and `cca_rda_biplot()`.
* Bacterial phylum names are shown with their validly published names
  (e.g. Proteobacteria -> Pseudomonadota, Firmicutes -> Bacillota) in every
  function that uses taxonomy, with a message listing the changes.
* Example 16S, 18S and Kraken2 datasets in `inst/extdata`, and vignettes in
  English and Spanish.
