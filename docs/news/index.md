# Changelog

## MicroBioMeta 0.99.0

NEW FEATURES

- Initial submission to Bioconductor.
- Functions for metabarcoding (16S/18S/ITS) and shotgun metagenomic
  (Kraken2/Bracken) tables with a taxonomy column, plus a metadata
  table:
  - Data preparation:
    [`merge_feature_taxonomy()`](https://steph0522.github.io/MicroBioMeta/reference/merge_feature_taxonomy.md),
    [`collapse_table()`](https://steph0522.github.io/MicroBioMeta/reference/collapse_table.md),
    [`from_phyloseq()`](https://steph0522.github.io/MicroBioMeta/reference/from_phyloseq.md)
    and
    [`from_tse()`](https://steph0522.github.io/MicroBioMeta/reference/from_tse.md).
  - Composition:
    [`abundance_bar_plot()`](https://steph0522.github.io/MicroBioMeta/reference/abundance_bar_plot.md),
    [`abundance_heatmap_plot()`](https://steph0522.github.io/MicroBioMeta/reference/abundance_heatmap_plot.md)
    and
    [`abundance_sankey_plot()`](https://steph0522.github.io/MicroBioMeta/reference/abundance_sankey_plot.md).
  - Alpha diversity:
    [`alpha_hill_corr_plot()`](https://steph0522.github.io/MicroBioMeta/reference/alpha_hill_corr_plot.md),
    [`alpha_hill_plot()`](https://steph0522.github.io/MicroBioMeta/reference/alpha_hill_plot.md),
    [`alpha_diversity_plot()`](https://steph0522.github.io/MicroBioMeta/reference/alpha_diversity_plot.md),
    [`alpha_decay_plot()`](https://steph0522.github.io/MicroBioMeta/reference/alpha_decay_plot.md)
    and
    [`venn_plot()`](https://steph0522.github.io/MicroBioMeta/reference/venn_plot.md).
  - Beta diversity:
    [`beta_ord_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_ord_plot.md),
    [`beta_test_table()`](https://steph0522.github.io/MicroBioMeta/reference/beta_test_table.md),
    [`beta_partition_ord_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_partition_ord_plot.md),
    [`beta_dissimilarity_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_dissimilarity_plot.md),
    [`beta_turnover_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_turnover_plot.md)
    and
    [`beta_decay_plot()`](https://steph0522.github.io/MicroBioMeta/reference/beta_decay_plot.md).
  - Differential abundance and important taxa:
    [`aldex_volcano_plot()`](https://steph0522.github.io/MicroBioMeta/reference/aldex_volcano_plot.md),
    [`aldex_heatmap_plot()`](https://steph0522.github.io/MicroBioMeta/reference/aldex_heatmap_plot.md),
    [`ancombc_plot()`](https://steph0522.github.io/MicroBioMeta/reference/ancombc_plot.md),
    [`random_forest_lollipop_plot()`](https://steph0522.github.io/MicroBioMeta/reference/random_forest_lollipop_plot.md)
    and
    [`ratios_bubble_plot()`](https://steph0522.github.io/MicroBioMeta/reference/ratios_bubble_plot.md).
  - Environmental variables:
    [`corr_env_abund_plot()`](https://steph0522.github.io/MicroBioMeta/reference/corr_env_abund_plot.md)
    and
    [`cca_rda_biplot()`](https://steph0522.github.io/MicroBioMeta/reference/cca_rda_biplot.md).
- Bacterial phylum names are shown with their validly published names
  (e.g. Proteobacteria -\> Pseudomonadota, Firmicutes -\> Bacillota) in
  every function that uses taxonomy, with a message listing the changes.
- Example 16S, 18S and Kraken2 datasets in `inst/extdata`, and vignettes
  in English and Spanish.
