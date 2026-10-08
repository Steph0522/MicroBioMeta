# Merge feature table and taxonomy

This function joins a feature table (counts per sample) with the
corresponding taxonomy data.

## Usage

``` r
merge_feature_taxonomy(
  table,
  taxonomy,
  save_table = FALSE,
  table_filename = "merged_feature_taxonomy.txt"
)
```

## Arguments

- table:

  A data frame with features in rows and samples in columns. Row names
  must hold the feature IDs (OTUs, ASVs, species...).

- taxonomy:

  A data frame or matrix with the taxonomy. Row names must match the
  feature IDs of `table`.

- save_table:

  Logical. If `TRUE`, saves the merged table as a tab-delimited file.
  Default `FALSE`.

- table_filename:

  Character. Name or path of the saved file (used when
  `save_table = TRUE`). Default `"merged_feature_taxonomy.txt"`.

## Value

A data frame with the counts of `table` and the taxonomy columns. A
single taxonomy column is renamed to `taxonomy`.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
full_table <- read.delim(table_path, row.names = 1, check.names = FALSE)

feature_table <- full_table[, setdiff(colnames(full_table), "taxonomy")]
taxonomy_df <- full_table[, "taxonomy", drop = FALSE]

merged <- merge_feature_taxonomy(
  table    = feature_table,
  taxonomy = taxonomy_df
)
```
