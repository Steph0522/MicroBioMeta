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

  A data frame or matrix with samples as columns and taxa (features) as
  rows. Row names must contain OTUIDs, ASVs, or species names.

- taxonomy:

  A data frame or matrix with taxonomy information. Row names must match
  the identifiers in the table.

- save_table:

  Logical. If `TRUE`, saves the merged table to disk. Default `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"merged_feature_taxonomy.txt"`.

## Value

A data frame with counts and taxonomy merged. If only one column in the
taxonomy is present, it will be renamed to 'taxonomy'.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
full_table <- read.delim(table_path, row.names = 1, check.names = FALSE)

# Split the combined table into a counts-only feature table and a
# separate taxonomy data frame, as merge_feature_taxonomy expects them
feature_table <- full_table[, setdiff(colnames(full_table), "taxonomy")]
taxonomy_df <- full_table[, "taxonomy", drop = FALSE]

merged <- merge_feature_taxonomy(
  table    = feature_table,
  taxonomy = taxonomy_df
)
```
