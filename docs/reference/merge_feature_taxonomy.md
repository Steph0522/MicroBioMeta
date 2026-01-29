# Merge feature table and taxonomy

This function joins a feature table (counts per sample) with the
corresponding taxonomy data.

## Usage

``` r
merge_feature_taxonomy(table, taxonomy)
```

## Arguments

- table:

  A data frame or matrix with samples as columns and taxa (features) as
  rows. Row names must contain OTUIDs, ASVs, or species names.

- taxonomy:

  A data frame or matrix with taxonomy information. Row names must match
  the identifiers in the table.

## Value

A data frame with counts and taxonomy merged. If only one column in the
taxonomy is present, it will be renamed to 'taxonomy'.

## Examples

``` r
library(vegan)
data(dune)
data(dune.env)
# Example usage (note: taxonomy not provided in this dataset)
alpha_hill_corplot(
    table = t(dune),
    metadata = dune.env %>% tibble::rownames_to_column("SampleID"),
    facet_orientation = "horizontal"
)
#> Error in alpha_hill_corplot(table = t(dune), metadata = dune.env %>% tibble::rownames_to_column("SampleID"),     facet_orientation = "horizontal"): could not find function "alpha_hill_corplot"
```
