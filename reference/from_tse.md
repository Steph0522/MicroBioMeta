# Convert a (Tree)SummarizedExperiment object to MicroBioMeta tables

Extracts one assay (the counts), the taxonomy (`rowData`) and the sample
metadata (`colData`) of a `TreeSummarizedExperiment` or any other
`SummarizedExperiment`, e.g. one built with the mia package, into the
plain data frames that every MicroBioMeta function takes.

## Usage

``` r
from_tse(tse, assay_name = "counts")
```

## Arguments

- tse:

  A `TreeSummarizedExperiment` or `SummarizedExperiment` object with
  taxa in rows and samples in columns.

- assay_name:

  Character. Name of the assay with the counts. Default `"counts"`.

## Value

A list with two data frames:

- `table`:

  Taxa in rows, samples in columns, and a last column `taxonomy` with
  the full taxonomic string (e.g.
  `"k__Bacteria; p__Firmicutes; c__Bacilli"`).

- `metadata`:

  Sample metadata whose first column, `SAMPLEID`, holds the sample names
  of `table`.

The taxonomy is built from the `rowData` columns named after taxonomic
ranks (Kingdom or Domain, Phylum, Class, Order, Family, Genus, Species;
any case). If `rowData` already has a `taxonomy` column, it is used as
is.

## See also

[`from_phyloseq`](https://steph0522.github.io/MicroBioMeta/reference/from_phyloseq.md)
for `phyloseq` objects.

## Examples

``` r
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    counts <- matrix(c(10, 0, 5, 3, 8, 1),
        nrow = 3,
        dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2"))
    )
    tax <- data.frame(
        Kingdom = "Bacteria",
        Phylum = c("Firmicutes", "Proteobacteria", NA),
        row.names = rownames(counts)
    )
    samples <- data.frame(
        Soil = c("Rhizosphere", "Bulk soil"),
        row.names = c("S1", "S2")
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts), rowData = tax, colData = samples
    )
    mbm <- from_tse(se)
    mbm$table
    mbm$metadata
}
#>   SAMPLEID        Soil
#> 1       S1 Rhizosphere
#> 2       S2   Bulk soil
```
