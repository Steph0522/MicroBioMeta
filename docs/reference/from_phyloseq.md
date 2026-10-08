# Convert a phyloseq object to MicroBioMeta tables

Extracts the count table, taxonomy and sample metadata of a `phyloseq`
object into the plain data frames that every MicroBioMeta function
takes.

## Usage

``` r
from_phyloseq(physeq)
```

## Arguments

- physeq:

  A `phyloseq` object with an OTU table and a taxonomy table and sample
  data.

## Value

A list with two data frames:

- `table`:

  Taxa in rows, samples in columns, and a last column `taxonomy` with
  the full taxonomic string (e.g.
  `"k__Bacteria; p__Firmicutes; c__Bacilli"`).

- `metadata`:

  Sample metadata whose first column, `SAMPLEID`, holds the sample names
  of `table`.

## See also

[`from_tse`](https://steph0522.github.io/MicroBioMeta/reference/from_tse.md)
for `(Tree)SummarizedExperiment` objects.

## Examples

``` r
if (requireNamespace("phyloseq", quietly = TRUE)) {
  counts <- matrix(c(10, 0, 5, 3, 8, 1), nrow = 3,
                   dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2")))
  tax <- matrix(c("Bacteria", "Bacteria", "Bacteria",
                  "Firmicutes", "Proteobacteria", NA),
                nrow = 3, dimnames = list(rownames(counts), c("Kingdom", "Phylum")))
  samples <- data.frame(Soil = c("Rhizosphere", "Bulk soil"),
                        row.names = c("S1", "S2"))
  ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE),
                           phyloseq::tax_table(tax),
                           phyloseq::sample_data(samples))
  mbm <- from_phyloseq(ps)
  mbm$table
  mbm$metadata
}
#>   SAMPLEID        Soil
#> 1       S1 Rhizosphere
#> 2       S2   Bulk soil
```
