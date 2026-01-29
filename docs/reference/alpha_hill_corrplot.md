# Alpha diversity correlation plot

This function generates a boxplot or barplot to visualize alpha
diversity Hill numbers (q = 0, 1, 2) for a given dataset, faceted by one
or two categorical variables (e.g., sample type or treatment). It
supports palette customization, faceting, and statistical comparison.

## Usage

``` r
alpha_hill_corrplot(
  table,
  facet_orientation = "horizontal",
  plot_title = "default"
)
```

## Arguments

- table:

  A data frame or matrix with samples as columns and taxa as rows. The
  first column must contain the OTUID, ASV, or species name.

- facet_orientation:

  Whether `facet_by` appears in columns ("horizontal", default) or rows
  ("vertical").

## Value

A ggplot object showing alpha diversity with Hill numbers.

## Examples

``` r
library(vegan)
data(dune)
data(dune.env)
alpha_hill_corplot(
    table = t(dune),
    metadata = dune.env %>% tibble::rownames_to_column("SampleID"),
    facet_orientation = "horizontal"
)
#> Error in alpha_hill_corplot(table = t(dune), metadata = dune.env %>% tibble::rownames_to_column("SampleID"),     facet_orientation = "horizontal"): could not find function "alpha_hill_corplot"
```
