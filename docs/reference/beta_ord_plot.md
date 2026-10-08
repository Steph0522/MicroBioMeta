# Beta diversity ordination plot

Computes the beta diversity between samples with several distance
metrics (including the compositional CLR/Aitchison distance via ALDEx2)
and plots a PCA, PCoA or NMDS ordination.

## Usage

``` r
beta_ord_plot(
  table,
  metadata,
  distance = "compositional",
  mc_samples = 1,
  ordination = "PCA",
  group_col = NULL,
  palette = "colorb",
  shape_col = NULL,
  legend_title = NULL,
  taxonomy_db = "silva",
  arrows_size = 10,
  top_n = 5,
  title = "auto",
  save_table = FALSE,
  table_filename = "ordination_scores.txt"
)
```

## Arguments

- table:

  A data frame with taxa in rows and samples in columns. The last column
  must be named `taxonomy`, containing full taxonomic strings.

- metadata:

  A data frame containing sample metadata. Its first column must hold
  the sample IDs (the column names of `table`).

- distance:

  Distance method: one of `"euclidean"`, `"bray"`, `"jaccard"`,
  `"sorensen"`, `"compositional"` (default; CLR/Aitchison via ALDEx2),
  `"aitchison"`, or `"robust.aitchison"`. Case-insensitive. Note:
  `ordination = "PCA"` requires `distance = "compositional"`.

- mc_samples:

  Number of ALDEx2 Monte Carlo instances used when
  `distance = "compositional"`. With `1` (default) the clr values of one
  random instance are used: fast, but the result changes between runs
  (use [`set.seed()`](https://rdrr.io/r/base/Random.html)). With more,
  the clr values are averaged across instances, which gives an almost
  identical result in every run; `128` (ALDEx2's default) is suggested
  for final analyses, and takes longer. Ignored for other distances.

- ordination:

  Ordination method: one of `"PCA"` (default), `"PCoA"`, or `"NMDS"`.
  Case-insensitive.

- group_col:

  Character. Name of the column in `metadata` used to color the points;
  if it is numeric, a continuous color scale is used. Optional; `NULL`
  (default) for no groups.

- palette:

  Either a palette **name** or a **vector of fixed colors**; which scale
  it produces depends on whether `group_col` is discrete or continuous.

  - Named, discrete `group_col`: one of `"colorb"` (default; qualitative
    colorblind-friendly palette), `"grey"`, `"viridis"`, or `"brewer"`
    (`"Set2"`).

  - Named, continuous `group_col`: `"viridis"` (default;
    `option = "cividis"`) or `"gradient"`

  - Vector of colors, discrete `group_col`: used as-is, one color per
    level (`scale_*_manual`).

  - Vector of colors, continuous `group_col`: used as gradient stops
    (`scale_*_gradientn`).

- shape_col:

  Character. Name of the column in `metadata` used for the point shapes.
  Optional; `NULL` (default) for one shape.

- legend_title:

  Character. Title of the legend. If `NULL` (default), the name of
  `group_col` is used.

- taxonomy_db:

  Character. Database the taxonomy strings come from: `"silva"`
  (default), `"gg2"` (Greengenes2, also `"gg"`), `"unite"` or
  `"Kraken2"` (also `"kraken"`). Case-insensitive. Only used when
  `ordination = "PCA"`.

- arrows_size:

  Numeric. Size/length scaling factor for biplot arrows. Default `10`.

- top_n:

  Integer. Number of taxa that contribute most to the PCA, drawn as
  arrows. Default `5`.

- title:

  Character. Plot title. `"auto"` (default) shows
  `"Ordination - <distance>"`; `NULL` shows no title; any other text is
  used as the title.

- save_table:

  Logical. If `TRUE`, saves the ordination scores and loadings (one
  table) as a tab-delimited file. Default `FALSE`.

- table_filename:

  Character. Name or path of the saved file (used when
  `save_table = TRUE`). Default `"ordination_scores.txt"`.

## Value

A ggplot object.

## Details

With `distance = "compositional"` and `mc_samples = 1`, the clr values
come from one random Monte Carlo instance of
[`ALDEx2::aldex.clr()`](https://rdrr.io/pkg/ALDEx2/man/aldex.clr.function.html);
call [`set.seed()`](https://rdrr.io/r/base/Random.html) before the
function to make the result reproducible, or use `mc_samples = 128`.

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)

beta_ord_plot(
  table      = table,
  metadata   = metadata,
  distance   = "aitchison",
  ordination = "NMDS",
  group_col  = "Location",
  top_n      = 5
)
#> Run 0 stress 0.1706715 
#> Run 1 stress 0.1858971 
#> Run 2 stress 0.1964941 
#> Run 3 stress 0.2008229 
#> Run 4 stress 0.2108036 
#> Run 5 stress 0.2110827 
#> Run 6 stress 0.2022634 
#> Run 7 stress 0.1852785 
#> Run 8 stress 0.2071095 
#> Run 9 stress 0.1928814 
#> Run 10 stress 0.1847932 
#> Run 11 stress 0.2006323 
#> Run 12 stress 0.1817891 
#> Run 13 stress 0.1812133 
#> Run 14 stress 0.1805684 
#> Run 15 stress 0.1911592 
#> Run 16 stress 0.1987114 
#> Run 17 stress 0.2251462 
#> Run 18 stress 0.1939527 
#> Run 19 stress 0.1890742 
#> Run 20 stress 0.1883558 
#> Run 21 stress 0.1938252 
#> Run 22 stress 0.2100218 
#> Run 23 stress 0.1886366 
#> Run 24 stress 0.1961417 
#> Run 25 stress 0.1752081 
#> Run 26 stress 0.1828438 
#> Run 27 stress 0.2031569 
#> Run 28 stress 0.1873056 
#> Run 29 stress 0.196631 
#> Run 30 stress 0.1971796 
#> Run 31 stress 0.1968274 
#> Run 32 stress 0.1890111 
#> Run 33 stress 0.1953745 
#> Run 34 stress 0.1822461 
#> Run 35 stress 0.1864966 
#> Run 36 stress 0.2013476 
#> Run 37 stress 0.1775038 
#> Run 38 stress 0.1702792 
#> ... New best solution
#> ... Procrustes: rmse 0.05496463  max resid 0.238398 
#> Run 39 stress 0.1840863 
#> Run 40 stress 0.1858768 
#> Run 41 stress 0.1826098 
#> Run 42 stress 0.1861802 
#> Run 43 stress 0.2099587 
#> Run 44 stress 0.207177 
#> Run 45 stress 0.178743 
#> Run 46 stress 0.1854874 
#> Run 47 stress 0.1827463 
#> Run 48 stress 0.186119 
#> Run 49 stress 0.2163555 
#> Run 50 stress 0.2081207 
#> Run 51 stress 0.1938853 
#> Run 52 stress 0.2120015 
#> Run 53 stress 0.1733722 
#> Run 54 stress 0.1733489 
#> Run 55 stress 0.1835213 
#> Run 56 stress 0.190806 
#> Run 57 stress 0.2269347 
#> Run 58 stress 0.2013745 
#> Run 59 stress 0.1934049 
#> Run 60 stress 0.2279811 
#> Run 61 stress 0.2033342 
#> Run 62 stress 0.2225681 
#> Run 63 stress 0.1816522 
#> Run 64 stress 0.1775682 
#> Run 65 stress 0.1994683 
#> Run 66 stress 0.1723301 
#> Run 67 stress 0.1955567 
#> Run 68 stress 0.2069365 
#> Run 69 stress 0.2084881 
#> Run 70 stress 0.1955151 
#> Run 71 stress 0.183886 
#> Run 72 stress 0.1966279 
#> Run 73 stress 0.2066014 
#> Run 74 stress 0.1855555 
#> Run 75 stress 0.180835 
#> Run 76 stress 0.1881507 
#> Run 77 stress 0.190056 
#> Run 78 stress 0.1993236 
#> Run 79 stress 0.1837502 
#> Run 80 stress 0.1737991 
#> Run 81 stress 0.179853 
#> Run 82 stress 0.1821538 
#> Run 83 stress 0.1792971 
#> Run 84 stress 0.1831922 
#> Run 85 stress 0.2154851 
#> Run 86 stress 0.1836506 
#> Run 87 stress 0.1747959 
#> Run 88 stress 0.2179816 
#> Run 89 stress 0.2042863 
#> Run 90 stress 0.1827904 
#> Run 91 stress 0.1756128 
#> Run 92 stress 0.1718132 
#> Run 93 stress 0.2176935 
#> Run 94 stress 0.2017589 
#> Run 95 stress 0.1962462 
#> Run 96 stress 0.1954244 
#> Run 97 stress 0.1764035 
#> Run 98 stress 0.2014986 
#> Run 99 stress 0.1926154 
#> Run 100 stress 0.2081464 
#> *** Best solution was not repeated -- monoMDS stopping criteria:
#>     19: no. of iterations >= maxit
#>     81: stress ratio > sratmax
#> Coordinate system already present.
#> ℹ Adding new coordinate system, which will replace the existing one.
```
