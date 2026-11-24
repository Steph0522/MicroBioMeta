# Pipe operator

See `magrittr::%>%` for details. Re-exported here so every function's
`@examples` (most of which chain steps with `%>%`) work with just
[`library(MicroBioMeta)`](https://github.com/Steph0522/MicroBioMeta),
without also requiring [`library(dplyr)`](https://dplyr.tidyverse.org).

## Usage

``` r
lhs %>% rhs
```

## Arguments

- lhs:

  A value.

- rhs:

  A function call using the magrittr semantics.

## Value

`rhs(lhs)`, i.e. `lhs` piped into `rhs`.

## Examples

``` r
c(1, 2, 3) %>% sum()
#> [1] 6
```
