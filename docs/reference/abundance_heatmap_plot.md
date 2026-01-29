# Heatmap of relative abundance

This function create a heatmap to visualize the relative abundance of
ASV's, features o bacterial groups.

## Usage

``` r
abundance_heatmap_plot(
  table,
  metadata,
  condition1 = NULL,
  condition2 = NULL,
  condition3 = NULL,
  colors_condition1 = NULL,
  colors_condition2 = NULL,
  colors_condition3 = NULL,
  name_legend_condition1 = NULL,
  name_legend_condition2 = NULL,
  name_legend_condition3 = NULL,
  top_n,
  cluster = TRUE,
  show_column_names = TRUE
)
```

## Arguments

- table:

  Data frame with taxonomy, where, the columns are the samples and rows
  are ASV's or taxa.

- metadata:

  Data frame of characteristics or important information of the samples

- condition1:

  Variable of the first horizontal annotation

- condition2:

  Variable of the second horizontal annotation

- condition3:

  Variable of the third horizontal annotation

- colors_condition1:

  Color vector for condition 1

- colors_condition2:

  Color vector for condition 2

- colors_condition3:

  Color vector for condition 3

- name_legend_condition1:

  Title assigned to legend of condition 1

- name_legend_condition2:

  Title assigned to legend of condition 2

- name_legend_condition3:

  Title assigned to legend of condition 3

- top_n:

  Number of features to plot.

- cluster:

  Logical indicating whether to cluster rows (TRUE) or order by
  abundance (FALSE)

- show_column_names:

  Logical indicating whether to show column names (TRUE) or not (FALSE)

## Value

A plot with the fifty (XX) taxonomic groups most abundant.

## Examples

``` r
abundance_heatmap_plot(table = table_taxonomy, 
                       metadata = metadata.gestacion.recto,
                       condition1 = "poblacion",
                       condition2 = "temporada",
                       condition3 = "sexo",
                       top_n = 50,
                       cluster = TRUE,
                       show_column_names = FALSE,
                       name_legend_condition1 = "Altitude",
                       name_legend_condition2 = "Season",
                       name_legend_condition3 = "Sex",
                       colors_condition1 = c("#676778","#D9D9C2"),
                       colors_condition2 = c("#0E6251", "#1B5E20"),
                       colors_condition3 = c("#5D3277","#AF6502"))
#> Error: object 'table_taxonomy' not found
```
