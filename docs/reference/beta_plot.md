# Box plot of beta diversity

Box plot of beta diversity

## Usage

``` r
beta_plot(
  table,
  metadata,
  comparison_condition1,
  comparison_condition2,
  condition1.x,
  condition1.y,
  condition2.x,
  condition2.y,
  color_facets_x,
  color_axis_x
)
```

## Arguments

- table:

  table Data frame where columns are samples and rows are ASVs or taxa.

- metadata:

  A data frame with sample metadata. The first column must match sample
  names in "table".

- comparison_condition1:

  Vector of comparisons 1

- comparison_condition2:

  Vector of comparisons 2

- condition1.x:

  Condition 1 of axis-x

- condition1.y:

  Condition 1 of axis-y

- condition2.x:

  Condition 2 of axis-x

- condition2.y:

  Condition 2 of axis-y

- color_facets_x:

  Color vector for facet

- color_axis_x:

  Color and conditions vector for axis-x

- title_axis_x:

  Axis-x title

- partition:

  Component of beta diversity, e.g. "turnover", "nestedness" or "shared"

- family:

  Family for turnover o nestednes, e.g. "sorensen" or "jaccard".

## Value

A plot with the values of beta diversity.

## Examples

``` r
   res <- beta_plot_flexible(
                       table = table,
                       metadata = metadata,
                       comparison_condition1 = c("Boca_vs_L.amniotico","Boca_vs_Membrana","Boca_vs_Yema","Boca_vs_Tracto.embrionario",
                                                 "Cloaca_vs_L.amniotico","Cloaca_vs_Membrana","Cloaca_vs_Yema","Cloaca_vs_Tracto.embrionario",
                                                 "Ileon_vs_L.amniotico","Ileon_vs_Membrana","Ileon_vs_Yema","Ileon_vs_Tracto.embrionario",
                                                 "Dorso_vs_L.amniotico", "Dorso_vs_Membrana", "Dorso_vs_Yema", "Dorso_vs_Tracto.embrionario"),
                       comparison_condition2 = c("3_vs_3", "7_vs_7", "16_vs_16", "8_vs_8", "12_vs_12"),
                       condition1.x = "Seccion.x",
                       condition1.y = "Seccion.y",
                       condition2.x = "ID.identificador.x",
                       condition2.y = "ID.identificador.y",
                       color_facets_x = c("#5D478B", "#8B668B", "#CDB5CD"),
                       color_axis_x = c("L.amniotico"="#2F4F4F", "Tracto.embrionario"="#698B69", "Membrana"="#458B74", "Yema"="#B4EEB4"),
                       partition = "turnover",
                       family = "jaccard",
                       title_axis_x = "Maternal samples"
                             )
#> Error in beta_plot_flexible(table = table, metadata = metadata, comparison_condition1 = c("Boca_vs_L.amniotico",     "Boca_vs_Membrana", "Boca_vs_Yema", "Boca_vs_Tracto.embrionario",     "Cloaca_vs_L.amniotico", "Cloaca_vs_Membrana", "Cloaca_vs_Yema",     "Cloaca_vs_Tracto.embrionario", "Ileon_vs_L.amniotico", "Ileon_vs_Membrana",     "Ileon_vs_Yema", "Ileon_vs_Tracto.embrionario", "Dorso_vs_L.amniotico",     "Dorso_vs_Membrana", "Dorso_vs_Yema", "Dorso_vs_Tracto.embrionario"),     comparison_condition2 = c("3_vs_3", "7_vs_7", "16_vs_16",         "8_vs_8", "12_vs_12"), condition1.x = "Seccion.x", condition1.y = "Seccion.y",     condition2.x = "ID.identificador.x", condition2.y = "ID.identificador.y",     color_facets_x = c("#5D478B", "#8B668B", "#CDB5CD"), color_axis_x = c(L.amniotico = "#2F4F4F",         Tracto.embrionario = "#698B69", Membrana = "#458B74",         Yema = "#B4EEB4"), partition = "turnover", family = "jaccard",     title_axis_x = "Maternal samples"): could not find function "beta_plot_flexible"

                  # Show plot
                    res$plot
#> Error: object 'res' not found


#FUNCION RECAMBIO
```
