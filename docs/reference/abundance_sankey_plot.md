# Genera un diagrama Sankey a partir de una tabla OTU con información taxonómica

Genera un diagrama Sankey a partir de una tabla OTU con información
taxonómica

## Usage

``` r
abundance_sankey_plot(
  table,
  output_file = "sankey.html",
  maxn = 25,
  taxRanks = c("D", "K", "P", "C", "O", "F", "G", "S"),
  taxonomy_db = "gg"
)
```

## Arguments

- table:

  Data frame con una columna de taxonomía y las demás columnas
  corresponden a muestras.

- output_file:

  Nombre del archivo HTML de salida (default: "sankey.html").

- maxn:

  Número máximo de taxones por nivel a incluir en el diagrama (default:
  25).

- taxRanks:

  Niveles taxonómicos a visualizar (default:
  c("D","K","P","C","O","F","G","S")).

- taxonomy_db:

  Base de datos taxonómica, ej: "gg" o "kraken2" (default: "gg").
