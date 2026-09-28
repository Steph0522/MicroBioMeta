# Revisión estática de rendimiento — MicroBioMeta 0.99.0

Alcance: los 25 archivos de `R/*.R` (≈8 200 líneas), `DESCRIPTION` y `NAMESPACE`, leídos completos.
Revisión **estática**: no se modificó ningún archivo del paquete. Las únicas mediciones son los tiempos de carga
de namespaces (R 4.6.1, proceso nuevo por paquete) de la sección 4.

Contexto del benchmark (`optimization/profiling/01_original_timings.R`, 6 045 features × 79 muestras):
`beta_test_table` ≈ 15 s (`method = "compositional"`, `Treatment*Type_of_soil`, 999 perm.),
`alpha_diversity_plot` ≈ 5.4 s, `abundance_heatmap_plot` ≈ 4.75 s.

Leyenda: **Impacto** alto / medio / bajo (sobre el tiempo de ejecución con datos reales).
**Esfuerzo**: S (< 30 min), M (unas horas), L (≥ 1 día / cambio de API).

---

## 0. Quick wins (top 10)

| # | Impacto | Esfuerzo | Archivo:línea | Qué hacer |
|---|---|---|---|---|
| 1 | **Alto** | S | `beta_test_table.R:148-151`, `beta_ord_plot.R:155-158` | `aldex.clr(mc.samples = 128)` genera 128 réplicas Monte Carlo y **solo se usa la 1ª** (`getMonteCarloSample(.., 1)`). Se desperdicia ~99 % del trabajo. Hay que calcular el CLR directamente (con pseudoconteo), o bajar a `mc.samples = 16` y usar la mediana de las réplicas. Además, la primera llamada paga unos 5.4 s solo por cargar el namespace de ALDEx2 (sección 4). |
| 2 | **Alto** | S | `abundance_heatmap_plot.R:384-396` | El heatmap **se dibuja dos veces**: una con `grid.grabExpr(draw())` y otra con `grid.newpage()` + `grid.draw()`. Además hace clustering y maquetación dentro de la función. Conviene devolver el objeto `Heatmap`/`HeatmapList` (o dibujarlo una sola vez). |
| 3 | **Alto** | M | `alpha_diversity_plot.R:300-481`, `alpha_hill_plot.R:284-491` | La ruta por defecto (`use_grid_compose`) construye 1 ggplot por celda (3 × niveles de `facet_by`; 12 en el benchmark). Luego **otra más** para extraer la leyenda (`build_cell(1,1)`, que repite `stat_compare_means`). `cowplot::get_legend` y `plot_grid` hacen `ggplotGrob` de todas **dentro** de la función. Mejor un solo ggplot facetado con las etiquetas A/B/C como `geom_text` (sección 1.2). Si no, al menos reutilizar `plots[[1]]` para la leyenda. |
| 4 | **Alto** | M | `beta_diversity.R:146-163` (`beta_turnover_plot`) | `hillR::hill_taxa_parti_pairwise` se llama 3 veces (q = 0, 1, 2). Hace un bucle en R sobre n(n−1)/2 pares (3 081 pares × 3 con 79 muestras). Se puede vectorizar con `tcrossprod` (sección 1.10). |
| 5 | **Alto** (con `pval_threshold`) | S | `corr_env_abund_plot.R:451-458, 463` | Doble bucle `for` con `cor.test()` por cada par (variable × taxón), es decir miles de llamadas, y la matriz `cor()` se recalcula completa. El p-valor se puede sacar vectorizado de `r` (`pt()`) y basta con subsetear la matriz ya calculada. |
| 6 | **Alto** | S | `ancombc_plot.R:177, 170-184` | `pseudo_sens = TRUE` está fijo y dispara el análisis de sensibilidad, que multiplica el tiempo de `ancombc2`. `n_cl` (paralelo) no está expuesto. Conviene exponer `pseudo_sens = FALSE` por defecto y `n_cl`. |
| 7 | **Medio-alto** | S | `random_forest_lollipop_plot.R:78-83` | `importance = TRUE` calcula la importancia por permutación (≈ 2× el coste), pero solo se usa `MeanDecreaseGini`. Poner `importance = FALSE`, pre-filtrar features raras, exponer `ntree` y fijar la semilla. |
| 8 | **Medio** | S | `abundance_heatmap_plot.R:110-129, 171-184` | `separate()` y 5 ramas de regex/`str_extract` se aplican a **todas** las features antes de `slice(top_n)`. Los 13 `case_when` de binning se hacen columna a columna sobre un data.frame. Mejor ordenar por abundancia primero, clasificar solo los candidatos y hacer el binning con `findInterval()` sobre la matriz. |
| 9 | **Medio** (instalación y primera llamada) | M | `DESCRIPTION` (Imports), `NAMESPACE` | Hay 35 Imports. ALDEx2 (5.4 s de carga), phyloseq (8.7 s), ANCOMBC (2.6 s + un árbol de instalación enorme), ComplexHeatmap, betapart, randomForest, ggvenn, ggVennDiagram, geosphere… Se pueden pasar a `Suggests` con `requireNamespace()`. |
| 10 | **Medio** | S | `collapse_table.R:83, 92`; `corr_env_abund_plot.R:185, 192`; `abundance_bar_plot.R:114-126`; `corr_env_abund_plot.R:130-142` | `vapply(strsplit())` recorre fila a fila (2 pasadas sobre la taxonomía), y hay 12 `dplyr::filter()` encadenados, cada uno copiando la tabla. Sustituir por regex vectorizada y un único `!taxonomy %in% c(...)`. |

---

## 1. Hallazgos por función (priorizados)

### 1.1 `beta_test_table` (`R/beta_test_table.R`), 15 s en el benchmark

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | S | 148-151 | `aldex.clr(t(table), mc.samples = 128)` genera 6 045 × 79 × 128 ≈ 61 M draws de Dirichlet más el CLR de cada réplica, y luego usa solo `getMonteCarloSample(aldex_obj, 1)`. Es casi seguro el principal consumidor de los 15 s. | CLR directo (ms), ver el snippet A. Si se quiere conservar la semántica de ALDEx2, exponer `mc.samples` (p. ej. 16) y promediar las réplicas: es más estable estadísticamente que usar 1 sola. |
| **Alto** (1ª llamada) | S | 148 | `ALDEx2::` no está cargado en el setup del benchmark. La primera llamada paga ≈ 5.4 s solo de `loadNamespace("ALDEx2")`. | Con el snippet A, ALDEx2 deja de ser necesario aquí. Si se mantiene, documentarlo y medir `first` aparte de `median`. |
| Medio | S | 160-164 | `adonis2(by = "terms")` con 999 permutaciones y un modelo de 3 términos (`A*B`). Es correcto, pero no expone `parallel`. | Añadir `parallel = getOption("mc.cores", 1)` (vegan lo acepta) y sugerir `permutations = 199` en exploración. |
| Medio | S | 145-147 | `set.seed` solo en la rama compositional. Las permutaciones de `adonis2`/`permutest` **no** están sembradas, así que los p-valores no son reproducibles. | Mover `restore_seed`/`set.seed(seed)` al principio de la función. |
| Bajo | S | 107-109 | `if (!all(num_cols)) { }`: bloque vacío. | Eliminar, o emitir un `warning`. |
| Bajo | S | 243-247 | Bucle `table_cell_font` por fila significativa (pocas filas). | Aceptable. |
| Bajo | S | 209-232 | Construye un `ggtexttable` (grob) aunque el usuario quiera la tabla. | Devolver también `tabla` (p. ej. como atributo `attr(tab, "data")`) o añadir un argumento `return = c("plot", "table")`. |

**Snippet A: CLR sin ALDEx2 (equivalente al `aitchison` de vegan con pseudoconteo)**
```r
# table: muestras x features (matriz numérica)
.mbm_clr <- function(x, pseudocount = 0.5) {
  l <- log(x + pseudocount)
  l - rowMeans(l)                 # centra por muestra
}
dist_matrix <- stats::dist(.mbm_clr(table))   # == vegan::vegdist(table, "aitchison", pseudocount = 0.5)
```
Opción ALDEx2 "barata y estable":
```r
aldex_obj <- ALDEx2::aldex.clr(t(table), mc.samples = mc.samples, denom = "all",
                               verbose = FALSE, useMC = FALSE)
inst  <- ALDEx2::getMonteCarloInstances(aldex_obj)       # lista por muestra: features x mc
clr_m <- t(vapply(inst, function(m) apply(m, 1, stats::median), numeric(nrow(inst[[1]]))))
dist_matrix <- stats::dist(clr_m)
```

### 1.2 `alpha_diversity_plot` (`R/alpha_diversity_plot.R`), 5.4 s en el benchmark

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | M | 203-204, 300-481 | Por defecto (`facet_by2 = NULL`, horizontal) entra en `use_grid_compose`: **12 ggplots** (4 niveles × 3 índices) combinados con `cowplot::plot_grid`, que llama `as_gtable()` sobre cada uno, es decir 12 `ggplot_build` + `ggplot_gtable` dentro de la función. | Un único ggplot con `facet_grid2(facet_by ~ Index)` (la rama "no grid" de las líneas 485-508 ya existe) y las etiquetas como capa. Ver el snippet B. El objeto vuelve a ser un ggplot "perezoso": se renderiza una vez, al imprimirlo. |
| **Alto** | S | 444 | `cowplot::get_legend(build_cell(1, 1) + ...)` **reconstruye** la celda 1 (incluido `stat_compare_means`) y la renderiza para extraer la leyenda. | `leg <- cowplot::get_legend(plots[[1]])` antes de quitar las leyendas (así lo hace `alpha_decay_plot.R:370`). |
| Medio | S | 554-608 | Ruta sin grid: `ggplotGrob(p)` + manipulación del gtable + `ggdraw()`, es decir, **render dentro de la función**, y se devuelve un objeto que ya no es ggplot (no se puede añadir `+ theme()`). | Mismo arreglo (snippet B). |
| Medio | S | 138 | `data.frame(t(table))`: transpone y convierte a data.frame (con `make.names` sobre 6 045 nombres de columna). | `x <- t(as.matrix(table))`: vegan acepta matrices. |
| Medio | S | 145-148 | `estimateR` + 2 × `diversity` recorren la matriz 3 veces (menor). | Aceptable. Si se quiere, calcular `p <- x / rowSums(x)` una vez y derivar Shannon y Simpson de ahí. |
| Bajo | S | 141 | `vegan::rrarefy` sin semilla: no reproducible. | Argumento `seed` y el patrón `.mbm_save_seed()`. |
| Bajo | S | 357-376 | Un `stat_compare_means` por celda: cada uno ejecuta su test al construir el plot, lo cual es correcto, pero con la leyenda se duplica (ver la fila de la línea 444). | Se resuelve con lo anterior. |

**Snippet B: una sola figura facetada con etiquetas A/B/C baratas**
```r
tags <- unique(results_largo[, c(facet_by, "Index")])
tags <- tags[order(tags[[facet_by]], tags$Index), ]
tags$.tag <- if (is.null(panel_labels)) LETTERS[seq_len(nrow(tags))] else panel_labels
p <- ggplot2::ggplot(results_largo, ggplot2::aes(.data[[x_col]], value, fill = .data[[fill_col]])) +
  capa_geom + capa_error + fill_scale + facet_config +
  ggplot2::geom_text(data = tags, ggplot2::aes(x = -Inf, y = Inf, label = .tag),
                     inherit.aes = FALSE, hjust = -0.4, vjust = 1.4,
                     family = "serif", fontface = "bold", size = 5)
# sin ggplotGrob / cowplot: se devuelve un ggplot normal
```
(Si hace falta que la etiqueta quede *fuera* del panel, se puede usar `ggh4x`/`patchwork::plot_annotation(tag_levels = "A")`, que también es perezoso.)

### 1.3 `alpha_hill_plot` (`R/alpha_hill_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | M | 183-184, 284-491 | El mismo patrón cowplot que 1.2, y además también para `vertical` sin `facet_by`. | Snippet B. Unas 400 líneas duplicadas con `alpha_diversity_plot` se pueden extraer a un helper común `.mbm_alpha_panels()`. |
| **Alto** | S | 453 | `get_legend(build_cell(1, 1))`: reconstrucción extra. | `get_legend(plots[[1]])`. |
| Medio | S | 126-131 | `data.frame(t(table))` y 3 llamadas a `hillR::hill_taxa`, cada una renormaliza la tabla. | Helper vectorizado (snippet C). |
| Medio | S | 572-630 | `ggplotGrob` + gtable + `ggdraw` dentro de la función. | Snippet B. |

**Snippet C: números de Hill vectorizados (reemplaza 3 × `hill_taxa`, en `alpha_hill_plot`, `alpha_hill_corr_plot` y `alpha_decay_plot`)**
```r
.mbm_hill <- function(x) {                    # x: matriz muestras x features
  p <- x / rowSums(x)
  plogp <- p * log(p); plogp[p == 0] <- 0
  data.frame(q0 = rowSums(x > 0),
             q1 = exp(-rowSums(plogp)),
             q2 = 1 / rowSums(p^2))
}
```

### 1.4 `abundance_heatmap_plot` (`R/abundance_heatmap_plot.R`), 4.75 s en el benchmark

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | S | 384-396 | `grid.grabExpr(ComplexHeatmap::draw(...))` renderiza en un dispositivo oculto y luego `grid.newpage(); grid.draw()` **renderiza otra vez**: 2 renders completos, con `rect_gp` por celda, por llamada. Además deja un efecto secundario (un `Rplots.pdf` en la raíz del repo en sesiones no interactivas). | Devolver `heats` (o la `HeatmapList`) con los parámetros de `draw` guardados. Si hace falta el grob, hacer **solo** `grid.grabExpr` y no dibujar: `out <- grid::grid.grabExpr(draw(...)); return(out)` y que el usuario haga `grid::grid.draw(out)`. O al revés: `draw()` una vez y `invisible(ht)`. |
| Medio | S | 110-129 | `mutate(abun = rowMeans(.))`, `left_join` y `separate()` en 7 columnas, más `case_when` con 10 `grepl` y 4 `str_extract` (todas las ramas se evalúan en **todas** las filas), **sobre las 6 045 features** antes de `slice(top_n)`. | Snippet D: ordenar por `rowMeans` sobre la matriz y clasificar solo las primeras `k` filas hasta reunir `top_n` clasificadas. |
| Medio | S | 132 | `mutate(across(everything(), trimws))` convierte **todas las columnas de muestra** a carácter, que luego se vuelven a convertir con `as.numeric` (líneas 170-171). | `trimws` solo en las columnas de texto. |
| Medio | S | 156-170 | Ida y vuelta data.frame → `t()` (matriz de **carácter**, porque hay columnas de texto) → data.frame → `inner_join` → `t()` → `as.numeric`. | Trabajar con la matriz numérica `ra[top_ids, samples]` y alinear la metadata con `match()`. |
| Medio | S | 171-184 | 13 ramas `case_when` por columna para discretizar. | `m[] <- findInterval(m, c(.001,.005,.01,.1,.2,1,2,5,10,25,50,75), left.open = TRUE)`, que da exactamente 0..12. |
| Bajo | S | 143-144 | `warning()` incondicional en **cada** llamada, aunque no se renombre ningún filo. | Emitirlo solo si hubo reemplazos (como en `random_forest_lollipop_plot.R:147`). |
| Bajo | S | 106 | `t(t(table_counts)/colSums(...))` sobre un data.frame (coerción implícita). | `ra <- sweep(m, 2, colSums(m), "/") * 100` con `m <- as.matrix(...)`. |

**Snippet D: filtrar antes de parsear**
```r
m   <- as.matrix(table[, setdiff(names(table), "taxonomy")])
ra  <- sweep(m, 2, colSums(m), "/") * 100
ord <- order(rowMeans(ra), decreasing = TRUE)
lab <- .mbm_tax_label(table$taxonomy[ord], level = "genus")   # helper vectorizado, sección 3
keep <- if (exclude_unclassified) ord[lab != "Unclassified"] else ord
keep <- head(keep, top_n)
heat <- ra[keep, metadata[[1]][metadata[[1]] %in% colnames(ra)], drop = FALSE]
heat[] <- findInterval(heat, c(.001,.005,.01,.1,.2,1,2,5,10,25,50,75), left.open = TRUE)
```
(Con el helper cacheado, basta con clasificar solo `head(ord, 5 * top_n)`.)

### 1.5 `abundance_bar_plot` (`R/abundance_bar_plot.R`), ya es rápida (0.25 s)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| Medio | S | 114-126 | 12 `dplyr::filter()` encadenados: 12 copias de la tabla completa. | `table <- table[!table$taxonomy %in% .mbm_uninformative_tax, ]` (constante compartida con `corr_env_abund_plot`). |
| Medio | S | 145-192 | Dos bloques casi idénticos (silva/Kraken2/gg2 y unite) de `sub()` por nivel. | Un único `switch(level, ...)` con el patrón de corte, o `.mbm_truncate_tax(tax, level_idx)` (sección 3). |
| Medio | S | 194-199, 214-234 | `group_by/summarise` sobre data.frame ancho y luego `pivot_longer` de **todas** las taxas colapsadas × muestras, antes de elegir el top-n. | Colapsar con `rowsum(m, tax)`, medias por grupo con una multiplicación de matrices y hacer `pivot_longer` solo del top-n. |
| Bajo | S | 247-250 | Con `add_remained`, `top_avg` **recalcula** desde `table_long` exactamente lo mismo que `avg_by_group` filtrado. | `top_avg <- dplyr::filter(avg_by_group, taxonomy %in% top_groups)`. |
| Bajo | M | 279-408 | 8 bloques `case_when` duplicados (se repiten en 5 archivos). Solo afectan a filas ya agregadas, así que el coste es bajo, pero dificultan el mantenimiento. | Helper `.mbm_tax_label()` (sección 3). |

### 1.6 `collapse_table` (`R/collapse_table.R`) y `merge_feature_taxonomy` (`R/merge_feature_taxonomy.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| Medio | S | `collapse_table.R:79-83` | `get_depth` con `vapply(strsplit)`: un closure de R por fila. | `depth <- lengths(regmatches(tax, gregexpr("[^;]*__[^;]*", tax)))` o `lengths(gregexpr("__", tax, fixed = TRUE))`, corrigiendo el caso `-1`. |
| Medio | S | `collapse_table.R:92-95` | Segunda pasada `vapply(strsplit + paste)` por fila para truncar. | `sub(sprintf("^((?:[^;]*;){%d}[^;]*).*$", level_idx - 1), "\\1", tax, perl = TRUE)`. |
| Medio | S | `collapse_table.R:98-105` | `group_by/summarise(across(where(is.numeric)))` sobre data.frame. | `rowsum(as.matrix(highres[, ordered_samples]), highres$taxonomy, reorder = FALSE)`, que va en C y es de 10 a 100 veces más rápido. |
| Bajo | S | `collapse_table.R:125-130` | Siempre calcula `long_format` (un `pivot_longer` completo), aunque no se use. | Argumento `long = FALSE`, o calcularlo de forma perezosa. |
| Bajo | S | `merge_feature_taxonomy.R:56-64` | `left_join` por una columna `OTUID` añadida, cuando ambas tablas ya están alineadas por `common_species`. | `cbind(table[common, , drop = FALSE], taxonomy[common, , drop = FALSE])`. |

### 1.7 `aldex_heatmap_plot` / `aldex_volcano_plot`

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | S | `aldex_heatmap_plot.R:138-147`, `aldex_volcano_plot.R:120` | `ALDEx2::aldex()` (CLR con 128 MC, t/Wilcoxon por feature y `aldex.effect`) sobre **todas** las features, sin pre-filtro de prevalencia. Features con 1-2 lecturas dominan el coste sin aportar nada. | Pre-filtro opcional (`min_prevalence = 0.1`, `min_count`) antes de `aldex()`. Exponer `mc.samples` y `useMC` (`aldex.clr`/`aldex.effect` aceptan `useMC`). |
| **Alto** (flujo) | M | ambas | Heatmap y volcano del mismo contraste **recalculan ALDEx2 dos veces**. | Aceptar un `aldex_res` precalculado (`if (inherits(table, "data.frame") && all(c("effect","wi.eBH") %in% names(table)))`) o una función pública `mbm_aldex()` que devuelva los resultados para ambos plots. |
| Medio | S | `aldex_heatmap_plot.R:375-387` | `ComplexHeatmap::draw()` dentro de la función: render como efecto secundario, y se devuelve `ht_list` (que al imprimirse se dibuja otra vez). | No dibujar; devolver `ht_list` con un método de impresión, o documentar que el render ocurre dentro. |
| Medio | S | `aldex_volcano_plot.R:120` | Sin semilla: resultados no reproducibles (Dirichlet MC). | `seed` y el patrón `.mbm_save_seed()`. |
| Bajo | S | `aldex_heatmap_plot.R:184-196`, `aldex_volcano_plot.R:131-137` | `case_when` con regex: todas las ramas se evalúan sobre todas las filas (solo las filtradas, pocas). | Helper común. |

### 1.8 `ancombc_plot` (`R/ancombc_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | S | 177 | `pseudo_sens = TRUE` fijo: `ancombc2` repite el ajuste para varios pseudoconteos (análisis de sensibilidad). Es el ajuste más costoso de ANCOM-BC2. | Parámetro `pseudo_sens = FALSE` por defecto (o `TRUE` solo en modo "publicación"). |
| **Alto** | S | 170-184 | No se expone `n_cl` (clusters paralelos de `ancombc2`) ni `iter_control`/`em_control`. | `n_cl = 1` como argumento y pasarlo; documentar `n_cl = parallel::detectCores() - 1`. |
| Medio | S | 147-151 | Construye un objeto **phyloseq** (≈ 8.7 s la primera vez que se carga el namespace), que `ancombc2` convierte internamente a TreeSummarizedExperiment. | Construir directamente un `TreeSummarizedExperiment` (ya es dependencia de ANCOMBC) y quitar phyloseq de Imports. |
| Medio | S | 152-154 | `apply(otu_table, 1, var)` por fila, solo para detectar filas constantes. | `m <- as.matrix(otumat); keep <- rowSums(m != m[, 1]) > 0`, o `matrixStats::rowVars`. |
| Bajo | S | 138 / `zzz.R:51` | `.mbm_parse_taxonomy`: `apply(taxonomy, 2, function(x) ifelse(...))` convierte a matriz de carácter y aplica `ifelse` por columna. | `taxonomy[taxonomy == ""] <- NA` sobre el data.frame. |
| Bajo | S | 406-407 | `apply(res_prim[, diff_t], 1, function(r) any(r %in% TRUE))`. | `sum(rowSums(as.matrix(res_prim[, diff_t]) == TRUE, na.rm = TRUE) > 0)`. |
| Bajo | S | 106-114 | `requireNamespace("microbiome")`: comprobar si `ancombc2` ≥ 2.x todavía lo necesita (hoy usa mia/TreeSummarizedExperiment). | Si no hace falta, quitar la comprobación (ahorra una instalación pesada). |

### 1.9 `beta_ord_plot` (`R/beta_ord_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | S | 155-158 | El mismo desperdicio de ALDEx2 que en 1.1 (128 MC, usa 1). Es la opción **por defecto** (`distance = "compositional"`). | Snippet A. |
| Medio | S | 185 | `metaMDS(trymax = 100)` sin `set.seed` (la semilla solo se fija en la rama compositional) y sin `parallel`. | Mover `set.seed(seed)` arriba y exponer `trymax`/`parallel`. Para las distancias precalculadas, `autotransform = FALSE`. |
| Bajo | S | 177 | `prcomp()` calcula todas las componentes (79) para usar 2. | `prcomp(pca_input, rank. = 2)`: el SVD sigue siendo completo, pero con n = 79 es barato. Con n grande, usar `irlba::prcomp_irlba(n = 2)`. |
| Bajo | S | 137 | `as.data.frame(lapply(abund_table, as.numeric))` columna a columna. | `m <- as.matrix(abund_table); storage.mode(m) <- "double"`. |
| Bajo | S | 311-391 | `extract_clean_label` se define en cada llamada y recorre fila a fila (solo top_n = 5, coste despreciable). | Mover a un helper interno. |

### 1.10 `beta_turnover_plot` (`R/beta_diversity.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | M | 146-163 | `hillR::hill_taxa_parti_pairwise(q)` × 3: bucle R sobre todos los pares, llamando `hill_taxa_parti` por par (O(n²·S) en R interpretado, con barra de progreso y `capture.output`). Ya se restringe a las muestras necesarias (bien). | Snippet E: versión matricial, equivalente a hillR con `rel_then_pool = TRUE` (validar con un test contra hillR). |
| Medio | S | 116-121 | `dplyr::filter(rowSums(across(where(is.numeric))) != 0)` y `as.data.frame(t())`. | `m <- as.matrix(table); m <- m[rowSums(m) > 0, ]; x <- t(m)`. |
| Bajo | S | 171-173 | Doble `inner_join` del resultado largo con la metadata completa (todas las columnas × 2). | Unir solo las columnas `condition1`/`condition2` con `match()`. |

**Snippet E: beta de Hill por pares vectorizada**
```r
.mbm_hill_beta_pairs <- function(x, q) {          # x: muestras x features (conteos)
  p <- x / rowSums(x); n <- nrow(p)
  if (q == 0) {
    b <- (p > 0) + 0; S <- rowSums(b)
    alpha <- outer(S, S, "+") / 2
    gamma <- outer(S, S, "+") - tcrossprod(b)            # |union|
  } else if (q == 2) {
    G <- tcrossprod(p); d <- diag(G)
    alpha <- 2 / outer(d, d, "+")
    gamma <- 4 / (outer(d, d, "+") + 2 * G)
  } else {                                               # q == 1
    xlx <- function(v) ifelse(v > 0, v * log(v), 0)
    H <- -rowSums(xlx(p))
    alpha <- exp(outer(H, H, "+") / 2)
    gamma <- matrix(NA_real_, n, n)
    for (i in seq_len(n - 1)) {                          # bucle sobre muestras, vectorizado en pares
      j <- (i + 1):n
      m <- (matrix(p[i, ], length(j), ncol(p), byrow = TRUE) + p[j, , drop = FALSE]) / 2
      gamma[i, j] <- exp(-rowSums(xlx(m)))
    }
  }
  beta <- gamma / alpha
  idx <- which(upper.tri(beta), arr.ind = TRUE)
  data.frame(site1 = rownames(x)[idx[, 1]], site2 = rownames(x)[idx[, 2]],
             TD_beta = beta[idx], q = q)
}
```

### 1.11 `beta_dissimilarity_plot` (`R/beta_dissimilarity_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| Medio | S | 89-92 | `table[table > 0] <- 1` sobre un **data.frame** (crea una matriz lógica completa y reasigna por columnas), y luego `as.data.frame(t())` + `lapply(as.numeric)`. | `m <- as.matrix(table) > 0; storage.mode(m) <- "double"; m <- t(m[rowSums(m) > 1, ])`. |
| Medio | S | 107-120 | `reshape2::melt` de la matriz **n × n completa** (ambos triángulos y la diagonal), `left_join` × 2 con toda la metadata y `distinct()`. Hace el doble de trabajo y, además, es un bug (ver B4). | Índices `lower.tri` (patrón ya usado en `beta_decay_plot.R:187-194`). |
| Bajo | S | 96 | `betapart.core()` completo solo para leer `$shared`. | `tcrossprod(m)` directamente. |

### 1.12 `beta_partition_ord_plot` (`R/beta_partition_ord_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| Medio | S | 89-94 | Pipeline tibble: `as_tibble(rownames)` → `arrange` (ordena **features** alfabéticamente, sin efecto útil) → `column_to_rownames` → `select_if` → `t()` → data.frame. | `pa <- t(as.matrix(table) > 0) * 1`. |
| Medio | S | 120-121, 186-188 | `env1` = **toda la matriz de presencia/ausencia (79 × 6 045)** + metadata, solo para alinear la metadata. Luego, en cada uno de los 3 plots, se vuelve a hacer `inner_join` con esa tabla gigante. | `env1 <- metadata[match(rownames(pa), metadata$SampleID), ]`. |
| Medio | S | 124-126 | 3 × `betadisper` (3 PCoA de n × n). Es necesario, pero se puede sacar `beta.pair` una vez (ya se hace). | Aceptable. |
| Bajo | S | 253-274 | `get_legend` más 3 plots con `plot_grid(align = "hv")`, que renderiza dentro de la función. | Aceptable (3 paneles), o usar patchwork. |
| Bajo | S | 67-277 | Toda la función va envuelta en `suppressWarnings({...})`, lo que oculta avisos reales (ver B14). | Quitarlo y silenciar solo la llamada que lo necesite. |

### 1.13 `beta_decay_plot` (`R/beta_decay_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| Medio | S | 231-232 | `vegan::mantel(permutations = 999)` una vez por grupo, sin `parallel` ni semilla. | Exponer `parallel` y `seed`. |
| Bajo | S | 172-194 | Bien vectorizado (`lower.tri`, `arr.ind`). | Es el patrón a copiar en `beta_dissimilarity_plot`. |
| Bajo | S | 141-142 | `requireNamespace()` sin comprobar el resultado (no hace nada útil si falta el paquete). | `if (!requireNamespace("geosphere", quietly = TRUE)) stop(...)`. |

### 1.14 `corr_env_abund_plot` (`R/corr_env_abund_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | S | 451-458 | `for (env) for (taxón) cor.test(...)`: E × T llamadas (5 × 1 000 = 5 000 `cor.test`, cada una con validaciones y `suppressWarnings`). | Snippet F. |
| **Alto** | S | 463 | `cor()` completa **recalculada** tras filtrar. | `corr_mat <- corr_mat[, signif_taxa, drop = FALSE]`. |
| Medio | M | 423, 480 | Sin `top_n`: con `pval_threshold = NULL` se correlacionan y dibujan **todas** las taxas colapsadas (miles de tiles y un `melt` enorme). | Argumento `top_n` (por abundancia media) antes de `cor()`. |
| Medio | S | 130-142, 181-204 | 12 `filter()` y el colapso por `strsplit` por fila **duplicados** de `collapse_table`. | Llamar al helper común (sección 3) o a `collapse_table()`. |
| Medio | S | 418-420 | `apply(env, 2, sd)`, `apply(abund, 1, sd)` y `sweep` sobre data.frame. | `m <- as.matrix(counts); ra <- sweep(m, 2, colSums(m), "/") * 100; ra <- ra[matrixStats::rowSds(ra) > 0, ]`, o bien `rowSums((ra - rowMeans(ra))^2) > 0`. |

**Snippet F: p-valores vectorizados**
```r
r  <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
n  <- crossprod(!is.na(env), !is.na(t(abund)))           # n por par (o nrow(env) si no hay NA)
tt <- r * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.eps))
pval_mat <- 2 * stats::pt(-abs(tt), df = n - 2)
```
(Pearson: idéntico a `cor.test`. Spearman: coincide con la aproximación t que `cor.test` ya usa cuando hay empates, que es lo habitual con conteos de microbioma; con n < 1290 y sin empates, `cor.test` usa el test exacto y pueden aparecer diferencias pequeñas.)

### 1.15 `random_forest_lollipop_plot` (`R/random_forest_lollipop_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| **Alto** | S | 78-83 | `importance = TRUE` calcula la importancia por permutación (OOB × p features × ntree). Solo se usa `MeanDecreaseGini`, que se calcula siempre. | `importance = FALSE`. |
| Medio | S | 78-83 | `ntree = 500` fijo, `x` es un data.frame con miles de features (muchas todo ceros tras filtrar muestras) y no hay semilla. | Quitar las columnas con varianza 0 y features raras. Pasar `x = as.matrix(...)`. Exponer `ntree` y `seed`. |
| Bajo | S | 111 | `across(everything(), trimws)` convierte las columnas numéricas a carácter (luego `as.numeric` en la línea 173). | `trimws` solo en las columnas de rango. |

### 1.16 `ratios_bubble_plot` (`R/ratios_bubble_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| Medio | S | 195-207 | `dplyr::rowwise()` sobre **todas** las taxas antes de `slice_max(top_n)`. `rowwise` es de 100 a 1 000 veces más lento que la versión vectorizada. | `mutate(MeanAbund = (.data[[A]] + .data[[B]]) / 2, Ratio = pmax(A,B)/pmin(A,B) - 1, Dominant = ifelse(A > B, A_name, B_name))` sin `rowwise`. |
| Medio | S | 88-170 | `case_when` anidado: las **tres** ramas internas (con 10+ `str_extract`) se evalúan sobre todas las filas aunque solo aplique una combinación db/level. | `if/else` fuera de `mutate` (o el helper común). |
| Bajo | S | 175-184 | `sweep` + `pivot_longer` + `group_by/summarise` para medias por condición. | `rowMeans(ra[, samples_A])` y `rowMeans(ra[, samples_B])` sobre la matriz. |

### 1.17 `abundance_sankey_plot` (`R/abundance_sankey_plot.R`)

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| Medio | S | 165-171 | `get_level_data()` × 7: cada una hace `unite`, `group_by/summarise` de **todas** las columnas de muestra, `t()`, `as.data.frame`, `colMeans`… sobre la tabla completa. | Calcular `ab <- rowMeans(ra)` **una vez** (la media de sumas es igual a la suma de medias) y luego `tapply(ab, lineage_k, sum)` por nivel: vector, no matriz. |
| Medio | S | 246 | `networkD3::saveNetwork()` **siempre** escribe HTML (selfcontained usa pandoc, que es lento) como efecto secundario. | `output_file = NULL` por defecto y guardar solo si se indica. |
| Bajo | S | 192-195 | `vapply` por link con `my_report$name == target` (O(L × N)). | `match(links$target, my_report$name)`. |

### 1.18 `alpha_hill_corr_plot` / `alpha_decay_plot` / `venn_plot` / `cca_rda_biplot`

| Impacto | Esfuerzo | Línea | Problema | Arreglo |
|---|---|---|---|---|
| Bajo | S | `alpha_hill_corr_plot.R:65-74`; `alpha_decay_plot.R:141-150` | `data.frame(t())` y 3 × `hill_taxa`. | Snippet C. |
| Bajo | S | `alpha_hill_corr_plot.R:99-196`; `alpha_decay_plot.R:366-385` | 3 paneles con `cowplot::plot_grid`: render dentro de la función (coste moderado con 3 paneles). | `patchwork`, que es perezoso, si se quiere devolver un objeto componible. |
| Bajo | S | `venn_plot.R:98-108` | Cada grupo subsetea el data.frame; `subset > 0` crea una copia lógica. | `m <- as.matrix(table)` una vez; `rowMeans(m[, s, drop = FALSE] > 0)`. |
| Bajo | S | `cca_rda_biplot.R:168` | `envfit` con 999 permutaciones fijas (no se exponen) ni semilla propia (usa la global ya fijada, bien). | Exponer `permutations`. |

---

## 2. Diseño de API: re-parseo en cada función

Todas las funciones exportadas reciben `table` (data.frame features × muestras con columna `taxonomy`) y `metadata`, y
**cada una** vuelve a:

1. detectar la columna de taxonomía con `grep()` (a veces asumiendo después que es la *última*, ver B8);
2. alinear muestras (`intersect`, `match`), unas veces con la 1ª columna de metadata y otras con `rownames`;
3. convertir data.frame ↔ matriz, transponer, y calcular abundancia relativa;
4. parsear o colapsar la taxonomía con regex (el mismo `case_when` de 5 a 7 ramas aparece en **6 archivos**:
   `abundance_bar_plot`, `abundance_heatmap_plot`, `aldex_heatmap_plot`, `corr_env_abund_plot`,
   `random_forest_lollipop_plot`, `ratios_bubble_plot`);
5. calcular distancias, CLR o diversidad que otra función ya había calculado (p. ej. `beta_ord_plot` +
   `beta_test_table` con `"compositional"` hacen **dos** `aldex.clr` de 128 MC sobre los mismos datos).

**Propuesta: objeto preprocesado compartido (compatible hacia atrás)**
```r
mbm_data <- function(table, metadata, taxonomy = NULL, tax_sep = ";\\s*") {
  tax_col <- grep("^(taxonomy|taxon|taxa)$", names(table), ignore.case = TRUE)
  counts  <- as.matrix(table[, setdiff(seq_along(table), tax_col), drop = FALSE])
  ids     <- intersect(as.character(metadata[[1]]), colnames(counts))
  counts  <- counts[, ids, drop = FALSE]
  tax     <- .mbm_split_tax(table[[tax_col]], tax_sep)   # matriz n x 7, parseada UNA vez
  obj <- list(counts = counts, tax = tax, raw_tax = table[[tax_col]],
              meta = metadata[match(ids, metadata[[1]]), , drop = FALSE],
              cache = new.env(parent = emptyenv()))       # distancias, CLR, rel.ab., hill...
  class(obj) <- "mbm_data"; obj
}
.mbm_get <- function(obj, key, fun) {                     # memoización por objeto
  if (is.null(obj$cache[[key]])) obj$cache[[key]] <- fun()
  obj$cache[[key]]
}
# en cada función exportada:
x <- if (inherits(table, "mbm_data")) table else mbm_data(table, metadata)
d <- .mbm_get(x, paste0("dist_", method), function() vegan::vegdist(t(x$counts), method))
```
Ventajas: el parseo de taxonomía y la alineación se hacen una sola vez, y las distancias, el CLR, `aldex()` y Hill se reutilizan entre
`beta_ord_plot`, `beta_test_table`, `aldex_*`, etc. Añadir métodos `as_mbm_data.phyloseq()` y
`as_mbm_data.TreeSummarizedExperiment()` (extrayendo `otu_table`/`assay`, `tax_table`/`rowData`, `sample_data`/`colData`),
con esos paquetes en `Suggests`, permite a usuarios de phyloseq y mia pasar sus objetos directamente.

---

## 3. Helpers compartidos recomendados (eliminan duplicación y coste)

| Helper | Sustituye | Idea |
|---|---|---|
| `.mbm_split_tax(tax)` | `separate()` en heatmap, sankey, RF y `.mbm_parse_taxonomy` | `do.call(rbind, lapply(strsplit(tax, ";\\s*"), \`length<-\`, 7))`: una sola división y una matriz n × 7. |
| `.mbm_tax_label(tax, level, db)` | Unos 15 bloques `case_when` (6 archivos) | Bucle **sobre los 7 rangos** (no sobre filas): se toma el rango pedido; si no es válido (`""`, `uncultured`, `Incertae_Sedis`), se prueba el superior con el prefijo `"other "`. Todo vectorizado con `ifelse`/`is.na`. |
| `.mbm_collapse(counts, tax, level)` | `collapse_table`, `corr_env_abund_plot` 172-215, `abundance_bar_plot` 145-196, `ratios_bubble_plot` 83-85 | `rowsum(counts, .mbm_truncate_tax(tax, level), reorder = FALSE)`. |
| `.mbm_uninformative_tax` | Los 12 `filter()` × 2 archivos | Un vector constante más `%in%`. |
| `.mbm_hill(x)` / `.mbm_hill_beta_pairs(x, q)` | `hillR::*` (4 archivos) | Snippets C y E. `hillR` puede pasar a Suggests (solo se usaría en tests de equivalencia). |
| `.mbm_clr(x)` | `aldex.clr(..., 128)` + 1 réplica (2 archivos) | Snippet A. |
| `.mbm_alpha_panels()` | Unas 400 líneas duplicadas `alpha_diversity_plot`/`alpha_hill_plot` | Snippet B. |

```r
.mbm_tax_label <- function(tax, level_idx = 6, bad = "uncultured|Incertae_Sedis|^$") {
  m <- sub("^[a-z]{1,2}__", "", .mbm_split_tax(tax))           # n x 7 sin prefijos
  ok <- !is.na(m) & !grepl(bad, m, ignore.case = TRUE)
  out <- ifelse(ok[, level_idx], m[, level_idx], NA_character_)
  for (k in seq(level_idx - 1, 2)) {                            # bucle sobre rangos (≤ 6 iteraciones)
    fill <- is.na(out) & ok[, k]
    out[fill] <- paste("other", m[fill, k])
  }
  out[is.na(out)] <- "Unclassified"
  out
}
```

---

## 4. Peso de dependencias (`DESCRIPTION`/`NAMESPACE`) y `zzz.R`

Tiempo de `requireNamespace()` en frío (R 4.6.1, un proceso nuevo por paquete; incluye las dependencias que arrastra cada uno):

| Paquete | Carga (s) | Dónde se usa | Recomendación |
|---|---|---|---|
| phyloseq | **8.72** | solo `ancombc_plot.R:147-154` | **Suggests**: construir un TSE directamente, o aceptar phyloseq como *entrada* opcional. |
| ALDEx2 | **5.40** | `aldex_*`, `beta_ord_plot`, `beta_test_table` (compositional) | **Suggests** + `requireNamespace()`. En `beta_*`, sustituir por `.mbm_clr()`. |
| ANCOMBC | 2.58 (instalación enorme: mia, lme4, CVXR, TreeSummarizedExperiment…) | solo `ancombc_plot` | **Suggests** (ya hace `requireNamespace`). |
| ComplexHeatmap (+ circlize) | 1.52 (+ 0.11) | 2 heatmaps | **Suggests**. Ojo: el valor por defecto de `aldex_heatmap_plot(effect_colors = circlize::colorRamp2(...))` se evalúa al llamar la función; cambiarlo a `NULL` y resolverlo dentro. |
| betapart | 0.97 | `beta_dissimilarity_plot`, `beta_partition_ord_plot` | Suggests, o reimplementar `beta.pair` (son 20 líneas con `tcrossprod`). |
| ggvenn / ggVennDiagram | 0.94 / 0.95 | `venn_plot` | **Suggests**. Además, `NAMESPACE` tiene `importFrom(ggvenn, ggvenn)` e `importFrom(ggVennDiagram, ggVennDiagram)`, que **se cargan en `library(MicroBioMeta)`** aunque nunca se use Venn. Quitar esos `importFrom` (el código ya usa `::`). |
| randomForest | 0.01 | `random_forest_lollipop_plot` | Suggests. Quitar `importFrom(randomForest, importance, randomForest)` de NAMESPACE. |
| geosphere | 0.16 | `beta_decay_plot` | Suggests (o haversine en 5 líneas de base R). |
| hillR | 0.00 | 4 archivos | Reemplazar por los snippets C/E y pasarlo a Suggests. |
| viridis | 0.75 | paletas | `viridisLite::viridis` (ligero, ya lo trae ggplot2) o `grDevices::hcl.colors(n, "Plasma")`. |
| reshape2 | 0.27 | 2 `melt` | Reemplazar por el índice `lower.tri` / `as.data.frame(as.table(m))` y quitarlo. |
| stringr | 0.17 | `str_extract`/`str_trim` | Base R (`regmatches`, `trimws`); ya se mezclan ambos estilos. |
| purrr, tidyselect | — | 3 usos | `lapply`/`Filter(Negate(is.null), ...)`; `tidyselect::starts_with` → `dplyr::starts_with`. |
| tools | — | 1 uso | Base R; se puede quitar de Imports (es un paquete base, pero se declara innecesariamente). |

Observaciones:

- `library(MicroBioMeta)` tarda hoy ≈ 0.85 s (el paquete instalado carga dplyr, ggplot2, ggvenn, ggVennDiagram y randomForest
  por los `importFrom`). Sin los `importFrom` de Venn y RF, baja a unos 0.4 s.
- El **coste real** está en: (a) la **instalación** (ANCOMBC y phyloseq arrastran cientos de MB de Bioconductor, un obstáculo
  para el público objetivo "beginners"); y (b) la **primera llamada** de cada función, que paga la carga del namespace
  (p. ej. unos 5.4 s de ALDEx2 dentro del primer `beta_test_table(method = "compositional")`). Al comparar con rbiom
  conviene separar `first` de `median`.
- `zzz.R` **no** define `.onLoad`/`.onAttach`, así que no añade coste de carga. Solo tiene dos `utils::globalVariables()`
  duplicados (líneas 24-27 y 282-393), que conviene unificar (cosmético).
- `Depends: R (>= 4.5.0)` es muy restrictivo. Para Bioconductor está bien, pero limita la adopción.
- `importFrom(dplyr, "%>%")` carga dplyr al cargar el paquete. `magrittr` es más ligero; también se puede migrar a `|>`.

Patrón para Suggests:
```r
.mbm_need <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required for this function. Install it with ",
         if (bioc) sprintf('BiocManager::install("%s")', pkg) else sprintf('install.packages("%s")', pkg),
         call. = FALSE)
}
```

---

## 5. Estadística: permutaciones, semillas, paralelismo

| Función | Línea | Estado actual | Recomendación |
|---|---|---|---|
| `beta_test_table` | 82, 160-174 | 999 perm.; semilla **solo** en compositional | Semilla global al inicio; `parallel`. Un solo `adonis2` con `by = "terms"` (correcto, **no** hay llamadas múltiples). |
| `beta_decay_plot` | 118, 231 | `mantel` 999 × grupos, sin semilla | `seed`, `parallel`. |
| `cca_rda_biplot` | 158, 168 | semilla OK; `envfit` 999 fijo | Exponer `permutations`. |
| `beta_ord_plot` | 185 | `metaMDS(trymax = 100)` sin semilla fuera de compositional | `set.seed` global y `parallel`. |
| `aldex_*` | 141 / 120 | `mc.samples = 128` fijo, sin semilla ni `useMC` | Exponer `mc.samples`, `useMC` y `seed`, y pre-filtrar. |
| `ancombc_plot` | 170-184 | `pseudo_sens = TRUE`, `n_cl` no expuesto, `lib_cut = 1000` fijo | `pseudo_sens = FALSE`, `n_cl`, `lib_cut` como argumentos. |
| `random_forest_lollipop_plot` | 78-83 | `importance = TRUE`, sin semilla | `importance = FALSE`, `seed`, `ntree`. |
| `alpha_diversity_plot` | 141 | `rrarefy` sin semilla | `seed`. |
| Tests pareados | — | No hay bucles de tests por pares: `stat_compare_means` hace **un** test global por panel (barato). | OK. |

---

## 6. Plots: render dentro de las funciones

| Función | Línea | Qué renderiza dentro | Recomendación |
|---|---|---|---|
| `abundance_heatmap_plot` | 384-396 | `grid.grabExpr(draw())` **y** `grid.draw()`: 2 renders | Devolver el objeto `Heatmap` (render al imprimir). |
| `aldex_heatmap_plot` | 375-386 | `ComplexHeatmap::draw()` | Igual. |
| `alpha_diversity_plot` / `alpha_hill_plot` | 444-481 / 453-491; 554-608 / 572-630 | `get_legend` + `plot_grid` (N + 1 builds) o `ggplotGrob` | Snippet B: devolver un ggplot. |
| `alpha_decay_plot`, `alpha_hill_corr_plot`, `beta_partition_ord_plot` | 366-404, 170-216, 253-276 | `plot_grid` de 3 paneles más la leyenda | Aceptable; patchwork si se quiere que sea perezoso. |
| `abundance_sankey_plot` | 246 | `saveNetwork()` (disco + pandoc) siempre | Guardar solo si `output_file` no es `NULL`. |
| Todas | — | `save_table = FALSE` por defecto | Correcto. |

Nota sobre el benchmark: phyloseq/rbiom devuelven ggplots **sin renderizar**, mientras que MicroBioMeta renderiza (y en el
heatmap lo hace dos veces) dentro de `system.time()`. Parte de la diferencia es metodológica, pero el doble render es
coste real para el usuario.

---

## 7. Bugs de corrección detectados de paso (breve)

| ID | Archivo:línea | Bug | Arreglo |
|---|---|---|---|
| B1 | `beta_test_table.R:97-98, 152-156` | Si `table` es un `dist` o una matriz de distancias, se convierte a matriz y **se vuelve a pasar por `vegdist`**. Resultado: distancia *entre filas de la matriz de distancias* (el ejemplo con `betadisper` de la documentación da resultados erróneos sin avisar). | `if (inherits(orig, "dist") \|\| (is.matrix(table) && isSymmetric(table))) dist_matrix <- as.dist(table) else ...`. |
| B2 | `beta_test_table.R:115-120` | No alinea las filas de `metadata` con los nombres de muestra: solo compara dimensiones. Si el orden difiere, PERMANOVA asigna grupos a muestras equivocadas. | `metadata <- metadata[match(rownames(table), metadata[[1]]), ]` y validar que no queden `NA`. |
| B3 | `random_forest_lollipop_plot.R:63-70` | `otu_filtered` sigue el orden de la tabla y `metadata_filtered` el de la metadata: la **respuesta queda desalineada** con las muestras si los órdenes difieren. | `metadata_filtered <- metadata[match(common_samples, metadata[[1]]), ]`. |
| B4 | `beta_dissimilarity_plot.R:107-120` | `melt` de la matriz simétrica completa: cada par aparece **dos veces** (A-B y B-A; `distinct(site1, site2)` no lo colapsa). Para `"shared"` se incluye la **diagonal** (auto-comparaciones = riqueza), y `value != 0` descarta ceros legítimos. Esto infla n y sesga medias y p-valores. | Solo `lower.tri`, sin el filtro `!= 0`. |
| B5 | `beta_dissimilarity_plot.R:152` | `if (!is.null(names(group_colors) == FALSE))` siempre es `TRUE`, así que **sobrescribe** los nombres de colores que haya dado el usuario. | `if (is.null(names(group_colors)))`. |
| B6 | `beta_ord_plot.R:421-423` | `coord_cartesian(xlim, ylim)` seguido de `coord_fixed()`: el segundo **reemplaza** al primero y se pierden los límites simétricos. | `ggplot2::coord_fixed(xlim = c(-r, r), ylim = c(-r, r))`. |
| B7 | `cca_rda_biplot.R:92-93` | `rownames(metadata) <- metadata$SampleID` se ejecuta aunque `metadata = NULL` (el valor por defecto), lo que produce un error. | Moverlo dentro del `if (!is.null(metadata))`. |
| B8 | `abundance_heatmap_plot.R:95`, `aldex_volcano_plot.R:96-104`, `beta_ord_plot.R:126-128`, `cca_rda_biplot.R:83-84`, `random_forest_lollipop_plot.R:53-55` | Detectan la columna de taxonomía con `grep` pero luego **asumen que es la última**. | Usar el índice de `grep` (`tax_col`). |
| B9 | `aldex_heatmap_plot.R:130, 180-181, 263`; `aldex_volcano_plot.R:116-117, 158-159, 263-265` | ALDEx2 calcula `diff.btw`/`effect` según los niveles de factor **ordenados** de `conditions`, pero las etiquetas "Higher/Lower in" y los colores usan el orden de `unique()`. El sentido puede quedar **invertido** si los niveles no están en orden alfabético. (Hay que verificarlo con la versión instalada de ALDEx2.) | `lv <- levels(factor(conditions))` y etiquetar respecto de `lv[1]`/`lv[2]`. |
| B10 | `aldex_volcano_plot.R:115-120` | No subsetea a `cond` y `other_cond`: con más de 2 grupos, `aldex(test = "t")` falla o compara mal. | Filtrar las muestras a esos dos grupos. |
| B11 | `alpha_hill_plot.R:158` | `nlevels()` sobre una columna no factor siempre da 0, así que la paleta de 2 grupos nunca se aplica. | `length(unique(results[[fill_col]])) == 2`. |
| B12 | `ancombc_plot.R:384-386` | `grep("^lfc_Location")` también captura `lfc_Location...:Treatment...` (interacciones), lo que puede mandar un término simple a la rama de heatmap. | Excluir las columnas que contienen `":"` en los términos simples. |
| B13 | `abundance_heatmap_plot.R:143` | `warning()` de renombrado de filos en **cada** llamada. | Condicional. |
| B14 | `beta_partition_ord_plot.R:67, 120-126` | `suppressWarnings` global. Si alguna muestra de la tabla falta en la metadata, `inner_join` la elimina de `env1` pero no de `jac`, y `betadisper` recibe grupos de otra longitud. | Alinear o subsetear la matriz de distancias con `match()`. |
| B15 | `corr_env_abund_plot.R:480-492` | Tras colapsar, varias filas comparten etiqueta (`"other X"`, `"Unclassified"`): los tiles **se solapan** en el mismo eje Y. | `make.unique()` o agrupar por etiqueta antes de `cor()`. |
| B16 | `ratios_bubble_plot.R:197-201` | `Ratio` da `Inf`/`NaN` si una media es 0. Además, `level == "specie"` no coincide con el `"species"` que usa el resto del paquete. | `pmax(A, B) / pmax(pmin(A, B), eps) - 1`; aceptar `"species"`. |
| B17 | `beta_decay_plot.R:172` | `vegdist(method = "jaccard")` sin `binary = TRUE` es Jaccard **cuantitativo**, no presencia/ausencia (lo que suele esperarse en distance-decay). | Documentarlo o añadir `binary = distance %in% c("jaccard", "sorensen")`. |
| B18 | Varias (sección 5) | Resultados no reproducibles: ALDEx2 (heatmap/volcano), `metaMDS`, `randomForest`, `mantel`, `adonis2` fuera de compositional, `rrarefy`. | Argumento `seed` y `.mbm_save_seed()` en todas. |
| B19 | `alpha_diversity_plot.R:3` | La documentación habla de "Hill numbers (q = 0, 1, 2)", pero se calcula Chao1/Shannon/Simpson. | Corregir la documentación. |

---

## 8. Orden de trabajo sugerido

1. **Día 1 (quick wins 1, 2, 5, 6, 7 y B1-B3):** CLR directo, heatmap sin doble render, p-valores vectorizados,
   `pseudo_sens`/`importance`, y los bugs de alineación. Previsiblemente: `beta_test_table` pasa de unos 15 s a < 2 s,
   `abundance_heatmap_plot` baja alrededor de un 50 %.
2. **Día 2 (quick win 3 y snippet B):** reescribir la composición de `alpha_diversity_plot`/`alpha_hill_plot` como un
   único ggplot. Previsiblemente: de 5.4 s a < 1 s en la llamada (el render se paga una sola vez al imprimir).
3. **Días 3-4:** helpers compartidos (sección 3), Hill vectorizado y `beta_turnover_plot`, con tests de equivalencia
   contra la salida actual y contra hillR/cor.test.
4. **Semana 2:** mover dependencias a Suggests (sección 4) y el objeto `mbm_data` con caché y métodos para phyloseq/TSE (sección 2).
