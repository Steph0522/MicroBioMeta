# Perfilado de rendimiento de MicroBioMeta frente a phyloseq y rbiom

**Entorno:** Windows 10, R 4.6.1, 6 núcleos. MicroBioMeta 0.99.0 cargado con `devtools::load_all()` (solo lectura). phyloseq (instalado) y rbiom 3.1.0 (instalado en la librería de usuario para este estudio).
**Datos:** `inst/extdata` con `Month == "2"`, preparados con el código de la colaboradora copiado literalmente (`00_setup_data.R`).
**Tamaño:** **79 muestras × 6045 features**; 3645 features tienen conteo > 0 en esas muestras. Hay 3 niveles de Treatment (TD, TED, TC) y 4 de Type_of_soil.

No se modificó ningún archivo del paquete. Todo lo generado está en `optimization/profiling/`. Para reproducirlo, ejecutar cada script desde esa carpeta con `Rscript <script>.R`.

> **Advertencias de método**
> - En las tablas, "1ª" es la primera ejecución en la sesión, que incluye carga perezosa y JIT. "Mediana" es la mediana de las 3 repeticiones siguientes (5 en algunos prototipos). Los tiempos son `elapsed` de `system.time`.
> - **rbiom guarda sus resultados en una caché en disco** (`tempdir()/rbiom/cache`). Si se repite la misma llamada, devuelve el resultado guardado en 0.01–0.2 s. Para medir de forma justa se usa un directorio de caché nuevo en cada repetición (`rbiom_fresh_cache()`).
> - En Windows, `Rprof` registra menos muestras de las reales (el tiempo total muestreado es aproximadamente 1/3 del `elapsed`). Por eso los hotspots se expresan en **porcentaje del tiempo muestreado**, no en segundos.

---

## 1. Reproducción de los tiempos originales (misma metodología: `system.time` de la llamada, sin renderizar)

Script: `01_original_timings.R`. Resultados en `results_01_original_timings.csv` (caché de rbiom limpia) y `results_01_original_timings_rbiom_cache_on.csv` (caché de rbiom activa).

| Análisis | Llamada | Colaboradora (s) | Reproducido 1ª (s) | Reproducido mediana (s) | rbiom con caché activa, mediana (s) |
|---|---|---:|---:|---:|---:|
| Barras | MBM `abundance_bar_plot` | 0.61 | 0.22 | **0.19** | |
| | phyloseq `plot_bar` + preparación | 2.67 | 2.72 | 1.67 | |
| | rbiom `taxa_stacked` | 0.25 | 1.31 | 1.30 | 0.18 |
| Heatmap | MBM `abundance_heatmap_plot` | 4.75 | 1.27 | **0.87** | |
| | phyloseq `plot_heatmap` | 2.11 | 1.55 | 1.34 | |
| | rbiom `taxa_heatmap` | 0.22 | 1.30 | 1.23 | 0.12 |
| Alfa | MBM `alpha_diversity_plot` | 5.42 | 3.55 | **3.56** | |
| | phyloseq `plot_richness` | 0.23 | 0.20 | 0.20 | |
| | rbiom `adiv_boxplot` | 0.21 | 0.92 | 0.89 | 0.15 |
| Ordenación | MBM `beta_ord_plot` | 0.22 | 0.17 | 0.18 | |
| | phyloseq `ordinate` + `plot_ordination` + `print` | 0.93 | 0.61 | 0.59 | |
| | rbiom `bdiv_ord_plot` (adonis2 999) | 0.24 | 1.32 | 1.33 | 0.14 |
| Estadística beta | MBM `beta_test_table` | 14.95 | 15.09 | **12.30** | |
| | rbiom `bdiv_stats` | 0.02 | 1.53 | 1.50 | 0.02 |
| | rbiom `bdiv_distmat` + `distmat_stats` | 0.19 | 0.11 | 0.11 | 0.01 |
| CCA | MBM `cca_rda_biplot` | 0.50 | 0.33 | 0.43 | |
| | phyloseq `ordinate(CCA)` + plot | 0.27 | 0.19 | 0.18 | |

Tiempo de preparación (`load_all` + librerías + lectura de datos): 12.3 s.

**Lectura de los resultados:**
- `beta_test_table` (≈12–15 s) y `alpha_diversity_plot` (≈3.5 s) **se reproducen**: son lentas de verdad.
- **`abundance_heatmap_plot` no se reproduce en su estado estable.** Tarda 0.7–0.9 s, más rápido que phyloseq (1.3 s) y rbiom (1.2 s). Las dos primeras llamadas de una sesión nueva tardan 1.4–3.0 s: la primera carga ComplexHeatmap de forma perezosa (1.5 s, sección 5) y la segunda también es lenta, por JIT y compilación de lambdas. Además, la función dibuja en el dispositivo activo (el panel Plots de RStudio). Los 4.75 s de la colaboradora corresponden a una llamada "en frío" en RStudio.
- **Los tiempos de rbiom de la colaboradora (0.13–0.25 s) coinciden con aciertos de su caché en disco.** Con la caché limpia, rbiom tarda 0.9–1.5 s por gráfico, igual o más que MicroBioMeta en barras, heatmap y ordenación.

---

## 2. Comparación justa: construir + renderizar

Script: `02_fair_timings.R`, con resultados en `results_02_fair_build_render.csv`.

phyloseq y rbiom devuelven objetos `ggplot` **sin renderizar**; el trabajo de dibujo (`ggplot_build` + `ggplot_gtable`) ocurre al imprimirlos. MicroBioMeta, en cambio, renderiza **dentro** de la función en dos casos:
- `alpha_diversity_plot`: `cowplot::plot_grid` convierte 13 ggplots en gtables.
- `abundance_heatmap_plot`: `grid.grabExpr(draw())` y después `grid.draw` otra vez.

Aquí se mide el tiempo de construcción y el de render en `pdf(NULL)`. Son medianas de 3 repeticiones tras un calentamiento, con la caché de rbiom limpia.

| Análisis | Llamada | Construir (s) | Render (s) | **Total (s)** |
|---|---|---:|---:|---:|
| Barras | MBM `abundance_bar_plot` | 0.22 | 0.30 | **0.52** |
| | phyloseq `plot_bar` (+preparación) | 1.54 | 0.25 | 1.78 |
| | rbiom `taxa_stacked` | 1.21 | 0.33 | 1.54 |
| Heatmap | MBM `abundance_heatmap_plot` (ya dibuja dentro) | 0.70 | 0.03 | **0.73** |
| | phyloseq `plot_heatmap` (+prune/glom) | 1.23 | 0.14 | 1.37 |
| | rbiom `taxa_heatmap` | 1.22 | 0.29 | 1.52 |
| Alfa | MBM `alpha_diversity_plot` | 2.92 | 0.36 | **3.28** |
| | phyloseq `plot_richness` | 0.19 | 0.45 | 0.64 |
| | rbiom `adiv_boxplot` | 0.87 | 0.81 | 1.67 |
| Ordenación | MBM `beta_ord_plot` | 0.18 | 0.17 | **0.36** |
| | phyloseq `ordinate` + `plot_ordination` | 0.44 | 0.18 | 0.62 |
| | rbiom `bdiv_ord_plot` | 1.30 | 0.25 | 1.55 |
| CCA | MBM `cca_rda_biplot` | 0.43 | 0.20 | **0.63** |
| | phyloseq CCA + plot | 0.19 | 0.16 | 0.35 |

Con un dispositivo png activo (que imita a la pantalla), el heatmap de MBM tarda 0.82 s. En una ventana `windows()` real, las dos primeras llamadas tardan 1.85 s y 3.06 s y las siguientes ≈0.85 s.

**Conclusión:** con una medición justa y sin caché, MicroBioMeta es **el más rápido** en barras, heatmap y ordenación. Las diferencias reales se reducen a dos funciones:
- `alpha_diversity_plot`: unas 5 veces más lento que phyloseq y 2 veces más que rbiom en total. Renderizar dentro de la función no explica la diferencia, pero es su causa principal (sección 3).
- `beta_test_table`: la diferencia viene de trabajo extra (ver abajo).

### 2.1 `beta_test_table`: ¿se compara el mismo trabajo?

Resultados en `results_02_beta_decomposition.csv`.

| Variante | Construir (s) | Render tabla (s) |
|---|---:|---:|
| MBM `compositional`, `~Treatment*Type_of_soil`, 999 perm (llamada original) | 12.92 | 0.12 |
| MBM `compositional`, `~Treatment`, 999 | 12.70 | 0.09 |
| MBM `bray`, `~Treatment*Type_of_soil`, 999 | 0.59 | 0.10 |
| **MBM `bray`, `~Treatment`, 999 (mismo trabajo que rbiom)** | **0.30** | 0.08 |
| solo `ALDEx2::aldex.clr(mc.samples = 128)` | **10.91** | – |
| solo `vegan::adonis2` bray `~T*S` 999, `by = "terms"` | 0.32 | – |
| solo `vegan::adonis2` bray `~Treatment` 999 | 0.07 | – |
| rbiom `bdiv_distmat` + `distmat_stats` (bray `~Treatment`, 999) | 0.10 | – |

La comparación original no era equivalente. MicroBioMeta hacía una transformación CLR composicional mediante ALDEx2 con 128 réplicas Monte Carlo de Dirichlet, un modelo con interacción de 3 términos (`by = "terms"`) y además construía una tabla gráfica con ggpubr. rbiom solo calculaba Bray-Curtis con un factor.

- Con el mismo trabajo, MBM tarda 0.30 s frente a 0.10 s de rbiom. Los ≈0.2 s de diferencia son la tabla `ggtexttable`.
- **Las permutaciones no son el problema** (adonis2 ≈0.3 s). El **85–95 % del tiempo es `aldex.clr`**.

---

## 3. Hotspots (Rprof con `line.profiling`)

Script: `03_profile_hotspots.R`. Salidas en `results_03_profile_summary.txt` y `rprof_*.out`. Se usó `keep.source = TRUE` para obtener archivo:línea.

| Función | Archivo:línea | % tiempo muestreado | Qué ocurre |
|---|---|---:|---|
| `beta_test_table` | **R/beta_test_table.R:148** | **92.9 %** | `ALDEx2::aldex.clr(t(table), mc.samples = 128, ...)`. Dentro: `rgamma` 42 % (self), `rdirichlet` 50 %, `apply`/`log2` sobre las 128 columnas ≈35 %, `t()`/`matrix` 54 % (valores inclusivos, se solapan). |
| | R/beta_test_table.R:149 | – | `getMonteCarloSample(aldex_obj, 1)`: **solo se usa la 1ª de las 128 instancias**; el 99 % del trabajo MC se descarta. |
| | R/beta_test_table.R:160 | 2.1 % | `vegan::adonis2(..., by = "terms")` con 999 permutaciones |
| `alpha_diversity_plot` | **R/alpha_diversity_plot.R:451** | **47.1 %** | `cowplot::plot_grid(plotlist = 12 plots)`, que llama a `align_plots`, `as_grob` y `ggplotGrob` sobre cada panel |
| | R/alpha_diversity_plot.R:435 | 20.5 % | `build_cell()` × 12 (4 suelos × 3 índices), cada uno con su propio `ggplot` y theme |
| | R/alpha_diversity_plot.R:332 | 19.3 % | construcción de `ggplot()` dentro de `build_cell` |
| | R/zzz.R:152 | 15.9 % | `.mbm_theme()` (`theme_bw` + `theme`) ejecutado 13 veces |
| | R/alpha_diversity_plot.R:465 | 10.0 % | 2º `cowplot::plot_grid(panel_grid, leg)` |
| | R/alpha_diversity_plot.R:444 | 6.8 % | `get_legend(build_cell(1, 1))`: construye un 13er ggplot solo para la leyenda |
| | R/alpha_diversity_plot.R:138,145–148 | ≈8 % | `data.frame(t(table))` de 6045 columnas + `vegan::estimateR`/`diversity`, ya vectorizados |
| `abundance_heatmap_plot` | **R/abundance_heatmap_plot.R:384** | **45.7 %** | `grid.grabExpr(ComplexHeatmap::draw(...))`: layout y dibujo fuera de pantalla |
| | **R/abundance_heatmap_plot.R:110** (bloque 110–142) | **27.9 %** | `separate` + `case_when` con `grepl`/`str_extract` sobre **las 6045 filas** antes de quedarse con `top_n = 20` |
| | R/abundance_heatmap_plot.R:165 (bloque 165–184) | 8.4 % | `case_when` de 13 ramas vía `across()` en 79 columnas. La lambda se recompila (JIT) en cada llamada: `cmpfun` suma ≈32 % del perfil. |
| | R/abundance_heatmap_plot.R:394–395 | 1.8 % en pdf nulo | `grid.newpage()` + `grid.draw()`: **segundo dibujo** en el dispositivo activo. Es más caro en RStudio. |
| `abundance_bar_plot` | R/abundance_bar_plot.R:447 | 19.0 % | construcción del `ggplot` |
| | R/zzz.R:152 | 14.7 % | `.mbm_theme()` |
| | R/abundance_bar_plot.R:114 (bloque 114–126) | 12.3 % | 12 `dplyr::filter()` encadenados |
| | R/abundance_bar_plot.R:194 | 8.6 % | `group_by(taxonomy)` + `summarise(across(...))` |

---

## 4. Prototipos

Todos son scripts independientes. Generan variantes a partir del código real de `R/` con `patch_fun()`: leen el archivo, aplican una sustitución literal y evalúan el resultado en el namespace del paquete, sin tocar el original. En todos se usa la misma semilla. Se comprueba que el resultado sea idéntico: `identical()` en tablas y distancias, y comparación píxel a píxel de un PNG renderizado.

### 4.1 `prototype_beta_test_table.R` (resultados en `results_prototype_beta_test_table.csv`)

| Variante | Mediana (s) | Aceleración | ¿Resultado idéntico? |
|---|---:|---:|---|
| `beta_test_table()` pública (original) | 12.32 | 1.0× | – |
| Núcleo original (`aldex.clr` 128 MC + dist + adonis2) | 11.72 | – | referencia |
| **P1 exacto:** mismos `rgamma(l*128, a)` por muestra (mismo flujo del RNG), pero solo se calcula la CLR de la 1ª instancia | **6.19–6.32** | **≈1.9×** | **Sí: `identical(dist)` TRUE y `identical(tabla adonis2)` TRUE (R², F y p idénticos)** |
| P2: 1 instancia MC (≈`mc.samples = 1`) | 0.48–0.50 | ≈25× | No: otra realización de Dirichlet. p de la interacción 0.616 → 0.297; correlación entre distancias 0.78. |
| P3: P1 + `adonis2(parallel = 4)` | 7.0–8.4 | más lento | Sí. Crear el clúster en Windows cuesta más de lo que ahorran unas permutaciones que duran 0.3 s. |
| **P4: CLR esperada determinista** (`digamma(a) − media`, que es E[log p] bajo Dirichlet, el límite de promediar infinitas instancias MC) | **0.46** | **≈27×** | No es idéntica a la implementación actual, pero la correlación con la **media de las 128 instancias** de ALDEx2 es 0.99992 y la tabla es prácticamente igual (R² del suelo 0.176 frente a 0.172; p 0.001, 0.535 frente a 0.531, 0.297 frente a 0.298). No depende de la semilla. |

**Hallazgo metodológico importante:** la implementación actual usa **una sola** instancia Monte Carlo, así que el resultado depende de la semilla. Con las semillas 1, 2, 3, 123 y 2024, el p de la interacción fue 0.567, 0.293, 0.618, 0.616 y 0.558, y el p de Treatment varió entre 0.465 y 0.793 (`results_prototype_beta_seed_sensitivity.csv`). Además, esa única instancia añade ruido que **reduce el R² de Type_of_soil** a 0.065, frente a 0.172 con la CLR promedio.

Recomendación: P1 si se exige una salida idéntica bit a bit (≈2× más rápido). P4 si se acepta cambiar el resultado a uno determinista y más estable (≈27× más rápido). Las dos decisiones corresponden a los autores.

### 4.2 `prototype_alpha_diversity_plot.R` (resultados en `results_prototype_alpha_diversity_plot.csv`)

| Variante | Construir (s) | Render (s) | Total (s) | ¿Idéntico? |
|---|---:|---:|---:|---|
| Original (paquete) | 2.94–3.06 | 0.36 | 3.30–3.44 | – |
| A: leyenda tomada de `plots[[1]]` en lugar de `build_cell(1, 1)` (línea 444) | 2.87–3.21 | 0.37 | ≈3.25–3.6 | PNG píxel a píxel idéntico |
| C: A + matriz en lugar de `data.frame(t(table))` (línea 138) | 2.91 | 0.39 | 3.33 | PNG píxel a píxel idéntico |
| **B: un único ggplot facetado** (la ruta `facet_grid2` + etiquetas por gtable que ya existe en la función, forzando `use_grid_compose = FALSE`) | **1.11–1.41** | 0.20 | **1.31–1.63** (≈2.3×) | No: diseño distinto (9 % de píxeles difieren; ver `alpha_*.png`). Los valores de diversidad son los mismos. |
| Solo el cálculo de diversidad (líneas 145–148) | 0.31 | – | – | – |
| Diversidad con matriz en lugar de data.frame | 0.13 | – | – | `all.equal` TRUE |
| phyloseq `plot_richness` (referencia) | 0.19 | 0.46 | 0.65 | – |

A y C son ganancias pequeñas, comparables al ruido de medición (≈0.1–0.2 s). **El coste está en el diseño de 12 ggplots independientes combinados con cowplot.** Solo cambiando a un único ggplot facetado se gana de forma clara (≈2.3×), a cambio de que el aspecto cambie.

### 4.3 `prototype_abundance_heatmap_plot.R` (resultados en `results_prototype_abundance_heatmap_plot.csv`)

| Variante | Construir, mediana (s) | Aceleración | ¿Idéntico? |
|---|---:|---:|---|
| Original, pdf(NULL) | 0.72 | 1.0× | – |
| Original, png activo | 0.86 | – | – |
| H1: `findInterval()` vectorizado en lugar del `case_when` de 13 ramas (líneas 170–184) | 0.61 | 1.2× | Tabla guardada y PNG idénticos |
| **H2: H1 + clasificación sobre taxonomías únicas y `separate` solo en las `top_n` filas (líneas 110–142)** | **0.40** | **1.8×** | Tabla guardada y PNG idénticos. También con `top_n = 200, cluster = FALSE`. |
| H3: H2 + no dibujar dentro (devolver el gTree) | 0.37–0.38 | 1.9× | Objeto y PNG idénticos. Cambia la API: no se autodibuja. |
| phyloseq `plot_heatmap` (referencia) | 1.36 | – | – |

El 46 % restante es `ComplexHeatmap::draw` dentro de `grid.grabExpr`, que es intrínseco si se quiere devolver un grob.

### 4.4 `prototype_abundance_bar_plot.R` (resultados en `results_prototype_abundance_bar_plot.csv`)

| Variante | Construir (s) | ¿Idéntico? |
|---|---:|---|
| Original | 0.19 | – |
| B1: un solo filtro `%in%` en lugar de 12 `dplyr::filter` | 0.16 | `p$data` y PNG idénticos |
| B2: B1 + `rowsum()` en lugar de `group_by`/`summarise` | 0.14 (≈1.35×) | `p$data` y PNG idénticos. También `all.equal` en genus, family y species. |

La función ya era la más rápida de su categoría, así que la mejora es marginal.

---

## 5. Coste de arranque

Script: `05_load_time.R`. Cada medición se hizo en un proceso R nuevo, con la mediana de 3 repeticiones. Resultados en `results_05_load_time.csv`.

| Carga | Mediana (s) | Namespaces cargados |
|---|---:|---:|
| `library(MicroBioMeta)` (instalado 0.99.0) | **0.78** | 33 |
| `devtools::load_all()` (lo que usó la colaboradora) | **8.75** | 198 |
| Todos los Imports con `loadNamespace` | 8.44 | 185 |
| `library(phyloseq)` (competidor) | 7.48 | 75 |
| `library(rbiom)` (competidor) | 2.22 | 77 |

Import más pesado cargado aislado (mediana en s):

| Paquete | s | | Paquete | s |
|---|---:|---|---|---:|
| phyloseq | 7.4–10.5 | | ggrepel | 1.58 |
| **ALDEx2** | **5.35** | | **ComplexHeatmap** | **1.54** |
| **ANCOMBC** | **2.81** | | ggplot2 | 1.50 |
| vegan | 1.47 | | cowplot | 1.19 |
| betapart | 1.05 | | ggpubr / ggh4x | 0.98 |
| randomForest | 0.01 | | | |

- `NAMESPACE` solo importa (`importFrom`) dplyr, ggVennDiagram, ggvenn y randomForest. El resto se usa con `pkg::`, así que **se carga de forma perezosa en la primera llamada**. Por eso `library()` tarda 0.78 s y `load_all()`, que carga todos los Imports, 8.75 s.
- En consecuencia, un usuario con `library()` pagará en su primera llamada:
  - ≈5.4 s adicionales por ALDEx2 en `beta_test_table(method = "compositional")`;
  - ≈1.5 s por ComplexHeatmap en el heatmap (1ª llamada medida: 2.98 s, 2ª: 2.78 s, 3ª: 0.73 s);
  - hasta ≈7–10 s por phyloseq, que solo usa `ancombc_plot`, junto con ANCOMBC.
- Mover phyloseq, ANCOMBC y ALDEx2 a `Suggests` con `requireNamespace()` no aceleraría `library()`, pero sí reduciría el tiempo de instalación y la cadena de dependencias Bioconductor.

---

## Archivos

| Archivo | Contenido |
|---|---|
| `00_setup_data.R` | Preparación de datos literal; utilidades `render_null`, `time_reps`, `patch_fun`, `render_png_pixels`, `rbiom_fresh_cache` |
| `01_original_timings.R` | Reproducción de los tiempos originales (`RBIOM_CACHE=1` activa la caché de rbiom) |
| `02_fair_timings.R` | Construir + renderizar; descomposición de `beta_test_table` |
| `03_profile_hotspots.R` | Rprof por función y línea |
| `05_load_time.R` | Tiempos de carga y de primera llamada |
| `prototype_*.R` | Prototipos con verificación y benchmark |
| `results_*.csv` / `results_*.txt` / `log_*.txt` / `rprof_*.out` / `alpha_*.png` | Salidas |
