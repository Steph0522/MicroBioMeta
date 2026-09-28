# MicroBioMeta: resumen de optimización y posicionamiento

Fecha: 2026-09-25. No se modificó ningún archivo del paquete; todo está en `optimization/`.
Recuerda agregar `^optimization$` a `.Rbuildignore` antes de construir el paquete.

Informes detallados:
- `profiling/PROFILING_REPORT.md`: tiempos reproducidos, comparación justa, zonas de mayor coste, prototipos
- `code_review/CODE_REVIEW.md`: revisión estática de las 25 funciones, errores de resultados
- `competitive_analysis/COMPETITIVE_ANALYSIS.md`: matriz de funciones frente a 11 paquetes, riesgos para Bioconductor, posicionamiento

Datos: 79 muestras × 6045 features (bacteria, Month == 2).

## 1. ¿Qué tan justa fue la comparación?

| Función | Word | Medido de nuevo | Construir + dibujar (MBM / phyloseq / rbiom) | Qué pasó |
|---|---:|---:|---|---|
| `beta_test_table` | 14.95 s | 12.3 s | — | **Real, pero no comparable.** Con el mismo trabajo que rbiom (Bray, `~Treatment`, 999 perm.) tarda 0.30 s frente a 0.10 s. El 93 % del tiempo se va en ALDEx2. |
| `alpha_diversity_plot` | 5.42 s | 3.56 s | 3.28 / 0.64 / 1.67 | **Real.** Arma 12 ggplots con cowplot y los dibuja dentro de la función. |
| `abundance_heatmap_plot` | 4.75 s | 0.72 s | 0.73 / 1.37 / 1.52 | **Primera llamada, sin calentar.** En uso normal, MBM es el más rápido. |
| `abundance_bar_plot` | 0.61 s | 0.52 s | 0.52 / 1.78 / 1.54 | MBM es el más rápido. |
| `beta_ord_plot` | 0.22 s | — | 0.36 / 0.62 / 1.55 | MBM es el más rápido. |

Problemas del método de medición:
1. phyloseq y rbiom devuelven un ggplot sin dibujar; MicroBioMeta sí dibuja dentro del tiempo medido.
2. Los tiempos de 0.13–0.25 s de rbiom salieron de su caché en disco. Sin caché tarda entre 0.9 y 1.5 s.
3. Una sola repetición y sin calentamiento, así que se mezcla la primera carga de ALDEx2 (5.4 s) y de ComplexHeatmap (1.5 s).
4. Se compararon cálculos distintos: una PERMANOVA composicional con interacción contra pruebas univariadas o de un solo factor.
5. No se midió cuánto tarda construir el objeto phyloseq ni limpiar su taxonomía (unas 11 líneas).
6. `devtools::load_all()` tarda 8.75 s, frente a 0.78 s de `library(MicroBioMeta)`.

## 2. Mejoras de velocidad, por prioridad

| # | Cambio | Dónde | Ganancia medida | Esfuerzo |
|---|---|---|---|---|
| 1 | Calcular solo la muestra Monte Carlo que se usa (misma semilla, resultado `identical()`) | `beta_test_table.R:148`, `beta_ord_plot.R:155` | 12.3 → 6.2 s | Bajo |
| 1b | O bien un CLR determinista (digamma o `log(x+0.5)` centrado); esto cambia el método | igual | 12.3 → 0.46 s | Bajo (decisión científica) |
| 2 | No dibujar el heatmap dos veces; devolver el objeto | `abundance_heatmap_plot.R:384-395` | alrededor de 2× | Bajo |
| 3 | Leer la taxonomía solo en las cadenas únicas o del top_n, no en las 6045 filas; `findInterval` en lugar de `case_when` | `abundance_heatmap_plot.R:110-184` | 0.72 → 0.40 s (idéntico) | Bajo |
| 4 | Un solo ggplot facetado en lugar de 12 paneles con cowplot; no dibujar dentro de la función | `alpha_diversity_plot.R:300-481`, `alpha_hill_plot.R` | 3.3 → 1.3 s (cambia el diseño) | Medio |
| 5 | Crear el tema una vez, no 13 | `zzz.R:152` | alrededor del 16 % de alpha | Bajo |
| 6 | Pasar phyloseq (7–10 s de carga), ALDEx2, ANCOMBC y randomForest a Suggests con `requireNamespace()` | `DESCRIPTION`, `NAMESPACE` | menos tiempo de instalación y de primer uso | Medio |
| 7 | Correlaciones con p-valores vectorizados en lugar de un doble bucle de `cor.test` | `corr_env_abund_plot.R:451-463` | no medido | Bajo |
| 8 | `beta_turnover_plot`: versión matricial (`tcrossprod`) en lugar de 3 llamadas a `hill_taxa_parti_pairwise` | `beta_diversity.R:146-163` | no medido | Medio |
| 9 | Funciones compartidas para leer la taxonomía (el mismo `case_when` está copiado en 6 archivos) | varios | mantenimiento | Medio |

**Decisión metodológica pendiente:** hoy, con una sola muestra Monte Carlo, la PERMANOVA composicional depende de la semilla. Entre 5 semillas, el p de la interacción varió de 0.29 a 0.62 y el R² del suelo de 0.172 a 0.065. Conviene usar un CLR determinista, o promediar las distancias o los resultados de las 128 muestras.

## 3. Errores de resultados (comprobados en el código)

- `beta_test_table.R:97-98,155`: si recibe un objeto `dist`, vuelve a pasarlo por `vegdist` y calcula distancias sobre la matriz de distancias. La PERMANOVA sale mal sin ningún aviso.
- `beta_dissimilarity_plot.R:119`: `distinct(site1, site2)` no elimina los pares espejo, así que (A,B) y (B,A) cuentan doble, lo que infla la n y sesga los p-valores. Con `partition = "shared"` (el valor por defecto) además entra la diagonal (cada muestra comparada consigo misma), y el filtro `value != 0` de la línea 108 elimina pares reales con disimilitud 0.
- Otros reportados en `CODE_REVIEW.md`, aún sin comprobar: metadata sin alinear (`beta_test_table.R:115-120`, `random_forest_lollipop_plot.R:63-70`), `cca_rda_biplot.R:93` con `metadata = NULL`, y posible inversión de las etiquetas de ALDEx2.

## 4. Ventajas frente a la competencia (para vender el paquete)

1. **Líneas de código:** unas 8 llamadas para todo el script, frente a unas 40 en phyloseq.
2. **Ecología espacial y partición de beta:** betapart (turnover/nestedness) y decaimiento con la distancia geográfica con Mantel, cada uno en una función. Ningún competidor revisado lo ofrece así.
3. **Números de Hill (q0/q1/q2) por defecto,** también sobre gradientes continuos.
4. **Enfoque composicional por defecto** (CLR/Aitchison), con ALDEx2 y ANCOM-BC2 ya graficados.
5. **Entrada directa de tablas de QIIME2 y Kraken2,** sin construir objetos. Una llamada entrega cálculo, estadística y figura lista para publicar.
6. **Documentación en español e inglés** y datos de ejemplo reales.

Lema sugerido: *"From QIIME2/Kraken2 tables to publication-ready microbial ecology in one line per figure."*

## 5. Lo que falta (sobre todo para Bioconductor)

- **Probable rechazo:** no acepta phyloseq, SummarizedExperiment/TreeSummarizedExperiment ni BIOM.
- **Prohibido por Bioconductor:** escribe archivos en el directorio de trabajo (`sankey.html`) y llama a `set.seed()` dentro de las funciones.
- **Resultados:** las tablas se guardan en disco en lugar de devolverse junto con el gráfico.
- **Estadística:** no hay árbol ni UniFrac/Faith, y no se corrigen las comparaciones múltiples en `corr_env_abund_plot` ni en las pruebas por pares.
- **Coherencia del paquete:** la posición de la columna de taxonomía varía entre funciones, y `alpha_hill_plot` casi duplica a `alpha_diversity_plot`.

## 6. Cómo reportar el benchmark en el artículo

Por tarea, reportar:
- número de llamadas y líneas de código hasta la figura guardada;
- qué incluye la figura (estadística, agrupación, etiquetas);
- tiempos con `bench::mark` que incluyan guardar la figura, con configuraciones equivalentes, 2–3 tamaños de datos y la caché de rbiom desactivada.

## 7. Devolver objetos ggplot modificables

Varias funciones dibujan la figura por dentro y devuelven una imagen ya armada. Conviene que devuelvan un objeto ggplot sin dibujar, como phyloseq y rbiom. Así la usuaria puede agregarle `+ theme()` o `+ labs()`, o guardarlo con `ggsave()`. También es lo que espera un revisor de Bioconductor, y hace justa la comparación de tiempos.

En la consola no se pierde nada: si la función solo devuelve el objeto, R lo muestra solo. Dentro de un `for` o de otra función hay que usar `print()`.

La velocidad total casi no cambia, porque el gráfico se dibuja igual cuando se muestra. Hay dos excepciones: el heatmap deja de dibujarse dos veces, y en alfa la ganancia (de 3.3 a 1.3 s) viene de armar un solo gráfico con facetas en lugar de 12.

| Función | Qué devuelve hoy | Cómo dejarlo como ggplot |
|---|---|---|
| `alpha_diversity_plot`, `alpha_hill_plot` | Imagen armada con cowplot (12 paneles); no admite `+ theme()` | Un solo ggplot con `facet_grid()`; letras A/B/C como `geom_text` o etiquetas de faceta |
| `alpha_decay_plot`, `alpha_hill_corr_plot`, `beta_partition_ord_plot` | Varios paneles unidos con `cowplot::plot_grid` | Facetas, o patchwork |
| `abundance_heatmap_plot`, `aldex_heatmap_plot` | Heatmap de ComplexHeatmap ya dibujado | Ver abajo: ComplexHeatmap no es ggplot2 |
| `beta_test_table` | Tabla como imagen (`ggtexttable`) | Devolver el data.frame con los resultados y dejar la imagen como opción |
| `abundance_sankey_plot` | HTML de networkD3 que se guarda en disco | No es ggplot; basta con no guardarlo por defecto |

Las demás (barras, ordenación, CCA/RDA, venn, volcano de ALDEx2, etc.) ya devuelven ggplot.

**Figuras con varios paneles distintos: patchwork.** Cuando no se pueden hacer con facetas, patchwork en lugar de cowplot:
- no dibuja nada hasta que se imprime;
- permite `p & theme_bw()` para todos los paneles, o `p[[2]] + labs(...)` para uno solo;
- pone las letras A/B/C con `plot_annotation(tag_levels = "A")`;
- funciona con `ggsave()`.

**Los heatmaps: dos caminos.**
- **Mantener ComplexHeatmap** y devolver el objeto sin dibujar: se dibuja al imprimirlo y conserva anotaciones y agrupamiento, pero no admite `+ theme()`.
- **Pasarlo a ggplot2** con `geom_tile`: consistente con el resto del paquete; las anotaciones se hacen con ggh4x o patchwork, con más trabajo.
