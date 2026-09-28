# Análisis competitivo y de posicionamiento — MicroBioMeta 0.99.0

Fecha: 2026-09-25 · Alcance: paquete local `C:\Users\HP\Documents\MicroBioMeta` (sin modificar ningún archivo del paquete) + investigación web de competidores (estado 2025-2026).

---

## 0. Resumen ejecutivo

- **MicroBioMeta no gana en velocidad y no debería venderse por velocidad.** Gana en *amplitud ecológica por llamada*: cada función hace cálculo + estadística + figura lista para publicar + (opcional) tabla de resultados, a partir de **tablas planas** (TSV/`data.frame`) tal como salen de QIIME2/Kraken2, sin construir objetos.
- **Nicho diferencial real:** ecología microbiana *ambiental/de suelo* para principiantes: números de Hill (q0/q1/q2) por defecto, partición de beta-diversidad (recambio/anidamiento, `betapart`), decaimiento con la distancia geográfica (Mantel + haversine), diversidad alfa a lo largo de gradientes, CCA/RDA con selección de variables por `envfit`, y distancia composicional (CLR/Aitchison) como opción por defecto. Ningún competidor revisado ofrece *todo eso en una llamada cada uno*.
- **Brechas principales:** (1) no acepta `phyloseq`, `SummarizedExperiment`/`TreeSummarizedExperiment` (bloqueante para Bioconductor), (2) sin árbol filogenético/UniFrac/Faith, (3) contrato de entrada inconsistente (columna `taxonomy` primera vs última, `SAMPLEID` fijo), (4) funciones devuelven sólo la figura (la tabla va a disco), (5) sin corrección por comparaciones múltiples en correlaciones/pruebas en gráficos, (6) ~37 dependencias en Imports.
- **El benchmark del colaborador no es justo** (salidas no equivalentes, render asimétrico, preprocesamiento excluido, una sola repetición). La métrica a reportar es **"líneas de código / llamadas hasta la figura publicable + análisis incluidos + tiempo con render"**.

---

## 1. Inventario de MicroBioMeta

### 1.1 Funciones exportadas (23 + `%>%`)

| Área | Función | Qué hace (cálculo + estadística + salida) |
|---|---|---|
| Preprocesado | `merge_feature_taxonomy()` | Une tabla de conteos + taxonomía (formato QIIME2). |
| Preprocesado | `collapse_table()` | Colapsa a un nivel taxonómico; opción abundancia relativa; devuelve `list(collapsed_table, long_format)`. |
| Composición | `abundance_bar_plot()` | Barras apiladas de abundancia relativa media por grupo; top-n + "Other"; limpieza de prefijos para SILVA/GG2/UNITE/Kraken2; facetas. |
| Composición | `abundance_heatmap_plot()` | Heatmap (ComplexHeatmap) de top-n ASV/taxa con hasta 3 anotaciones de condición + anotación de filo. |
| Composición | `abundance_sankey_plot()` | Sankey interactivo (networkD3, HTML) del flujo taxonómico D→S. |
| Composición | `venn_plot()` | Venn de taxa compartidos (ggvenn o ggVennDiagram), filtro de prevalencia. |
| Composición | `ratios_bubble_plot()` | Razones de abundancia A/B por taxón en burbujas. |
| Alfa | `alpha_diversity_plot()` / `alpha_hill_plot()` | Números de Hill q0/q1/q2 (hillR), box/bar plot, facetado doble, prueba (`stat_compare_means`), etiquetas de panel A/B/C, rarefacción opcional (sólo `alpha_diversity_plot`). *Casi duplicadas.* |
| Alfa | `alpha_hill_corr_plot()` | Correlación entre q0, q1, q2. |
| Alfa | `alpha_decay_plot()` | Hill q0–q2 vs variable continua (gradiente), Spearman/Pearson + R²/pendiente, por grupo. |
| Beta | `beta_ord_plot()` | PCA/PCoA/NMDS; distancias Bray, Jaccard, Sørensen, euclídea, **composicional (CLR ALDEx2, por defecto)**, Aitchison, robust-Aitchison; biplot de taxa en PCA; estrés NMDS; escala continua/discreta automática. |
| Beta | `beta_test_table()` | PERMANOVA (`adonis2`, `by="terms"`, fórmula con interacción, `strata`) o `betadisper`; tabla formateada. |
| Beta | `beta_partition_ord_plot()` | Partición `betapart` (Jaccard/Sørensen): total, recambio, anidamiento; tres ordenaciones con `betadisper`. |
| Beta | `beta_dissimilarity_plot()` | Boxplots de disimilitud compartida/recambio/anidamiento entre grupos + prueba. |
| Beta | `beta_turnover_plot()` | Recambio intra- vs inter-grupo por pares. |
| Beta | `beta_decay_plot()` | **Distance-decay**: disimilitud vs distancia geográfica (km, `geosphere`), prueba de Mantel, R², pendiente, por grupo. |
| Abundancia diferencial | `aldex_volcano_plot()` / `aldex_heatmap_plot()` | ALDEx2 (CLR, BH) → volcano/effect plot o heatmap anotado. |
| Abundancia diferencial | `ancombc_plot()` | ANCOM-BC2 (efectos fijos y aleatorios, ajuste p, nivel taxonómico) → barras LFC o heatmap. |
| Aprendizaje automático | `random_forest_lollipop_plot()` | Random forest → importancia de variables (lollipop). |
| Ambiente | `corr_env_abund_plot()` | Correlación variables ambientales × taxa (tile/círculo, clustering, umbral p). |
| Ambiente | `cca_rda_biplot()` | CCA/RDA (Hellinger por defecto) + `envfit` → sólo vectores significativos, grupos, flechas escaladas. |

### 1.2 Entradas aceptadas
- `data.frame` con taxa en filas y muestras en columnas + columna `taxonomy` (cadena completa `d__;p__;...`) — **posición inconsistente** según la función (la doc de `abundance_bar_plot` dice "primera columna", `beta_ord_plot`/`ancombc_plot` dicen "última").
- Metadatos `data.frame` con `SAMPLEID` (o primera columna = IDs).
- Tabla ambiental separada (`env_data`) para CCA/RDA y correlaciones.
- Formatos de taxonomía: SILVA, Greengenes2, UNITE, Kraken2 (metagenómica).
- `beta_test_table()` acepta además un objeto `dist`.
- **No acepta** `phyloseq`, `SummarizedExperiment`, `TreeSummarizedExperiment`, BIOM, `.qza` (la viñeta remite a `qiime2R`). `phyloseq` se importa sólo internamente para ANCOM-BC2.

### 1.3 Estadística integrada en figuras
Wilcoxon/Kruskal/ANOVA en alfa y disimilitud (ggpubr), Mantel en decay, Spearman/Pearson + R² en gradientes, `envfit` en CCA/RDA, estrés NMDS, BH en ALDEx2, `p_adj_method` en ANCOM-BC2. **Ausente:** corrección por comparaciones múltiples en `corr_env_abund_plot()` y en comparaciones por pares de ggpubr (no aparece `p.adjust` en el código).

### 1.4 Publicación y público
- Paleta daltónica (Okabe-Ito) por defecto, cursiva en géneros/especies, etiquetas A/B/C de panel, tema homogéneo, ángulo de ejes, facetas.
- `save_table`/`table_filename` homogéneo en todas las funciones → tabla de resultados reproducible.
- Viñetas en **inglés y español** (metabarcoding y metagenómica), datos reales de suelo (16S, ITS, Kraken2) incluidos.
- Público declarado: principiantes con poco R ("minimal coding"). Muchos parámetros por función (hasta 30), pero con valores por defecto sensatos.

---

## 2. Competidores (estado 2025-2026)

| Paquete | Estado / versión | Entrada | Filosofía |
|---|---|---|---|
| **phyloseq** | Bioc 3.23, v1.56.0; en ene-2026 Bioconductor advirtió deprecación por fallos de build en release y devel (issue #1796) | objeto `phyloseq` (importadores BIOM/QIIME/mothur) | Contenedor + gráficos básicos; estadística fuera (vegan, DESeq2) |
| **rbiom** | CRAN 3.1.0 (may-2026) | BIOM, `phyloseq` (`as_rbiom`), matrices | Rápido (C, UniFrac multihilo); gráficos con `stat.by` integrado |
| **microeco** (+file2meco, mecoturn) | CRAN 2.4.0 (sep-2026); microeco 2 en *iMeta* 2026 | `microtable` (R6); file2meco: QIIME2, phyloseq↔, TSE↔, HUMAnN, MetaPhlAn | Clases modulares `trans_*`; muy amplio (redes, null models, 15+ métodos DA) |
| **MicrobiotaProcess** | Bioc 3.23, v1.24.0 | `MPSE` (hereda de TSE), phyloseq, TSE | Gramática tidy (`mp_*`) |
| **microViz** | r-universe/GitHub (manual mar-2026) | `phyloseq` | Visualización + PERMANOVA paralela, Shiny `ord_explore` |
| **amplicon / EasyAmplicon** | GitHub; EasyAmplicon 2 (*Adv. Sci.* 2026) | tablas planas (otutab, metadata, taxonomy) | Pipeline + funciones R de una llamada con estadística (el más parecido en estilo) |
| **mia / miaViz** | Bioc 3.23, v1.20.0; libro OMA + preprint bioRxiv oct-2025 | `TreeSummarizedExperiment` | Ecosistema oficial Bioconductor; escalable; interoperable |
| **microbiomeMarker** | **Eliminado de Bioc 3.20** (mantenedor sin respuesta); GitHub | `phyloseq` | Wrapper de DA (LEfSe, ALDEx2, ANCOM, ANCOMBC, DESeq2, edgeR…) |
| **MicrobiomeStat** | CRAN (may-2026) | objeto lista propio; convertidores phyloseq, QIIME2, DADA2, BIOM, SE/MAE | LinDA/linda2, `generate_*` (plot+stat), **reportes automáticos**, longitudinal |
| **animalcules** | Bioc (R ≥ 4.3) | `MultiAssayExperiment` | Shiny interactivo, biomarcadores |
| **MicrobiomeAnalyst 2.0** | Web (NAR 2023) | tablas/BIOM/mothur subidas | Sin código; no reproducible por script |
| **QIIME 2** (amplicon distribution 2025.x) | Activo; q2-boots (2025) | artefactos `.qza` | Emperor interactivo; `core-metrics`; figuras no "de revista" |

---

## 3. Matriz de funcionalidades

Leyenda: ✅ = sí, en una llamada · 🟡 = posible pero en varios pasos / parcial / vía otro paquete · ❌ = no · ? = no verificado

| Funcionalidad | **MicroBioMeta** | phyloseq | rbiom | microeco | MicrobiotaProcess | microViz | amplicon | mia/miaViz | MicrobiomeStat | microbiomeMarker | MicrobiomeAnalyst |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Entrada: tablas planas (TSV/data.frame) | ✅ | 🟡 | 🟡 | 🟡 | 🟡 | ❌ | ✅ | 🟡 | 🟡 | ❌ | ✅ (web) |
| Entrada: phyloseq | ❌ | ✅ | ✅ | ✅ (file2meco) | ✅ | ✅ | ❌ | ✅ (convert) | ✅ | ✅ | ❌ |
| Entrada: SE/TSE (Bioconductor) | ❌ | ❌ | ❌ | ✅ (file2meco) | ✅ | ❌ | ❌ | ✅ | ✅ | ❌ | ❌ |
| Árbol / UniFrac / Faith | ❌ | ✅ | ✅ (rápido) | ✅ | ✅ | ✅ | 🟡 | ✅ | 🟡 | 🟡 | ✅ |
| Barras abundancia relativa por grupo (top-n + Other) | ✅ | 🟡 (~13 pasos) | 🟡 | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ |
| Heatmap con anotaciones | ✅ | 🟡 | ✅ | ✅ | 🟡 | ✅ | 🟡 | 🟡 | ✅ | 🟡 | ✅ |
| Sankey taxonómico | ✅ | ❌ | ❌ | 🟡 | ❌ | ❌ | 🟡 | ❌ | ❌ | ❌ | ❌ |
| Venn/UpSet | ✅ (Venn) | ❌ | ❌ | ✅ | ✅ | ❌ | ✅ | ❌ | ❌ | ❌ | 🟡 |
| **Números de Hill q0/q1/q2 por defecto** | ✅ | ❌ | ❌ | 🟡 | ❌ | ❌ | ❌ | ❌ (índices clásicos) | ❌ | ❌ | ❌ |
| Alfa + prueba estadística en figura | ✅ | ❌ | ✅ | ✅ | ✅ | 🟡 | ✅ | 🟡 (ggpubr manual) | ✅ | ❌ | ✅ |
| **Alfa a lo largo de gradiente continuo (+R², rho)** | ✅ | ❌ | 🟡 (`adiv_corrplot`) | 🟡 | ❌ | ❌ | ❌ | ❌ | 🟡 | ❌ | ❌ |
| Ordenación PCA/PCoA/NMDS | ✅ | 🟡 (2 pasos) | ✅ (+UMAP/t-SNE) | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ |
| Distancia composicional CLR/Aitchison | ✅ (por defecto) | ❌ | ✅ | ✅ | ✅ | ✅ | ❌ | ✅ | ✅ | 🟡 | ✅ |
| PERMANOVA multifactorial (interacción, strata) | ✅ | ❌ | 🟡 | ✅ | ✅ | ✅ (paralelo) | 🟡 | 🟡 | ✅ | ❌ | 🟡 |
| betadisper (dispersión) | ✅ | ❌ | ❌ | ✅ | ✅ | ❌ | ❌ | 🟡 | ? | ❌ | ❌ |
| **Partición betapart (recambio/anidamiento)** | ✅ (3 funciones) | ❌ | ❌ | 🟡 (mecoturn / manual) | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| **Distance-decay geográfico + Mantel** | ✅ | ❌ | ❌ | 🟡 (`plot_scatterfit` + matriz manual) | 🟡 (`mp_mantel`) | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ |
| CCA/RDA con variables ambientales | ✅ (+envfit) | 🟡 (sin flechas) | ❌ | ✅ | ✅ | ✅ | ✅ (CPCoA) | ✅ (`plotRDA`) | ❌ | ❌ | ❌ |
| Correlación ambiente × taxa | ✅ | ❌ | 🟡 | ✅ | 🟡 | ✅ | 🟡 | 🟡 | 🟡 | ❌ | ✅ |
| ALDEx2 | ✅ (volcano + heatmap) | ❌ | ❌ | ✅ | ❌ | ❌ | ❌ | 🟡 (OMA) | ❌ | ✅ | ❌ |
| ANCOM-BC2 | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ | ❌ | 🟡 (OMA) | ❌ | ✅ (ANCOMBC 1) | ✅ |
| Otros DA (LEfSe, DESeq2, LinDA, MaAsLin) | ❌ | 🟡 | 🟡 (Wilcoxon/KW) | ✅ | ✅ (LEfSe-like) | 🟡 | ✅ | 🟡 | ✅ (LinDA) | ✅ | ✅ |
| Random forest (importancia) | ✅ | ❌ | ❌ | ✅ | ❌ | ❌ | 🟡 | ❌ | ❌ | ✅ | ✅ |
| Tabla de resultados exportable homogénea | ✅ (a disco) | ❌ | 🟡 (`$stats`) | ✅ (en objeto) | ✅ (tidy) | 🟡 | 🟡 | ✅ (colData) | ✅ | ✅ | ✅ |
| Reporte automático | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ | ❌ | ✅ |
| Longitudinal / efectos aleatorios | 🟡 (ANCOM-BC2 `rand_formula`, `strata`) | ❌ | ✅ (emmeans) | 🟡 | ❌ | ❌ | ❌ | 🟡 | ✅ | ❌ | 🟡 |
| Redes / null models | ❌ | 🟡 | ❌ | ✅ | ❌ | ❌ | 🟡 | 🟡 | ❌ | ❌ | ✅ |
| Paralelización / C | ❌ | ❌ | ✅ | ❌ | ❌ | ✅ (PERMANOVA) | ❌ | ✅ (DelayedArray/BiocParallel) | ❌ | ❌ | n/a |
| Docs en español | ✅ | ❌ | ❌ | 🟡 (chino/inglés) | ❌ | ❌ | 🟡 (chino) | ❌ | ❌ | ❌ | ❌ |
| Metagenómica (Kraken2) tratada explícitamente | ✅ | ❌ | ❌ | ✅ (file2meco) | 🟡 | ❌ | ❌ | 🟡 | 🟡 | ❌ | 🟡 |

> Nota: las celdas de competidores se basan en documentación pública (ver fuentes); las marcadas 🟡/? conviene verificarlas antes de citarlas en el artículo.

---

## 4. Ventajas, brechas y riesgos Bioconductor

### 4.a Ventajas únicas o raras
1. **Módulo de ecología espacial/β completo en un solo paquete:** `betapart` (3 funciones: ordenación de componentes, boxplots por grupo, recambio intra/inter) + **distance-decay con coordenadas lat/lon** (haversine + Mantel + R² + por grupo). Ningún competidor lo ofrece como función de una llamada.
2. **Hill numbers como estándar** (q0, q1, q2) en boxplots, correlación entre órdenes y **gradientes continuos** (`alpha_decay_plot`). Alineado con la recomendación actual (Chao, Jost) y ausente por defecto en phyloseq, rbiom, mia, MicrobiotaProcess, MicrobiomeStat.
3. **Composicional por defecto** (CLR de ALDEx2 en ordenación y PERMANOVA) + ALDEx2 y ANCOM-BC2 con figura lista.
4. **Entrada en tablas planas QIIME2/Kraken2** con limpieza de prefijos para 4 bases de datos: cero construcción de objetos — la barrera nº1 del principiante.
5. **CCA/RDA con `envfit`** que dibuja sólo los vectores significativos (phyloseq necesita ~5 pasos y no filtra).
6. **Consistencia de publicación:** paleta daltónica, cursivas, paneles A/B/C, `save_table` uniforme.
7. **Bilingüe (ES/EN)** — mercado hispanohablante (Latinoamérica/España) desatendido.
8. Frente a microbiomeMarker (fuera de Bioc) y phyloseq (en riesgo de deprecación), hay hueco para un paquete *amigable y mantenido*.

### 4.b Brechas
1. **Interoperabilidad:** no acepta `phyloseq`/`TSE`/BIOM; no devuelve objetos reutilizables.
2. **Sin filogenia**: ni UniFrac ni Faith PD ni árbol.
3. **Salida sólo gráfica**: las tablas se escriben a disco (`save_table`) en lugar de devolverse (p. ej. atributo o `list(plot, table)`); PERMANOVA es la excepción.
4. **Contrato de entrada inconsistente** (posición de `taxonomy`, `SAMPLEID` vs "primera columna").
5. **Estadística:** sin corrección de comparaciones múltiples en `corr_env_abund_plot()` ni en pruebas por pares de alfa/beta; la distancia "compositional" usa **una sola** réplica Monte Carlo de ALDEx2 (se calculan 128 y se usa 1: coste ×128 sin beneficio y resultado dependiente de la semilla).
6. **Rendimiento:** sin paralelización de permutaciones; renderizado interno de grobs (cowplot/`ggplotGrob`/`grid.grabExpr`) en cada llamada.
7. **Faltan:** métodos DA adicionales (LinDA, MaAsLin2, DESeq2), redes, longitudinal, UpSet, rarefacción/curvas de acumulación (iNEXT sería natural con Hill), reporte automático.
8. **Redundancia:** `alpha_hill_plot` ≈ `alpha_diversity_plot` (difieren sólo en `rarefy_depth`); `beta_turnover_plot` vive en `beta_diversity.R`.

### 4.c Riesgos para aceptación en Bioconductor
| Riesgo | Evidencia | Gravedad |
|---|---|---|
| No interoperabilidad con clases Bioc (SE/TSE) | Guía Bioc: "package submissions are generally not accepted unless they demonstrate such interoperability" | **Alta (probable rechazo)** |
| Escritura en el directorio de trabajo | `abundance_sankey_plot(output_file = "sankey.html")` por defecto; `table_filename` relativo | Alta (regla explícita Bioc) |
| `set.seed()` dentro de funciones | `beta_ord_plot`, `beta_test_table`, `cca_rda_biplot` (aunque restauran el RNG) | Media (regla explícita) |
| `suppressWarnings()` ×5 | BiocCheck ya lo reporta | Baja |
| Dependencias pesadas en Imports (~37: ANCOMBC, ALDEx2, ComplexHeatmap, randomForest, phyloseq, geosphere…) | DESCRIPTION | Media (mover a Suggests + `requireNamespace`) |
| Funciones muy largas / `sapply` / bucles | 8 222 líneas en R/, funciones de >600 líneas | Baja-media |
| Mensajes/argumentos en español dentro del código (`decimales`, "No hay variables ambientales…") | `beta_test_table`, `cca_rda_biplot` | Baja |
| Duplicación de funciones | `alpha_hill_plot`/`alpha_diversity_plot` | Baja |
| biocView `LongRead` sin soporte específico; instalación sólo vía GitHub en README | DESCRIPTION/README | Baja |

---

## 5. Crítica del benchmark del colaborador

Archivo: `optimization/original_comparison_table.txt` + `original_comparison_script.R`.

1. **Salidas no equivalentes (el problema mayor).**
   - *PERMANOVA*: MicroBioMeta (14.95 s) = CLR ALDEx2 con 128 réplicas Monte Carlo + PERMANOVA de **dos factores con interacción**, 999 permutaciones, `by="terms"`. rbiom `bdiv_stats` (0.02 s) es **univariado** (comparación de distancias por pares), y `distmat_stats` (0.19 s) es **Bray-Curtis con un solo factor**. Se compara una ecuación distinta sobre una distancia distinta. Casi todo el tiempo de MicroBioMeta es `aldex.clr(mc.samples = 128)` del que sólo se usa 1 réplica.
   - *Heatmap*: MicroBioMeta = top-20 ASV con 2 anotaciones + filo; phyloseq = filos sin anotaciones.
   - *Barras*: rbiom muestra top-10 sin promediar réplicas; MicroBioMeta top-15 promediado por grupo con facetas.
   - *Alfa*: MicroBioMeta = 3 números de Hill en paneles compuestos con etiquetas; phyloseq = índices clásicos sin estadística.
2. **Render asimétrico.** Las funciones de MicroBioMeta **dibujan dentro de la llamada** (`ggplotGrob` + `cowplot::plot_grid` en alfa; `grid.grabExpr(ComplexHeatmap::draw())` + `grid.draw` en el heatmap), mientras `plot_richness`, `taxa_stacked`, `adiv_boxplot` devuelven un ggplot **perezoso** que aún no se ha renderizado: el `system.time` no incluye su coste de dibujo. Sólo el bloque de ordenación de phyloseq incluye `print()`. Explica gran parte de 5.42 s vs 0.2 s (alfa) y 4.75 s vs 0.2 s (heatmap).
3. **Preprocesamiento excluido.** No se cronometra la construcción del objeto phyloseq (≈5 sentencias) ni la limpieza de taxonomía (≈6-7 sentencias) ni `as_rbiom()`, que MicroBioMeta no necesita.
4. **Una sola medición, sin calentamiento.** `system.time()` único; la primera llamada paga carga de namespaces (ALDEx2, ComplexHeatmap, vegan). Usar `bench::mark()` con ≥10 iteraciones, mediana e IQR, y sesión fresca.
5. **Un solo conjunto de datos pequeño**: no hay curva de escalado (n muestras × n ASV), que es lo que importa a un revisor.
6. **Errores de redacción en la tabla**: la observación de `beta_ord_plot` dice "Performs and visualizes CCA/RDA" (copiada de `cca_rda_biplot`); texto mezclado ES/EN.
7. **Relevancia práctica**: diferencias de 0.2–5 s son irrelevantes para el usuario principiante frente a los minutos/horas de escribir y depurar 13 líneas de `tax_glom`/`merge_samples`.

### 5.1 Recuento de líneas de código (del propio script del colaborador)

Sentencias necesarias para una salida comparable (excluyendo los 3 `read.delim` comunes). Aproximado.

| Tarea | MicroBioMeta | phyloseq (+vegan) | rbiom | Comentario |
|---|---|---|---|---|
| Preparar objeto | 1 (`merge_feature_taxonomy`) | ~11 (4 constructores + ~7 limpieza de taxonomía) | 1 (`as_rbiom`, desde phyloseq) | |
| Barras (top-15 + Other, media por grupo, facetas) | **1** | 13 | 2 (sin media, top-10) | |
| Heatmap con anotaciones | **1** | 3 (sin anotaciones) | 1 | rbiom empata |
| Alfa + prueba | **1** (Hill) | 2 (sin prueba) | 1 | rbiom empata |
| Ordenación PCoA | **1** | 2-3 | 1 (+adonis2) | rbiom empata |
| PERMANOVA 2 factores + interacción | **1** | 3 (vegan) | 3 (1 factor) | |
| CCA + flechas ambientales | **1** (+1 tabla ambiental) | 5 (sin `envfit`) | n/a | |
| **Total** | **≈8** | **≈40** | **≈9 (sin CCA)** | |

Lectura honesta: frente a **phyloseq**, MicroBioMeta reduce ~5× el código; frente a **rbiom** empata en concisión y pierde en velocidad → la diferencia con rbiom está en *qué análisis existen* (Hill, betapart, decay, CCA/RDA, ALDEx2, ANCOM-BC2, composicional), no en líneas.

### 5.2 Benchmark recomendado para el artículo
Tres ejes por tarea, con salida **equivalente** definida a priori:
1. **Esfuerzo**: nº de llamadas a función, nº de sentencias y nº de paquetes a cargar hasta el PNG final (`ggsave` incluido para todos).
2. **Cobertura**: checklist de lo que incluye la figura (estadística, ajuste p, etiquetas de panel, paleta daltónica, tabla exportada).
3. **Tiempo** con `bench::mark` (≥10 iteraciones, mediana), **incluyendo render** (`ggsave` a archivo temporal), en 2–3 tamaños de datos (p. ej. 50/200/1000 muestras) y con parámetros igualados (misma distancia, mismo nº de permutaciones, mismo diseño).
Opcional: prueba con usuarios (estudiantes) midiendo tiempo hasta la figura y errores — muy convincente para un paquete "for beginners".

---

## 6. Posicionamiento y marketing

### 6.1 Taglines (propuestas)
- **"From QIIME2/Kraken2 tables to publication-ready microbial ecology in one line per figure."**
- "Hill numbers, beta partitioning and distance-decay — no objects, no boilerplate."
- ES: "De la tabla de QIIME2 a la figura de la revista, una línea por análisis."

### 6.2 Qué destacar en el artículo
1. **Ecología de comunidades de primera clase**: Hill numbers, betapart, distance-decay, gradientes, CCA/RDA con envfit — con una figura-resumen tipo "workflow" (ver viñeta de suelo).
2. **Composicionalidad por defecto** (CLR/Aitchison, ALDEx2, ANCOM-BC2): buenas prácticas sin que el usuario las conozca.
3. **Una llamada = cálculo + estadística + figura + tabla**: tabla de comparación de sentencias (5.1) frente a phyloseq/mia.
4. **Accesibilidad**: tablas planas, docs bilingües, paleta daltónica, defaults de revista.
5. **Casos reales**: 16S e ITS de suelo, metagenómica Kraken2.
6. Posicionarlo como **complemento** (capa de "última milla" sobre vegan/betapart/hillR/ALDEx2/ANCOMBC), no como sustituto de mia/phyloseq. Citar honestamente que rbiom es más rápido y mia más escalable.

### 6.3 Hoja de ruta (prioridad)
**P0 — necesario para Bioconductor**
1. Aceptar `TreeSummarizedExperiment`/`SummarizedExperiment` (y `phyloseq`) en todas las funciones mediante un conversor interno único (`.mbm_as_tables(x, assay.type, rank)`), manteniendo tablas planas como vía principiante. Añadir `mbm_from_qiime2()` / `mbm_to_tse()`.
2. No escribir en el directorio de trabajo: `output_file`/`table_filename` sin valor por defecto en wd (usar `tempfile()` o exigir ruta).
3. Devolver resultados: `list(plot = , table = )` o atributo `attr(p, "table")`; `save_table` queda como atajo.
4. Unificar contrato de entrada (columna `taxonomy` por nombre, no por posición; `sample_col` configurable).
5. Mover dependencias pesadas a Suggests; eliminar `set.seed()` interno (usar `withr::with_seed` o dejar la semilla al usuario).

**P1 — rigor y rendimiento**
6. `p_adjust` (BH por defecto) en `corr_env_abund_plot()` y comparaciones por pares.
7. Distancia composicional: `mc.samples` configurable (1 si sólo se usa 1 réplica, o promedio de la mediana CLR) — elimina ~128× coste en `beta_test_table`/`beta_ord_plot`.
8. Permutaciones en paralelo (`parallel` en `adonis2`/`mantel`, o `BiocParallel`) y opción `render = FALSE` (devolver objetos sin componer grobs).
9. Fusionar `alpha_hill_plot` en `alpha_diversity_plot`.

**P2 — competitividad**
10. Árbol filogenético opcional: UniFrac (vía `rbiom`/`phyloseq`), Faith PD, Hill filogenéticos (`hillR::hill_phylo`) → otro diferencial raro.
11. Curvas de rarefacción/extrapolación de Hill (iNEXT) y cobertura.
12. DA adicionales (LinDA, MaAsLin2) + figura de consenso entre métodos.
13. UpSet para >4 grupos; función `mbm_report()` (Rmd/Quarto automático) como MicrobiomeStat.
14. Shiny ligero opcional para principiantes absolutos.

---

## 7. Fuentes

- phyloseq (Bioconductor): https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html
- phyloseq issue #1796 (riesgo de deprecación, ene-2026): https://github.com/joey711/phyloseq/issues/1796
- rbiom (CRAN): https://cran.r-project.org/package=rbiom · manual: https://cran.r-project.org/web/packages/rbiom/refman/rbiom.html · sitio: https://cmmr.github.io/rbiom/
- microeco tutorial: https://chiliubio.github.io/microeco_tutorial/ · trans_diff: https://chiliubio.github.io/microeco_tutorial/model-based-class.html · trans_env: https://rdrr.io/cran/microeco/man/trans_env.html · trans_beta: https://rdrr.io/cran/microeco/man/trans_beta.html
- microeco 2 (iMeta 2026, benchmark de pasos/tiempo): https://pmc.ncbi.nlm.nih.gov/articles/PMC13377413/
- file2meco (phyloseq/TSE): https://cran.r-project.org/web/packages/file2meco/refman/file2meco.html
- MicrobiotaProcess (Bioc): https://bioconductor.org/packages/release/bioc/html/MicrobiotaProcess.html · artículo: https://pmc.ncbi.nlm.nih.gov/articles/PMC9988672/
- microViz: https://david-barnett.github.io/microViz/ · dist_permanova (paralelo): https://david-barnett.github.io/microViz/reference/dist_permanova.html
- amplicon: https://github.com/microbiota/amplicon · EasyAmplicon (iMeta 2023): https://onlinelibrary.wiley.com/doi/full/10.1002/imt2.83 · EasyAmplicon 2: https://pmc.ncbi.nlm.nih.gov/articles/PMC12767004/
- mia (Bioc): https://bioconductor.org/packages/release/bioc/html/mia.html · miaViz: https://www.bioconductor.org/packages/release/bioc/html/miaViz.html
- OMA, diversidad alfa: https://microbiome.github.io/OMA/docs/devel/pages/alpha_diversity.html · beta: https://microbiome.github.io/OMA/docs/devel/pages/beta_diversity.html · preprint: https://www.biorxiv.org/content/10.1101/2025.10.29.685036v1
- microbiomeMarker eliminado en Bioc 3.20: https://www.bioconductor.org/packages/microbiomeMarker · https://support.bioconductor.org/p/9160253/
- MicrobiomeStat: https://cran.r-project.org/package=MicrobiomeStat · https://cafferychen777.github.io/MicrobiomeStat/
- animalcules: https://bioconductor.org/packages/release/bioc/html/animalcules.html · https://link.springer.com/article/10.1186/s40168-021-01013-0
- MicrobiomeAnalyst 2.0: https://academic.oup.com/nar/article/51/W1/W310/7160190
- QIIME 2 amplicon docs (diversity): https://amplicon-docs.qiime2.org/en/2025.4/references/plugins/diversity.html · q2-boots: https://f1000research.com/articles/14-87/v1
- Guías Bioconductor — reutilización de clases: https://contributions.bioconductor.org/reusebioc.html · código R (archivos, set.seed): https://contributions.bioconductor.org/r-code.html
- betapart (CRAN): https://mirror.linux.duke.edu/pub/cran/web/packages/betapart/index.html
