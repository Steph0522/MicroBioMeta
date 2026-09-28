# prototype_alpha_diversity_plot.R
# Hotspots (Rprof, ver results_03_profile_summary.txt):
#   R/alpha_diversity_plot.R#451  cowplot::plot_grid(plotlist = 12 ggplots)   ~47%  (ggplotGrob de cada panel)
#   R/alpha_diversity_plot.R#435  build_cell() x 12 (4 suelos x 3 indices)    ~20%
#   R/alpha_diversity_plot.R#465  2o cowplot::plot_grid(panel_grid, leg)       ~10%
#   R/alpha_diversity_plot.R#444  get_legend(build_cell(1,1)) -> 13a ggplot     ~7%
#   R/zzz.R#152 (.mbm_theme -> theme_bw) llamado 13 veces                      ~16% (incluido arriba)
#   Calculo de diversidad (#145-148, vegan::estimateR/diversity ya vectorizado) ~6%
# Es decir, el coste es RENDERIZAR (convertir a gtable) 13 ggplots dentro de la funcion,
# no el calculo. phyloseq/rbiom devuelven UN ggplot sin renderizar.
#
# Variantes (generadas con patch_fun() sobre el codigo real de R/, sin modificarlo):
#   A: reutiliza plots[[1]] para la leyenda en vez de construir un 13er ggplot (salida identica).
#   C: A + usar matriz en vez de data.frame(t(table)) en #138 (salida identica).
#   B: usa la ruta ya existente de un unico ggplot facetado (facet_grid2 + etiquetas A/B/C via gtable),
#      forzando use_grid_compose = FALSE. Un solo ggplotGrob. Salida visual similar, NO identica.
# Uso: cd optimization/profiling; Rscript prototype_alpha_diversity_plot.R
source("00_setup_data.R")
N <- as.integer(Sys.getenv("NREPS", "4"))

f_orig <- patch_fun("alpha_diversity_plot.R", "alpha_diversity_plot")
f_A <- patch_fun("alpha_diversity_plot.R", "alpha_diversity_plot", list(c(
  "leg <- cowplot::get_legend(build_cell(1, 1) + ggplot2::theme(legend.position = legend_position))",
  "leg <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = legend_position))")))
f_C <- patch_fun("alpha_diversity_plot.R", "alpha_diversity_plot", list(c(
  "leg <- cowplot::get_legend(build_cell(1, 1) + ggplot2::theme(legend.position = legend_position))",
  "leg <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = legend_position))"),
  c("table <- data.frame(t(table))", "table <- t(as.matrix(table))")))
f_B <- patch_fun("alpha_diversity_plot.R", "alpha_diversity_plot", list(c(
  "use_grid_compose <- is.null(facet_by2) &&", "use_grid_compose <- FALSE && is.null(facet_by2) &&")))

args <- list(table = table_bac, metadata = metadata_bacteria, x_col = "Treatment", fill_col = "Treatment",
             facet_by = "Type_of_soil", facet_orientation = "horizontal", save_table = FALSE)
call_f <- function(f) suppressWarnings(suppressMessages(do.call(f, args)))

# ---- Verificacion ------------------------------------------------------------------------------
p_pkg  <- call_f(alpha_diversity_plot)
p_orig <- call_f(f_orig)
p_A    <- call_f(f_A)
p_B    <- call_f(f_B)
p_C    <- call_f(f_C)
px_pkg  <- render_png_pixels(p_pkg,  "alpha_original.png")
px_orig <- render_png_pixels(p_orig)
px_A    <- render_png_pixels(p_A,    "alpha_variantA.png")
px_B    <- render_png_pixels(p_B,    "alpha_variantB_facet_unico.png")
px_C    <- render_png_pixels(p_C)
cat("Pixeles identicos paquete vs copia patch_fun sin cambios:", identical(px_pkg, px_orig), "\n")
cat("Pixeles identicos original vs A (leyenda reutilizada):  ", identical(px_pkg, px_A), "\n")
cat("Pixeles identicos original vs C (A + matriz en #138):  ", identical(px_pkg, px_C), "
")
cat("Fraccion de pixeles distintos original vs B:            ", mean(px_pkg != px_B), "(esperado: distinto layout)\n")

# Los valores de diversidad son los mismos en ambas rutas (mismo codigo #134-181)
div_vals <- function() {
  tb <- table_bac[, intersect(colnames(table_bac), metadata_bacteria$SAMPLEID)]
  tb <- data.frame(t(tb))
  c(vegan::estimateR(tb)["S.chao1", ], vegan::diversity(tb, "shannon"), vegan::diversity(tb, "simpson"))
}

# ---- Benchmark (build = llamada; render = print en pdf(NULL)) -----------------------------------
bm <- function(f, n = N) {
  tb <- tr <- numeric(n)
  for (i in seq_len(n)) { gc(FALSE)
    tb[i] <- system.time(p <- call_f(f))[["elapsed"]]
    tr[i] <- system.time(render_null(p))[["elapsed"]] }
  c(build_first = tb[1], build_med = median(tb[-1]), render_med = median(tr[-1]), total_med = median(tb[-1] + tr[-1]))
}
res <- rbind(
  "original (paquete)"                 = bm(alpha_diversity_plot),
  "A: leyenda de plots[[1]]"           = bm(f_A),
  "C: A + matriz en #138"              = bm(f_C),
  "B: un solo ggplot facetado"         = bm(f_B),
  "solo calculo diversidad (#145-148)" = { t <- replicate(N, system.time(div_vals())[["elapsed"]])
                                           c(t[1], median(t[-1]), NA, median(t[-1])) },
  "phyloseq plot_richness (ref.)"      = bm(function(...) plot_richness(ps_bacteria, x = "Type_of_soil",
                                          color = "Treatment", measures = c("Chao1", "Shannon", "Simpson")) +
                                          facet_grid(Type_of_soil ~ variable))
)
print(round(res, 3))

# ---- Sub-hotspot: #138 data.frame(t(table)) de 79 x 6045 columnas vs matriz ---------------------
div_df  <- function() { tb <- data.frame(t(table_bac[, metadata_bacteria$SAMPLEID]))
  list(vegan::estimateR(tb)["S.chao1", ], vegan::diversity(tb, "shannon"), vegan::diversity(tb, "simpson")) }
div_mat <- function() { tb <- t(as.matrix(table_bac[, metadata_bacteria$SAMPLEID]))
  list(vegan::estimateR(tb)["S.chao1", ], vegan::diversity(tb, "shannon"), vegan::diversity(tb, "simpson")) }
cat("\nDiversidad data.frame vs matriz, all.equal:", isTRUE(all.equal(lapply(div_df(), unname), lapply(div_mat(), unname))), "\n")
sub_res <- rbind(
  "data.frame(t(table)) + vegan (actual)" = median(replicate(5, system.time(div_df())[["elapsed"]])),
  "t(as.matrix(table)) + vegan"           = median(replicate(5, system.time(div_mat())[["elapsed"]])))
print(sub_res)
res <- rbind(res, "sub: diversidad con data.frame (#138)" = c(NA, sub_res[1], NA, sub_res[1]),
                  "sub: diversidad con matriz"            = c(NA, sub_res[2], NA, sub_res[2]))
write.csv(data.frame(variant = rownames(res), round(res, 3)), "results_prototype_alpha_diversity_plot.csv", row.names = FALSE)
