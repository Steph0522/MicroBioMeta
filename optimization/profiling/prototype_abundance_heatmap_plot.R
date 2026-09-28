# prototype_abundance_heatmap_plot.R
# Hotspots (Rprof):
#   R/abundance_heatmap_plot.R#384-392 grid.grabExpr(ComplexHeatmap::draw(...))  ~46%  (layout + dibujo offscreen)
#   R/abundance_heatmap_plot.R#110-142 pipeline de taxonomia sobre TODAS las filas  ~28%
#        (tidyr::separate + case_when con grepl/str_extract sobre 6045 ASVs, antes de quedarse con top_n = 20)
#   R/abundance_heatmap_plot.R#165-184 case_when de 13 ramas via across() en 79 columnas ~8%
#        (+ JIT/compilacion de la lambda en cada llamada: 'cmpfun' ~30% del perfil)
#   R/abundance_heatmap_plot.R#394-395 grid.newpage()+grid.draw(): se DIBUJA OTRA VEZ en el dispositivo activo.
#
# Variantes (patch_fun sobre el codigo real de R/, sin modificarlo):
#   H1: binning con findInterval() vectorizado sobre la matriz (en vez de case_when x 79 columnas).
#   H2: H1 + clasificacion taxonomica sobre taxonomias UNICAS y tidyr::separate solo en las top_n filas.
#   H3: H2 + no dibujar dentro de la funcion (devolver el gTree; el usuario lo dibuja al imprimir).
# Verificacion: tabla guardada (save_table) identica y render PNG pixel-identico del objeto devuelto.
# Uso: cd optimization/profiling; Rscript prototype_abundance_heatmap_plot.R
source("00_setup_data.R")
N <- as.integer(Sys.getenv("NREPS", "6"))

src_lines <- readLines(file.path(PKG_ROOT, "R", "abundance_heatmap_plot.R"))
grab <- function(from_pat, to_pat) {             # bloque de texto entre dos anclas (inclusive)
  i <- grep(from_pat, src_lines, fixed = TRUE)[1]
  j <- i - 1 + grep(to_pat, src_lines[i:length(src_lines)], fixed = TRUE)[1]
  paste(src_lines[i:j], collapse = "\n")
}
# --- H1: bloque de binning (#170-184) --------------------------------------------------------------
old_bin <- grab("dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(.))) %>%", ". >  75.00 ~ 12)))")
new_bin <- paste(
  "    as.matrix()",
  "  heatmap <- matrix(as.numeric(heatmap), nrow(heatmap), dimnames = dimnames(heatmap))",
  "  heatmap[] <- as.numeric(findInterval(heatmap, c(0.001, 0.005, 0.01, 0.10, 0.20, 1, 2, 5, 10, 25, 50, 75), left.open = TRUE))",
  sep = "\n")
# --- H2: pipeline de taxonomia (#110-142) ----------------------------------------------------------
old_tax <- grab("table_abundance <- phy.ra.complete %>%", "dplyr::slice(seq_len(top_n)) %>%")
new_tax <- paste(
  "  table_abundance <- phy.ra.complete %>%",
  "    dplyr::mutate(abun = rowMeans(.)) %>%",
  "    tibble::rownames_to_column(var=\"OTUID\") %>%",
  "    dplyr::left_join(table_tax, by = \"OTUID\")  %>%",
  "    dplyr::mutate(taxonomy2 = taxonomy)",
  "  .u <- unique(table_abundance$taxonomy)",
  "  .uc <- dplyr::case_when(",
  "      grepl(\"g__[^;]*\", .u) & !grepl(\"g__uncultured|g__$\", .u) ~ sub(\".*g__([^;]*).*\", \"\\\\1\", .u),",
  "      grepl(\"f__[^;]*\", .u) & !grepl(\"f__uncultured|f__$\", .u) ~ paste0(\"other \", stringr::str_extract(.u, \"f__[^;]*\") %>% sub(\"f__\", \"\", .)),",
  "      grepl(\"o__[^;]*\", .u) & !grepl(\"o__uncultured|o__$\", .u) ~ paste0(\"other \", stringr::str_extract(.u, \"o__[^;]*\") %>% sub(\"o__\", \"\", .)),",
  "      grepl(\"c__[^;]*\", .u) & !grepl(\"c__uncultured|c__$\", .u) ~ paste0(\"other \", stringr::str_extract(.u, \"c__[^;]*\") %>% sub(\"c__\", \"\", .)),",
  "      grepl(\"p__[^;]*\", .u) & !grepl(\"p__uncultured|p__$\", .u) ~ paste0(\"other \", stringr::str_extract(.u, \"p__[^;]*\") %>% sub(\"p__\", \"\", .)),",
  "      TRUE ~ \"Unclassified\")",
  "  table_abundance$taxonomy <- .uc[match(table_abundance$taxonomy, .u)]",
  "  table_abundance <- table_abundance %>%",
  "    { if (exclude_unclassified) dplyr::filter(., taxonomy != \"Unclassified\") else . } %>%",
  "    dplyr::arrange(-abun) %>%",
  "    dplyr::slice(seq_len(top_n)) %>%",
  "    tidyr::separate(taxonomy2, into = c(\"dominio\",\"phylum\",\"clase\",\"orden\",\"familia\",\"genero\",\"especie\"),",
  "                    sep = \";\", fill = \"right\", extra = \"merge\") %>%",
  sep = "\n")
# --- H3: sin dibujo interno --------------------------------------------------------------------------
old_draw <- "  grid::grid.newpage()\n  grid::grid.draw(heatmap_output)\n  return(invisible(heatmap_output))"
new_draw <- "  return(heatmap_output)"

f_orig <- patch_fun("abundance_heatmap_plot.R", "abundance_heatmap_plot")
f_H1 <- patch_fun("abundance_heatmap_plot.R", "abundance_heatmap_plot", list(c(old_bin, new_bin)))
f_H2 <- patch_fun("abundance_heatmap_plot.R", "abundance_heatmap_plot", list(c(old_bin, new_bin), c(old_tax, new_tax)))
f_H3 <- patch_fun("abundance_heatmap_plot.R", "abundance_heatmap_plot", list(c(old_bin, new_bin), c(old_tax, new_tax),
                                                                              c(old_draw, new_draw)))
args <- list(table = table_bac, metadata = metadata_bacteria, top_n = 20, show_column_names = FALSE,
             condition1 = "Treatment", condition2 = "Type_of_soil", feature_prefix = "ASV")
call_f <- function(f, ...) suppressWarnings(suppressMessages(do.call(f, utils::modifyList(args, list(...)))))

# ---- Verificacion ------------------------------------------------------------------------------
tabs <- lapply(list(pkg = abundance_heatmap_plot, orig = f_orig, H1 = f_H1, H2 = f_H2, H3 = f_H3), function(f) {
  tf <- tempfile(fileext = ".txt"); g <- call_f(f, save_table = TRUE, table_filename = tf)
  list(tab = readLines(tf), px = render_png_pixels(g))
})
for (nm in c("orig", "H1", "H2", "H3"))
  cat(sprintf("%-5s tabla guardada identica: %s | PNG pixel-identico: %s\n", nm,
              identical(tabs$pkg$tab, tabs[[nm]]$tab), identical(tabs$pkg$px, tabs[[nm]]$px)))
# Matriz binned identica (H1 vs original) sobre mas top_n
bin_check <- function(f) render_png_pixels(call_f(f, top_n = 200, cluster = FALSE), width = 1400, height = 4000)
cat("top_n = 200, cluster = FALSE, PNG identico original vs H2:", identical(bin_check(f_orig), bin_check(f_H2)), "\n")

# ---- Benchmark ---------------------------------------------------------------------------------
bm <- function(f, n = N, device = c("null", "png")) {
  device <- match.arg(device); t <- tr <- numeric(n)
  for (i in seq_len(n)) { gc(FALSE)
    if (device == "png") { tf <- tempfile(fileext = ".png"); grDevices::png(tf, 1400, 1000, res = 110) }
    t[i] <- system.time(g <- call_f(f))[["elapsed"]]
    if (device == "png") { grDevices::dev.off(); unlink(tf) }
    tr[i] <- system.time(render_null(g))[["elapsed"]] }
  c(first = t[1], build_med = median(t[-1]), render_med = median(tr[-1]))
}
res <- rbind(
  "original (paquete), pdf(NULL)"       = bm(abundance_heatmap_plot),
  "original (paquete), png activo"      = bm(abundance_heatmap_plot, device = "png"),
  "H1 findInterval"                     = bm(f_H1),
  "H2 H1 + taxonomia unica/top_n"       = bm(f_H2),
  "H3 H2 + sin dibujo interno"          = bm(f_H3),
  "H3, png activo"                      = bm(f_H3, device = "png"),
  "phyloseq plot_heatmap (ref.)"        = bm(function(...) { ps2 <- prune_taxa(taxa_sums(ps_bacteria) > 0, ps_bacteria)
                                            plot_heatmap(tax_glom(ps2, taxrank = "Phylum"), taxa.label = "Phylum",
                                                         sample.label = "Type_of_soil") })
)
print(round(res, 3))
write.csv(data.frame(variant = rownames(res), round(res, 3)), "results_prototype_abundance_heatmap_plot.csv", row.names = FALSE)
