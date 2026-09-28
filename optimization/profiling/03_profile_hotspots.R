# 03_profile_hotspots.R
# Perfilado con Rprof(line.profiling = TRUE) de las funciones lentas.
# options(keep.source = TRUE) ANTES de load_all para que haya srcref (file:line) de R/*.R.
# Uso: cd optimization/profiling; Rscript 03_profile_hotspots.R
options(keep.source = TRUE, keep.source.pkgs = TRUE)
source("00_setup_data.R")

# warm-up (carga perezosa de namespaces: ALDEx2, ComplexHeatmap, ggpubr, cowplot...)
invisible(suppressWarnings(suppressMessages({
  abundance_bar_plot(table = table_bac, metadata = metadata_bacteria, taxonomy_db = "silva", level = "phylum",
                     top_n = 15, x_col = "Type_of_soil", facet_col = "Treatment", add_remained = TRUE,
                     label = "Phylum", x_label_angle = 45)
  abundance_heatmap_plot(table = table_bac, metadata = metadata_bacteria, top_n = 20, show_column_names = FALSE,
                         condition1 = "Treatment", condition2 = "Type_of_soil", feature_prefix = "ASV")
  beta_test_table(table = table_bac, metadata = metadata_bacteria, formula_str = "Treatment",
                  method = "bray", permutations = 9)
})))

targets <- list(
  beta_test_table = function() beta_test_table(table = table_bac, metadata = metadata_bacteria,
      formula_str = "Treatment*Type_of_soil", method = "compositional", test = "permanova", permutations = 999),
  alpha_diversity_plot = function() for (i in 1:3) alpha_diversity_plot(table = table_bac, metadata = metadata_bacteria,
      x_col = "Treatment", fill_col = "Treatment", facet_by = "Type_of_soil",
      facet_orientation = "horizontal", save_table = FALSE),
  abundance_heatmap_plot = function() for (i in 1:5) abundance_heatmap_plot(table = table_bac, metadata = metadata_bacteria,
      top_n = 20, show_column_names = FALSE, condition1 = "Treatment", condition2 = "Type_of_soil",
      feature_prefix = "ASV"),
  abundance_bar_plot = function() for (i in 1:10) abundance_bar_plot(table = table_bac, metadata = metadata_bacteria,
      taxonomy_db = "silva", level = "phylum", top_n = 15, x_col = "Type_of_soil", facet_col = "Treatment",
      add_remained = TRUE, label = "Phylum", x_label_angle = 45, save_table = FALSE)
)

sink_file <- "results_03_profile_summary.txt"
cat("", file = sink_file)
out <- function(...) cat(..., "\n", file = sink_file, append = TRUE)

for (nm in names(targets)) {
  prof <- sprintf("rprof_%s.out", nm)
  Rprof(prof, interval = 0.005, line.profiling = TRUE, memory.profiling = FALSE)
  el <- system.time(suppressWarnings(suppressMessages(targets[[nm]]())))[["elapsed"]]
  Rprof(NULL)
  s_fun  <- summaryRprof(prof)
  s_line <- summaryRprof(prof, lines = "show")
  s_both <- summaryRprof(prof, lines = "both")
  out("==========================================================")
  out("FUNCION:", nm, "| elapsed total:", el, "s")
  out("--- by.total (funciones, top 40) ---")
  capture.output(print(head(s_fun$by.total, 40)), file = sink_file, append = TRUE)
  out("--- by.self (funciones, top 25) ---")
  capture.output(print(head(s_fun$by.self, 25)), file = sink_file, append = TRUE)
  out("--- lineas (by.line, solo R/*.R del paquete, ordenado por total) ---")
  bl <- s_both$by.total
  bl <- bl[grepl("^[A-Za-z_]+\\.R#[0-9]+$", rownames(bl)) &
             sub("#.*", "", rownames(bl)) %in% basename(list.files(file.path(PKG_ROOT, "R"))), , drop = FALSE]
  capture.output(print(head(bl, 30)), file = sink_file, append = TRUE)
  cat(nm, "done", el, "s\n")
}
cat("Resumen escrito en", sink_file, "\n")
