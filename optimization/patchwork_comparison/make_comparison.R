# Genera, para cada funcion pasada a patchwork, las figuras de sus ejemplos
# de la ayuda con la version vieja (cowplot, *_antes.png) y la nueva
# (patchwork, *_despues.png). Correr desde la raiz del paquete:
#   Rscript optimization/patchwork_comparison/make_comparison.R

suppressMessages(devtools::load_all(".", quiet = TRUE))
out <- "optimization/patchwork_comparison"
fns <- c("alpha_hill_plot", "alpha_diversity_plot", "alpha_hill_corr_plot",
         "alpha_decay_plot", "beta_partition_ord_plot")

for (f in fns) sys.source(file.path("optimization/old_versions", paste0(f, "_old.R")),
                          envir = globalenv())

save_plot <- function(p, file) {
  ggplot2::ggsave(file.path(out, file), p, width = 11, height = 7, dpi = 90)
}

# Runs the example code of `fn`'s help page, saving every plot it returns.
run_examples <- function(fn, version) {
  rd <- tools::Rd_db(dir = ".")[[paste0(fn, ".Rd")]]
  code_file <- tempfile(fileext = ".R")
  tools::Rd2ex(rd, code_file)
  exprs <- parse(code_file)
  env <- new.env(parent = globalenv())
  if (version == "antes") assign(fn, get(paste0(fn, "_old")), envir = env)
  i <- 0
  for (e in exprs) {
    val <- withVisible(suppressMessages(suppressWarnings(eval(e, env))))
    if (val$visible && inherits(val$value, "gg")) {
      i <- i + 1
      save_plot(val$value, sprintf("%s_ej%d_%s.png", fn, i, version))
    }
  }
  i
}

for (f in fns) for (v in c("antes", "despues")) {
  n <- run_examples(f, v)
  cat(sprintf("%-24s %-8s %d figuras\n", f, v, n))
}

# Modos que no salen en los ejemplos: facet_by2 y facet_orientation = "vertical"
tab <- read.delim(system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta"),
                  row.names = 1, check.names = FALSE)
meta <- read.delim(system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta"),
                   check.names = FALSE)
extra <- list(
  facet_by2 = list(x_col = "Treatment", fill_col = "Treatment",
                   facet_by = "Location", facet_by2 = "Month"),
  vertical  = list(x_col = "Treatment", fill_col = "Treatment",
                   facet_by = "Location", facet_orientation = "vertical")
)
for (f in c("alpha_hill_plot", "alpha_diversity_plot")) for (m in names(extra)) {
  args <- c(list(table = tab, metadata = meta), extra[[m]])
  for (v in c("antes", "despues")) {
    fun <- if (v == "antes") get(paste0(f, "_old")) else get(f)
    p <- suppressMessages(suppressWarnings(do.call(fun, args)))
    save_plot(p, sprintf("%s_%s_%s.png", f, m, v))
  }
  cat(sprintf("%-24s %-9s listo\n", f, m))
}
