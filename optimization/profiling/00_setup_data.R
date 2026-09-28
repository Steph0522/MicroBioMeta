# 00_setup_data.R
# Preparacion de datos IDENTICA (verbatim) al script original de la colaboradora
# (optimization/original_comparison_script.R). Se "sourcea" desde los demas scripts.
# Ejecutar con el directorio de trabajo = optimization/profiling (PKG_ROOT = ../..).
# Se abre un dispositivo pdf(NULL) para que nada se escriba en Rplots.pdf del paquete.
PKG_ROOT <- normalizePath(file.path("..", ".."))
if (!file.exists(file.path(PKG_ROOT, "DESCRIPTION"))) stop("Ejecutar desde optimization/profiling")
grDevices::pdf(NULL)

# Cache en disco de rbiom: rbiom guarda resultados en tempdir() y las repeticiones
# devolverian el valor cacheado. rbiom_fresh_cache() apunta a un directorio vacio nuevo
# antes de cada repeticion (desactivarlo con "FALSE" rompe varias funciones de rbiom 3.1.0).
# RBIOM_CACHE=1 deja la cache normal activada.
rbiom_fresh_cache <- function() {
  if (Sys.getenv("RBIOM_CACHE", "0") != "1") options(rbiom.cache_dir = tempfile("rbiomcache"))
}
rbiom_fresh_cache()
HAS_RBIOM <- requireNamespace("rbiom", quietly = TRUE)
# Mismo orden de carga que el script original: load_all, rbiom, phyloseq, tidyverse
# (asi phyloseq::taxa_sums enmascara a rbiom::taxa_sums, como en su sesion).
suppressPackageStartupMessages({
  devtools::load_all(PKG_ROOT, quiet = TRUE)
  if (HAS_RBIOM) library(rbiom)
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

###IMPORTAR LOS DATOS ####
bacteria_table <- read.delim(
  system.file("extdata", "table_bacteria.txt", package = "MicroBioMeta"),
  row.names = 1, check.names = FALSE
)

taxonomy_bacteria <- read.delim(
  system.file("extdata", "taxonomy_bacteria.txt", package = "MicroBioMeta"),
  check.names = FALSE
) %>%
  rename(taxonomy = Taxon) %>%
  dplyr::select(-Confidence) %>%
  column_to_rownames(var = "Feature.ID")

bacteria_metadata <- read.delim(
  system.file("extdata", "metadata_bacterias.txt", package = "MicroBioMeta"),
  check.names = FALSE
) %>%
  filter(Month == "2")

samples_bac <- bacteria_metadata$SAMPLEID[
  bacteria_metadata$SAMPLEID %in% colnames(bacteria_table)]

table_bacteria <- bacteria_table[, samples_bac]

metadata_bacteria <- bacteria_metadata[
  bacteria_metadata$SAMPLEID %in% samples_bac, ]

table_bac <- merge_feature_taxonomy(table = table_bacteria,
                                    taxonomy = taxonomy_bacteria)

###CREAR ARCHIVO PHYLOSEQ ####
otu_bac <- otu_table(as.matrix(table_bacteria), taxa_are_rows = TRUE)
tax_bac <- tax_table(as.matrix(taxonomy_bacteria))
meta_bac <- sample_data(data.frame(metadata_bacteria,
                                   row.names = metadata_bacteria$SAMPLEID))
ps_bacteria <- phyloseq(otu_bac, tax_bac, meta_bac)

####HACER COMPATIBLE LA TAXONOMIA ####
tax_df <- as.data.frame(as(tax_table(ps_bacteria), "matrix"))
nombre_columna <- colnames(tax_df)[1]
tax_limpia <- tax_df %>%
  separate(!!sym(nombre_columna),
           into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";\\s*",
           fill = "right") %>%
  mutate(across(everything(), ~ gsub("^[a-z]__", "", .)))
tax_limpia[tax_limpia == ""] <- NA
rownames(tax_limpia) <- rownames(tax_df)
tax_table(ps_bacteria) <- tax_table(as.matrix(tax_limpia))

if (HAS_RBIOM) data <- as_rbiom(ps_bacteria)

env_table_bac <- metadata_bacteria %>%
  dplyr::select(SAMPLEID, pH:Arbus_per) %>%
  remove_rownames() %>%
  column_to_rownames(var = "SAMPLEID")

# ---- Utilidades de medicion ----------------------------------------------
# Renderizado "justo": dibujar el objeto en un dispositivo pdf nulo.
render_null <- function(obj) {
  grDevices::pdf(NULL, width = 10, height = 8)
  on.exit(grDevices::dev.off())
  if (inherits(obj, c("gg", "ggplot"))) print(obj)
  else if (inherits(obj, "grob") || inherits(obj, "gTree")) { grid::grid.newpage(); grid::grid.draw(obj) }
  else if (methods::is(obj, "HeatmapList") || methods::is(obj, "Heatmap")) ComplexHeatmap::draw(obj)
  else print(obj)
  invisible(NULL)
}

# Repite expr n veces, devuelve primera ejecucion y mediana de las restantes
time_reps <- function(expr_fun, n = 5) {
  t <- numeric(n)
  for (i in seq_len(n)) {
    gc(FALSE); rbiom_fresh_cache()
    t[i] <- system.time(expr_fun())[["elapsed"]]
  }
  c(first = t[1], median_rest = stats::median(t[-1]), min = min(t[-1]), max = max(t[-1]))
}

# Crea una VARIANTE de una funcion del paquete leyendo R/<file> (sin modificarlo), aplicando
# sustituciones de texto literales (cada una debe coincidir exactamente 1 vez) y evaluandola
# con el namespace del paquete como entorno padre. Asi los prototipos siguen el codigo real.
patch_fun <- function(file, fname, subs = list()) {
  src <- paste(readLines(file.path(PKG_ROOT, "R", file), warn = FALSE), collapse = "\n")
  for (s in subs) {
    n <- lengths(regmatches(src, gregexpr(s[1], src, fixed = TRUE)))
    if (n != 1) stop("Sustitucion no unica (", n, "): ", substr(s[1], 1, 60))
    src <- sub(s[1], s[2], src, fixed = TRUE)
  }
  env <- new.env(parent = asNamespace("MicroBioMeta"))
  eval(parse(text = src, keep.source = FALSE), envir = env)
  get(fname, envir = env)
}

# Renderiza a PNG y devuelve la matriz de pixeles (para comparar salidas graficas)
render_png_pixels <- function(obj, file = tempfile(fileext = ".png"), width = 1400, height = 1000) {
  grDevices::png(file, width = width, height = height, res = 110)
  if (inherits(obj, c("gg", "ggplot"))) print(obj) else { grid::grid.newpage(); grid::grid.draw(obj) }
  grDevices::dev.off()
  png::readPNG(file)
}
