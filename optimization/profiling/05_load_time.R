# 05_load_time.R
# Coste de arranque: cada medicion se hace en un proceso R NUEVO (Rscript), 3 repeticiones, mediana.
# Mide: library(MicroBioMeta) instalado, devtools::load_all(), y loadNamespace() de cada Import pesado.
# Uso: cd optimization/profiling; Rscript 05_load_time.R
PKG_ROOT <- normalizePath(file.path("..", ".."), winslash = "/")
RSCRIPT <- file.path(R.home("bin"), "Rscript.exe")
NREP <- as.integer(Sys.getenv("NREPS", "3"))

run_child <- function(code) {
  f <- tempfile(fileext = ".R")
  writeLines(c("t0 <- proc.time()[['elapsed']]",
               "suppressWarnings(suppressPackageStartupMessages({", code, "}))",
               "cat('ELAPSED=', proc.time()[['elapsed']] - t0, '\\n', sep = '')",
               "cat('NNS=', length(loadedNamespaces()), '\\n', sep = '')"), f)
  out <- suppressWarnings(system2(RSCRIPT, shQuote(f), stdout = TRUE, stderr = TRUE))
  unlink(f)
  el  <- as.numeric(sub("ELAPSED=", "", grep("^ELAPSED=", out, value = TRUE)))
  nns <- as.integer(sub("NNS=", "", grep("^NNS=", out, value = TRUE)))
  if (!length(el)) { message(paste(tail(out, 5), collapse = "\n")); return(c(NA, NA)) }
  c(el, nns)
}
measure <- function(code) {
  r <- replicate(NREP, run_child(code))
  c(median_s = median(r[1, ]), min_s = min(r[1, ]), n_namespaces = r[2, 1])
}

desc <- read.dcf(file.path(PKG_ROOT, "DESCRIPTION"))
imports <- trimws(sub("\\(.*", "", strsplit(desc[, "Imports"], ",")[[1]]))
imports <- imports[imports != ""]
base_pk <- c("utils", "stats", "grDevices", "grid", "tools")

cases <- c(
  "R vacio (baseline)"                    = "invisible(NULL)",
  "library(MicroBioMeta) [instalado]"     = "library(MicroBioMeta)",
  "loadNamespace('devtools')"             = "loadNamespace('devtools')",
  "devtools::load_all() [como la colaboradora]" = sprintf("devtools::load_all('%s', quiet = TRUE)", PKG_ROOT),
  "todos los Imports (loadNamespace)"     = paste(sprintf("loadNamespace('%s')", setdiff(imports, base_pk)), collapse = "; ")
)
for (p in setdiff(imports, base_pk)) cases[paste0("  import: ", p)] <- sprintf("loadNamespace('%s')", p)
cases["library(phyloseq) [competidor]"] <- "library(phyloseq)"
if (requireNamespace("rbiom", quietly = TRUE)) cases["library(rbiom) [competidor]"] <- "library(rbiom)"

res <- do.call(rbind, lapply(names(cases), function(nm) {
  r <- measure(cases[[nm]]); cat(sprintf("%-45s %6.2f s  (ns=%d)\n", nm, r[1], as.integer(r[3])))
  data.frame(case = nm, t(round(r, 2)), check.names = FALSE)
}))
base <- res$median_s[1]
res$minus_baseline <- round(res$median_s - base, 2)
res <- res[c(seq_len(5), 5 + order(-res$median_s[-(1:5)])), ]
print(res, row.names = FALSE)
cat("\nVersion instalada:", as.character(packageVersion("MicroBioMeta")), " | DESCRIPTION fuente:", desc[, "Version"], "\n")
write.csv(res, "results_05_load_time.csv", row.names = FALSE)

# Primera llamada tras library() (lazy-load de sub-dependencias, JIT) vs load_all
first_call <- function(loader) sprintf('%s
t <- read.delim(system.file("extdata","table_bacteria.txt",package="MicroBioMeta"), row.names=1, check.names=FALSE)
tx <- read.delim(system.file("extdata","taxonomy_bacteria.txt",package="MicroBioMeta"), check.names=FALSE)
m <- read.delim(system.file("extdata","metadata_bacterias.txt",package="MicroBioMeta"), check.names=FALSE)
m <- m[m$Month == "2", ]; s <- m$SAMPLEID[m$SAMPLEID %%in%% colnames(t)]; t <- t[, s]; m <- m[m$SAMPLEID %%in%% s, ]
tx <- data.frame(taxonomy = tx$Taxon, row.names = tx$Feature.ID)
tb <- merge_feature_taxonomy(table = t, taxonomy = tx)
grDevices::pdf(NULL)
t1 <- system.time(abundance_heatmap_plot(table = tb, metadata = m, top_n = 20, show_column_names = FALSE, condition1 = "Treatment", condition2 = "Type_of_soil", feature_prefix = "ASV"))[["elapsed"]]
t2 <- system.time(abundance_heatmap_plot(table = tb, metadata = m, top_n = 20, show_column_names = FALSE, condition1 = "Treatment", condition2 = "Type_of_soil", feature_prefix = "ASV"))[["elapsed"]]
t3 <- system.time(abundance_heatmap_plot(table = tb, metadata = m, top_n = 20, show_column_names = FALSE, condition1 = "Treatment", condition2 = "Type_of_soil", feature_prefix = "ASV"))[["elapsed"]]
cat("HM=", t1, ",", t2, ",", t3, "\\n", sep = "")', loader)
fc <- function(loader) { f <- tempfile(fileext = ".R"); writeLines(first_call(loader), f)
  o <- suppressWarnings(system2(RSCRIPT, shQuote(f), stdout = TRUE, stderr = TRUE)); unlink(f)
  grep("^HM=", o, value = TRUE) }
cat("\nabundance_heatmap_plot llamadas 1,2,3 en sesion nueva:\n")
cat("  con library(MicroBioMeta) (byte-compilado):", fc("suppressPackageStartupMessages(library(MicroBioMeta))"), "\n")
cat("  con devtools::load_all() (sin byte-compilar):", fc(sprintf("suppressPackageStartupMessages(devtools::load_all('%s', quiet = TRUE))", PKG_ROOT)), "\n")
