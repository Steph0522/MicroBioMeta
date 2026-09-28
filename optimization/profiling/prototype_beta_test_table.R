# prototype_beta_test_table.R
# Hotspot: R/beta_test_table.R#148-150
#   aldex_obj <- ALDEx2::aldex.clr(t(table), mc.samples = 128, ...)   # 128 instancias Monte Carlo
#   clr_samples <- t(ALDEx2::getMonteCarloSample(aldex_obj, 1))        # ...pero solo se usa la 1a
# ~90-95% del tiempo de beta_test_table esta en aldex.clr (rgamma de 79 x 3645 x 128 = 36.9 M
# valores, log2 + apply(clr) sobre las 128 columnas, t(), etc.). adonis2 (999 perm, by="terms") ~2%.
#
# Prototipos:
#  P1 "exacto":  genera EXACTAMENTE los mismos numeros aleatorios que ALDEx2 (misma semilla, mismo
#                orden: rgamma(l*128, a) por muestra, relleno byrow) pero solo calcula la CLR de la
#                1a instancia MC. Distancia bit-identica y MISMO estado del RNG al llegar a adonis2
#                -> tabla identica (incluidos p-valores).
#  P2 "1 instancia": solo genera l valores gamma por muestra (equivale a mc.samples = 1). Resultado
#                estadisticamente equivalente pero NO identico (otra realizacion de Dirichlet).
#  P3 adonis2 con parallel = 4 (solo para mostrar que no compensa: la permutacion no es el cuello).
#
# Uso: cd optimization/profiling; Rscript prototype_beta_test_table.R
source("00_setup_data.R")
N <- as.integer(Sys.getenv("NREPS", "4"))

# ---- replica de la parte de calculo de beta_test_table (sin la tabla ggpubr) ----------------
prep_counts <- function(table, metadata) {
  tax_cols <- grep("taxonomy|taxon|Taxonomy|Taxa", names(table))
  if (length(tax_cols) > 0) table <- table[, -tax_cols[1], drop = FALSE]
  num_cols <- vapply(table, is.numeric, logical(1))
  table <- as.matrix(table[, num_cols, drop = FALSE])
  if (ncol(table) == nrow(metadata)) table <- t(table)
  table                       # muestras x features
}

# CLR "exacta" = 1a instancia MC de ALDEx2::aldex.clr(mc.samples = 128, denom = "all")
clr_first_mc_exact <- function(counts_sxf, mc.samples = 128, prior = 0.5) {
  reads <- t(counts_sxf)                                  # features x muestras (como ALDEx2)
  reads <- reads[rowSums(reads) > 0, , drop = FALSE]      # ALDEx2 quita filas con suma 0
  l <- nrow(reads)
  out <- matrix(NA_real_, ncol(reads), l, dimnames = list(colnames(reads), rownames(reads)))
  for (j in seq_len(ncol(reads))) {
    a <- reads[, j] + prior
    g <- rgamma(l * mc.samples, a)                        # consume el RNG igual que ALDEx2:::rdirichlet
    x <- g[seq_len(l)]                                    # fila 1 (byrow = TRUE) = 1a instancia MC
    p <- x / sum(x)
    lp <- log2(p)
    out[j, ] <- lp - mean(lp)
  }
  out
}
clr_one_mc <- function(counts_sxf, prior = 0.5) clr_first_mc_exact(counts_sxf, mc.samples = 1, prior = prior)

run_adonis <- function(dist_matrix, metadata, formula_str, permutations, parallel = getOption("mc.cores")) {
  as.data.frame(vegan::adonis2(as.formula(paste("dist_matrix ~", formula_str)), data = metadata,
                               permutations = permutations, by = "terms", parallel = parallel))
}

original_core <- function(table, metadata, formula_str, permutations = 999, seed = 123) {
  counts <- prep_counts(table, metadata)
  set.seed(seed)
  aldex_obj <- ALDEx2::aldex.clr(t(counts), mc.samples = 128, denom = "all", verbose = FALSE, useMC = FALSE)
  clr_samples <- t(ALDEx2::getMonteCarloSample(aldex_obj, 1))
  d <- stats::dist(clr_samples, method = "euclidean")
  list(dist = d, tab = run_adonis(d, metadata, formula_str, permutations))
}
p1_core <- function(table, metadata, formula_str, permutations = 999, seed = 123) {
  counts <- prep_counts(table, metadata)
  set.seed(seed)
  d <- stats::dist(clr_first_mc_exact(counts), method = "euclidean")
  list(dist = d, tab = run_adonis(d, metadata, formula_str, permutations))
}
p2_core <- function(table, metadata, formula_str, permutations = 999, seed = 123) {
  counts <- prep_counts(table, metadata)
  set.seed(seed)
  d <- stats::dist(clr_one_mc(counts), method = "euclidean")
  list(dist = d, tab = run_adonis(d, metadata, formula_str, permutations))
}
p3_core <- function(table, metadata, formula_str, permutations = 999, seed = 123) {
  counts <- prep_counts(table, metadata)
  set.seed(seed)
  d <- stats::dist(clr_first_mc_exact(counts), method = "euclidean")
  list(dist = d, tab = run_adonis(d, metadata, formula_str, permutations, parallel = 4))
}

fs <- "Treatment*Type_of_soil"
# ---- Verificacion -----------------------------------------------------------------------------
o  <- suppressMessages(original_core(table_bac, metadata_bacteria, fs))
r1 <- p1_core(table_bac, metadata_bacteria, fs)
r2 <- p2_core(table_bac, metadata_bacteria, fs)
cat("\n== Verificacion P1 (exacto) ==\n")
cat("identical(dist):", identical(as.vector(o$dist), as.vector(r1$dist)),
    "| max|diff| dist:", max(abs(o$dist - r1$dist)), "\n")
cat("identical(tabla adonis2):", identical(o$tab, r1$tab), "\n")
cat("all.equal(tabla adonis2):", isTRUE(all.equal(o$tab, r1$tab)), "\n")
print(cbind(original = o$tab[, c("R2", "F", "Pr(>F)")], P1 = r1$tab[, c("R2", "F", "Pr(>F)")]))
cat("\n== P2 (1 instancia MC; NO identico, solo comparacion) ==\n")
print(cbind(original = o$tab[, c("R2", "F", "Pr(>F)")], P2 = r2$tab[, c("R2", "F", "Pr(>F)")]))
cat("cor(dist original, dist P2):", cor(as.vector(o$dist), as.vector(r2$dist)), "\n")

# ---- P4: CLR determinista = valor esperado de log2(p) bajo Dirichlet(a) --------------------------
# E[log p_i] = digamma(a_i) - digamma(sum a). Es el limite de promediar infinitas instancias MC:
# no depende de la semilla y cuesta ~ lo mismo que una CLR con pseudoconteo.
clr_expected <- function(counts_sxf, prior = 0.5) {
  reads <- t(counts_sxf); reads <- reads[rowSums(reads) > 0, , drop = FALSE]
  el <- digamma(reads + prior)                              # la constante digamma(sum a) se cancela en la CLR
  el <- sweep(el, 2, colMeans(el)) / log(2)
  t(el)
}
p4_core <- function(table, metadata, formula_str, permutations = 999, seed = 123) {
  counts <- prep_counts(table, metadata)
  d <- stats::dist(clr_expected(counts), method = "euclidean")
  set.seed(seed)
  list(dist = d, tab = run_adonis(d, metadata, formula_str, permutations))
}
r4 <- p4_core(table_bac, metadata_bacteria, fs)
# Comparar con el promedio de las 128 instancias MC de ALDEx2 (lo que P4 aproxima)
set.seed(123)
ax <- suppressMessages(ALDEx2::aldex.clr(t(prep_counts(table_bac, metadata_bacteria)), mc.samples = 128,
                                         denom = "all", verbose = FALSE, useMC = FALSE))
clr_mean128 <- t(sapply(ALDEx2::getMonteCarloInstances(ax), rowMeans))
d_mean128 <- stats::dist(clr_mean128)
cat("\n== P4 (CLR esperada, determinista) vs media de 128 instancias MC ==\n")
cat("cor(dist P4, dist media-128MC):", cor(as.vector(r4$dist), as.vector(d_mean128)),
    "| cor(dist P4, dist 1a-instancia original):", cor(as.vector(r4$dist), as.vector(o$dist)), "\n")
set.seed(123); tab_mean128 <- run_adonis(d_mean128, metadata_bacteria, fs, 999)
print(cbind(media128 = tab_mean128[, c("R2", "F", "Pr(>F)")], P4 = r4$tab[, c("R2", "F", "Pr(>F)")]))
# Variabilidad de la implementacion actual: p-valor de la interaccion con distintas semillas
cat("\n== Sensibilidad a la semilla de la implementacion actual (1 instancia de 128) ==\n")
sens <- t(sapply(c(1, 2, 3, 123, 2024), function(s) {
  r <- p1_core(table_bac, metadata_bacteria, fs, seed = s)$tab
  c(seed = s, R2_int = r["Treatment:Type_of_soil", "R2"], p_Treat = r["Treatment", "Pr(>F)"],
    p_Soil = r["Type_of_soil", "Pr(>F)"], p_int = r["Treatment:Type_of_soil", "Pr(>F)"])
}))
print(sens)
write.csv(sens, "results_prototype_beta_seed_sensitivity.csv", row.names = FALSE)

# ---- Benchmark ---------------------------------------------------------------------------------
bm <- function(f, n = N) { t <- numeric(n); for (i in seq_len(n)) { gc(FALSE)
  t[i] <- system.time(suppressMessages(f()))[["elapsed"]] }; c(first = t[1], median = median(t[-1])) }
res <- rbind(
  "beta_test_table() publica (original)"   = bm(function() beta_test_table(table_bac, metadata_bacteria, fs,
                                               method = "compositional", permutations = 999)),
  "nucleo original (aldex.clr 128 MC)"     = bm(function() original_core(table_bac, metadata_bacteria, fs)),
  "P1 exacto (misma RNG, solo 1a CLR)"     = bm(function() p1_core(table_bac, metadata_bacteria, fs)),
  "P2 1 instancia MC (no identico)"        = bm(function() p2_core(table_bac, metadata_bacteria, fs)),
  "P3 = P1 + adonis2(parallel = 4)"        = bm(function() p3_core(table_bac, metadata_bacteria, fs)),
  "P4 CLR esperada digamma (determinista)" = bm(function() p4_core(table_bac, metadata_bacteria, fs)),
  "solo aldex.clr 128 MC"                  = bm(function() { set.seed(123); ALDEx2::aldex.clr(t(prep_counts(table_bac, metadata_bacteria)),
                                               mc.samples = 128, denom = "all", verbose = FALSE, useMC = FALSE) }),
  "solo rgamma(l*128) por muestra"         = bm(function() { set.seed(123); clr_first_mc_exact(prep_counts(table_bac, metadata_bacteria)) }),
  "solo adonis2 ~T*S 999 by=terms"         = bm(function() run_adonis(r1$dist, metadata_bacteria, fs, 999))
)
print(round(res, 3))
write.csv(data.frame(variant = rownames(res), round(res, 3)), "results_prototype_beta_test_table.csv", row.names = FALSE)
