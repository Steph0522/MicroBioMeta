# 02_fair_timings.R
# Comparacion "justa": tiempo de CONSTRUIR el objeto + tiempo de RENDERIZARLO
# (print/grid.draw en un dispositivo pdf(NULL)) para todos los paquetes.
# Tambien descompone beta_test_table en trabajo equivalente al de rbiom.
# Uso: cd optimization/profiling; Rscript 02_fair_timings.R
source("00_setup_data.R")
N <- as.integer(Sys.getenv("NREPS", "4"))  # 1 warm-up + 3 reps

build_render <- function(build_fun, n = N, render = TRUE) {
  tb <- tr <- numeric(n)
  for (i in seq_len(n)) {
    gc(FALSE); rbiom_fresh_cache()
    tb[i] <- system.time(obj <- suppressWarnings(suppressMessages(build_fun())))[["elapsed"]]
    tr[i] <- if (render) system.time(render_null(obj))[["elapsed"]] else NA
  }
  c(build_first = tb[1], render_first = tr[1],
    build_med = median(tb[-1]), render_med = median(tr[-1]), total_med = median(tb[-1] + tr[-1]))
}

plots <- list(
  "bar | MicroBioMeta abundance_bar_plot" = function() abundance_bar_plot(table = table_bac, metadata = metadata_bacteria,
      taxonomy_db = "silva", level = "phylum", top_n = 15, x_col = "Type_of_soil", facet_col = "Treatment",
      add_remained = TRUE, label = "Phylum", x_label_angle = 45, save_table = FALSE),
  "bar | phyloseq plot_bar (+prep)" = function() {
    ps_filo <- tax_glom(ps_bacteria, taxrank = "Phylum")
    top15_filos <- names(sort(taxa_sums(ps_filo), decreasing = TRUE)[1:15])
    tax_mat <- as(tax_table(ps_filo), "matrix")
    tax_mat[!(rownames(tax_mat) %in% top15_filos), "Phylum"] <- "Otros"
    tax_table(ps_filo) <- tax_table(tax_mat)
    sample_data(ps_filo)$Group <- paste(sample_data(ps_filo)$Type_of_soil, sample_data(ps_filo)$Treatment, sep = "_")
    ps_merged <- merge_samples(ps_filo, "Group")
    df_meta <- data.frame(Group = sample_names(ps_merged))
    df_meta <- separate(df_meta, Group, into = c("Type_of_soil", "Treatment"), sep = "_", remove = FALSE)
    rownames(df_meta) <- df_meta$Group
    sample_data(ps_merged) <- sample_data(df_meta)
    ps_rel <- transform_sample_counts(ps_merged, function(x) (x / sum(x)))
    plot_bar(ps_rel, x = "Type_of_soil", fill = "Phylum", facet_grid = ~Treatment)
  },
  "bar | rbiom taxa_stacked" = function() rbiom::taxa_stacked(rbiom::biom_relativize(data), rank = "Phylum",
      taxa = 10, label.by = "Type_of_soil", facet.by = "Treatment"),
  "heatmap | MicroBioMeta abundance_heatmap_plot (dibuja DENTRO)" = function() abundance_heatmap_plot(
      table = table_bac, metadata = metadata_bacteria, top_n = 20, show_column_names = FALSE,
      condition1 = "Treatment", condition2 = "Type_of_soil", feature_prefix = "ASV"),
  "heatmap | phyloseq plot_heatmap (+prune/glom)" = function() {
    ps2 <- prune_taxa(taxa_sums(ps_bacteria) > 0, ps_bacteria)
    plot_heatmap(tax_glom(ps2, taxrank = "Phylum"), taxa.label = "Phylum", sample.label = "Type_of_soil")
  },
  "heatmap | rbiom taxa_heatmap" = function() rbiom::taxa_heatmap(data, taxa = 20,
      tracks = c("Treatment", "Type_of_soil"), grid = blues9),
  "alpha | MicroBioMeta alpha_diversity_plot" = function() alpha_diversity_plot(table = table_bac,
      metadata = metadata_bacteria, x_col = "Treatment", fill_col = "Treatment", facet_by = "Type_of_soil",
      facet_orientation = "horizontal", save_table = FALSE),
  "alpha | phyloseq plot_richness" = function() plot_richness(ps_bacteria, x = "Type_of_soil", color = "Treatment",
      measures = c("Chao1", "Shannon", "Simpson")) + facet_grid(Type_of_soil ~ variable),
  "alpha | rbiom adiv_boxplot" = function() rbiom::adiv_boxplot(data, x = "Treatment", facet.by = "Type_of_soil",
      adiv = c("otu", "shannon", "simpson")),
  "ord | MicroBioMeta beta_ord_plot" = function() beta_ord_plot(table = table_bac, metadata = metadata_bacteria,
      distance = "bray", ordination = "PCoA", group_col = "Treatment", shape_col = "Type_of_soil"),
  "ord | phyloseq ordinate+plot_ordination (sin print)" = function() plot_ordination(ps_bacteria,
      ordinate(ps_bacteria, "PCoA", "bray"), type = "SAMPLEID", color = "Treatment", title = "PCoA",
      shape = "Type_of_soil"),
  "ord | rbiom bdiv_ord_plot (adonis2 999)" = function() rbiom::bdiv_ord_plot(data, bdiv = "bray", ord = "PCoA",
      layers = "petm", stat.by = "Treatment", facet.by = NULL, colors = TRUE, shapes = TRUE, tree = NULL,
      test = "adonis2", seed = 0, permutations = 999),
  "cca | MicroBioMeta cca_rda_biplot" = function() cca_rda_biplot(table = table_bac, env_data = env_table_bac,
      metadata = metadata_bacteria, env_vars = c("pH", "TN", "WHC", "EC", "Clay"), analysis = "CCA",
      show_all_env_vectors = TRUE, group_col = "Type_of_soil", scale_arrows = 3),
  "cca | phyloseq ordinate CCA + plot" = function() {
    cr <- ordinate(ps_bacteria, method = "CCA", formula = ~ pH + TN + WHC + EC + Clay)
    fl <- as.data.frame(cr$CCA$biplot); fl$Variable <- rownames(fl)
    plot_ordination(ps_bacteria, cr, type = "SAMPLEID", color = "Type_of_soil") + theme_minimal() +
      geom_point(size = 3) +
      geom_segment(data = fl, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), inherit.aes = FALSE) +
      geom_text(data = fl, aes(x = CCA1, y = CCA2, label = Variable), inherit.aes = FALSE)
  }
)
if (!HAS_RBIOM) plots <- plots[!grepl("rbiom", names(plots))]

res_plots <- do.call(rbind, lapply(names(plots), function(nm) {
  r <- tryCatch(build_render(plots[[nm]]), error = function(e) { message(nm, ": ", conditionMessage(e)); rep(NA, 5) })
  cat(sprintf("%-62s build=%5.2f render=%5.2f total=%5.2f\n", nm, r[3], r[4], r[5]))
  data.frame(call = nm, t(round(r, 3)), check.names = FALSE)
}))
write.csv(res_plots, "results_02_fair_build_render.csv", row.names = FALSE)

# --- Heatmap MicroBioMeta: coste del dibujo interno en un dispositivo "real" (png) ----------
# En RStudio el grid.newpage()+grid.draw() de abundance_heatmap_plot.R#394-395 va al panel Plots.
dev_cost <- function(fun, n = N) {
  t <- numeric(n)
  for (i in seq_len(n)) {
    f <- tempfile(fileext = ".png"); grDevices::png(f, width = 1400, height = 1000, res = 110)
    t[i] <- system.time(suppressWarnings(suppressMessages(fun())))[["elapsed"]]
    grDevices::dev.off(); unlink(f)
  }
  c(first = t[1], median = median(t[-1]))
}
hm_png <- dev_cost(plots[["heatmap | MicroBioMeta abundance_heatmap_plot (dibuja DENTRO)"]])
cat("Heatmap MBM con dispositivo png activo: first=", hm_png[1], " median=", hm_png[2], "\n")

# --- beta_test_table: descomposicion del trabajo --------------------------------------------
dm_bray <- vegan::vegdist(t(table_bacteria), "bray")
beta_calls <- list(
  "MBM beta_test_table compositional ~Treatment*Type_of_soil 999 (original)" = function()
    beta_test_table(table_bac, metadata_bacteria, "Treatment*Type_of_soil", method = "compositional", permutations = 999),
  "MBM beta_test_table compositional ~Treatment 999" = function()
    beta_test_table(table_bac, metadata_bacteria, "Treatment", method = "compositional", permutations = 999),
  "MBM beta_test_table bray ~Treatment*Type_of_soil 999" = function()
    beta_test_table(table_bac, metadata_bacteria, "Treatment*Type_of_soil", method = "bray", permutations = 999),
  "MBM beta_test_table bray ~Treatment 999 (= trabajo rbiom)" = function()
    beta_test_table(table_bac, metadata_bacteria, "Treatment", method = "bray", permutations = 999),
  "solo ALDEx2::aldex.clr(mc.samples=128)" = function() { set.seed(123)
    tb <- as.matrix(table_bac[, colnames(table_bac) != "taxonomy"])
    ALDEx2::aldex.clr(tb, mc.samples = 128, denom = "all", verbose = FALSE, useMC = FALSE) },
  "solo vegan::adonis2 bray ~Treatment*Type_of_soil 999 by=terms" = function()
    vegan::adonis2(dm_bray ~ Treatment * Type_of_soil, data = metadata_bacteria, permutations = 999, by = "terms"),
  "solo vegan::adonis2 bray ~Treatment 999" = function()
    vegan::adonis2(dm_bray ~ Treatment, data = metadata_bacteria, permutations = 999),
  "rbiom distmat_stats (bdiv_distmat bray + adonis2 ~Treatment 999)" = function() {
    dm <- rbiom::bdiv_distmat(data, bdiv = "bray")
    rbiom::distmat_stats(dm, groups = dplyr::pull(data, "Treatment"), test = "adonis2") }
)
if (!HAS_RBIOM) beta_calls <- beta_calls[!grepl("rbiom", names(beta_calls))]
res_beta <- do.call(rbind, lapply(names(beta_calls), function(nm) {
  is_tab <- grepl("^MBM", nm)
  r <- build_render(beta_calls[[nm]], render = is_tab)
  cat(sprintf("%-75s build=%6.2f render=%5.2f\n", nm, r[3], r[4]))
  data.frame(call = nm, t(round(r, 3)), check.names = FALSE)
}))
write.csv(res_beta, "results_02_beta_decomposition.csv", row.names = FALSE)
writeLines(c(sprintf("heatmap_png_first=%.3f", hm_png[1]), sprintf("heatmap_png_median=%.3f", hm_png[2])),
           "results_02_heatmap_png.txt")
