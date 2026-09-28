# 01_original_timings.R
# Reproduce las mediciones de la colaboradora (system.time sobre la llamada,
# SIN renderizar el objeto devuelto). Primera ejecucion separada + mediana de N reps.
# Uso: cd optimization/profiling; Rscript 01_original_timings.R
#      RBIOM_CACHE=1 Rscript 01_original_timings.R  (deja activa la cache en disco de rbiom, como en una sesion normal)
t_load <- system.time(source("00_setup_data.R"))[["elapsed"]]
cat("Tiempo setup (load_all + librerias + datos):", t_load, "s\n")
cat("Dimensiones table_bac:", nrow(table_bac), "features x", ncol(table_bac) - 1, "muestras\n")
cat("Muestras metadata:", nrow(metadata_bacteria), "\n")
cat("Niveles Treatment:", paste(unique(metadata_bacteria$Treatment), collapse = ","),
    "| Type_of_soil:", paste(unique(metadata_bacteria$Type_of_soil), collapse = ","), "\n")
cat("Features con conteo > 0:", sum(rowSums(table_bacteria) > 0), "\n")

N <- as.integer(Sys.getenv("NREPS", "4"))  # 1 primera + 3 reps
calls <- list(
  mbm_abundance_bar_plot = function() abundance_bar_plot(table = table_bac, metadata = metadata_bacteria,
      taxonomy_db = "silva", level = "phylum", top_n = 15, x_col = "Type_of_soil",
      facet_col = "Treatment", add_remained = TRUE, label = "Phylum", x_label_angle = 45, save_table = FALSE),
  phyloseq_plot_bar = function() {
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
    p <- plot_bar(ps_rel, x = "Type_of_soil", fill = "Phylum", facet_grid = ~Treatment)
    p
  },
  rbiom_taxa_stacked = function() { data3 <- rbiom::biom_relativize(data)
    rbiom::taxa_stacked(data3, rank = "Phylum", taxa = 10, label.by = "Type_of_soil", facet.by = "Treatment") },
  mbm_abundance_heatmap_plot = function() abundance_heatmap_plot(table = table_bac, metadata = metadata_bacteria,
      top_n = 20, show_column_names = FALSE, condition1 = "Treatment", condition2 = "Type_of_soil",
      feature_prefix = "ASV"),
  phyloseq_plot_heatmap = function() {
    ps_bacteria2 <- prune_taxa(taxa_sums(ps_bacteria) > 0, ps_bacteria)
    ps_bacteria3 <- tax_glom(ps_bacteria2, taxrank = "Phylum")
    plot_heatmap(ps_bacteria3, taxa.label = "Phylum", sample.label = "Type_of_soil")
  },
  rbiom_taxa_heatmap = function() rbiom::taxa_heatmap(data, taxa = 20, tracks = c("Treatment", "Type_of_soil"), grid = blues9),
  mbm_alpha_diversity_plot = function() alpha_diversity_plot(table = table_bac, metadata = metadata_bacteria,
      x_col = "Treatment", fill_col = "Treatment", facet_by = "Type_of_soil",
      facet_orientation = "horizontal", save_table = FALSE),
  phyloseq_plot_richness = function() {
    p <- plot_richness(ps_bacteria, x = "Type_of_soil", color = "Treatment", measures = c("Chao1", "Shannon", "Simpson"))
    p + facet_grid(Type_of_soil ~ variable)
  },
  rbiom_adiv_boxplot = function() rbiom::adiv_boxplot(data, x = "Treatment", facet.by = "Type_of_soil",
      adiv = c("otu", "shannon", "simpson")),
  mbm_beta_ord_plot = function() beta_ord_plot(table = table_bac, metadata = metadata_bacteria,
      distance = "bray", ordination = "PCoA", group_col = "Treatment", shape_col = "Type_of_soil"),
  phyloseq_ordinate_plot = function() { GP.ord <- ordinate(ps_bacteria, "PCoA", "bray")
    p1 <- plot_ordination(ps_bacteria, GP.ord, type = "SAMPLEID", color = "Treatment", title = "PCoA", shape = "Type_of_soil")
    print(p1) },  # (la colaboradora hizo print() aqui)
  rbiom_bdiv_ord_plot = function() rbiom::bdiv_ord_plot(data, bdiv = "bray", ord = "PCoA", layers = "petm",
      stat.by = "Treatment", facet.by = NULL, colors = TRUE, shapes = TRUE, tree = NULL,
      test = "adonis2", seed = 0, permutations = 999),
  mbm_beta_test_table = function() beta_test_table(table = table_bac, metadata = metadata_bacteria,
      formula_str = "Treatment*Type_of_soil", method = "compositional", test = "permanova", permutations = 999),
  rbiom_bdiv_stats = function() rbiom::bdiv_stats(data, stat.by = "Treatment"),
  rbiom_distmat_stats = function() { dm <- rbiom::bdiv_distmat(data, bdiv = "bray")
    sex_vector <- dplyr::pull(data, "Treatment")
    rbiom::distmat_stats(dm, groups = sex_vector, test = "adonis2") },
  mbm_cca_rda_biplot = function() cca_rda_biplot(table = table_bac, env_data = env_table_bac,
      metadata = metadata_bacteria, env_vars = c("pH", "TN", "WHC", "EC", "Clay"), analysis = "CCA",
      show_all_env_vectors = TRUE, group_col = "Type_of_soil", scale_arrows = 3),
  phyloseq_cca = function() {
    cca_resultado <- ordinate(ps_bacteria, method = "CCA", formula = ~ pH + TN + WHC + EC + Clay)
    p <- plot_ordination(ps_bacteria, cca_resultado, type = "SAMPLEID", color = "Type_of_soil") +
      theme_minimal() + geom_point(size = 3)
    flechas <- as.data.frame(cca_resultado$CCA$biplot); flechas$Variable <- rownames(flechas)
    p + geom_segment(data = flechas, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), inherit.aes = FALSE) +
      geom_text(data = flechas, aes(x = CCA1, y = CCA2, label = Variable), inherit.aes = FALSE)
  }
)
if (!HAS_RBIOM) calls <- calls[!grepl("^rbiom", names(calls))]

res <- do.call(rbind, lapply(names(calls), function(nm) {
  n <- if (nm == "mbm_beta_test_table") min(N, 4) else N
  r <- tryCatch(suppressWarnings(suppressMessages(time_reps(calls[[nm]], n))),
                error = function(e) { message(nm, ": ", conditionMessage(e)); c(first = NA, median_rest = NA, min = NA, max = NA) })
  cat(sprintf("%-28s first=%6.2f  median=%6.2f\n", nm, r[1], r[2]))
  data.frame(call = nm, t(r))
}))
print(res, row.names = FALSE)
write.csv(res, if (Sys.getenv("RBIOM_CACHE", "0") == "1") "results_01_original_timings_rbiom_cache_on.csv" else "results_01_original_timings.csv", row.names = FALSE)
