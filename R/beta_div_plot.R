#' Beta diversity plot with multiple distance and ordination methods
#'
#' This function computes beta diversity using several distance metrics and ordination methods
#' (PCA, PCoA, NMDS). It requires an abundance table with taxonomy, metadata, and allows customization
#' of color and shape aesthetics. It also supports compositional transformation via ALDEx2.
#'
#' @param table A data frame with abundances. The last column must contain taxonomy information, regardless of its name.
#' @param metadata A data frame with sample metadata. The first column must contain the sample IDs, regardless of the column name.
#' @param distance Distance method: one of "euclidean", "bray", "jaccard", "sorensen", or "compositional" (default).
#' @param ordination Ordination method: one of "PCA" (default), "PCoA", or "NMDS".
#' @param color_by Name of the column in `metadata` used to color points.
#' @param shape_by (Optional) Name of the column in `metadata` used to shape points.
#' @param n_taxa Number of top contributing taxa to display as arrows in PCA (default = 5).
#'
#' @return A `ggplot2` object representing the beta diversity ordination plot.
#' @export
#'
#' @importFrom vegan vegdist metaMDS
#' @importFrom stats cmdscale prcomp dist
#' @importFrom ALDEx2 aldex.clr getMonteCarloSample
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_text arrow theme_minimal ggtitle
#' @importFrom ggrepel geom_label_repel
#' @importFrom dplyr inner_join mutate rename
#' @importFrom tidyr everything
#' @importFrom stringr str_extract
#' @importFrom tibble rownames_to_column

beta_div_plot <- function(table, metadata, distance = "compositional",
                          ordination = "PCA", color_by, shape_by = NULL,
                          n_taxa = 5) {
  
  requireNamespace("vegan")
  requireNamespace("ggplot2")
  requireNamespace("ggrepel")
  requireNamespace("ALDEx2")
  requireNamespace("dplyr")
  requireNamespace("stringr")
  
  if (ordination == "PCA" && distance != "compositional") {
    stop("PCA solo está disponible con distancia 'compositional' (Aitchison). Usa otra combinación o cambia a 'PCoA'.")
  }
  
  tax_col <- names(table)[ncol(table)]
  taxonomy <- table[[tax_col]]
  abund_table <- table[, -ncol(table)]
  
  if ("Feature.ID" %in% names(abund_table)) {
    feature_ids <- abund_table$Feature.ID
    abund_table <- abund_table[, !(names(abund_table) == "Feature.ID")]
  } else {
    feature_ids <- rownames(abund_table)
  }
  rownames(abund_table) <- feature_ids
  
  otu_table <- as.data.frame(lapply(abund_table, as.numeric))
  
  metadata_ids <- trimws(as.character(metadata[[1]]))
  sample_ids <- trimws(colnames(otu_table))
  
  muestras_tabla_no_en_metadata <- setdiff(sample_ids, metadata_ids)
  muestras_metadata_no_en_tabla <- setdiff(metadata_ids, sample_ids)
  
  if(length(muestras_tabla_no_en_metadata) > 0 | length(muestras_metadata_no_en_tabla) > 0) {
    warning("Diferencias en nombres de muestras detectadas:")
    if(length(muestras_tabla_no_en_metadata) > 0) {
      warning(paste("Muestras en tabla no en metadata:", paste(muestras_tabla_no_en_metadata, collapse = ", ")))
    }
    if(length(muestras_metadata_no_en_tabla) > 0) {
      warning(paste("Muestras en metadata no en tabla:", paste(muestras_metadata_no_en_tabla, collapse = ", ")))
    }
  }
  
  colnames(otu_table) <- sample_ids
  metadata[[1]] <- metadata_ids
  
  common_samples <- intersect(sample_ids, metadata_ids)
  if (length(common_samples) == 0) stop("No matching sample names between table and metadata.")
  otu_table <- otu_table[, common_samples]
  metadata <- metadata[metadata_ids %in% common_samples, ]
  
  if (distance == "compositional") {
    set.seed(123)
    aldex_obj <- ALDEx2::aldex.clr(t(otu_table), mc.samples = 128,
                                   denom = "all", verbose = FALSE, useMC = FALSE)
    otu_trans <- t(ALDEx2::getMonteCarloSample(aldex_obj, 1))
    dist_matrix <- dist(otu_trans, method = "euclidean")
  } else {
    dist_matrix <- vegan::vegdist(t(otu_table), method = distance)
    otu_trans <- NULL
  }
  
  expl_var <- NULL
  ord_res <- switch(ordination,
                    "PCA" = {
                      pca_input <- if (!is.null(otu_trans)) otu_trans else t(otu_table)
                      res <- prcomp(pca_input)
                      expl_var <<- round(100 * summary(res)$importance[2, 1:2], 1)
                      res
                    },
                    "PCoA" = {
                      res <- cmdscale(dist_matrix, eig = TRUE, k = 2)
                      eigs <- res$eig
                      expl_var <<- round(100 * eigs[1:2] / sum(eigs[eigs > 0]), 1)
                      res
                    },
                    "NMDS" = vegan::metaMDS(dist_matrix, k = 2, trymax = 100),
                    stop("Invalid ordination method."))
  
  ord_df <- switch(ordination,
                   "PCA" = as.data.frame(ord_res$x),
                   "PCoA" = as.data.frame(ord_res$points),
                   "NMDS" = as.data.frame(ord_res$points))
  ord_df$SampleID <- rownames(ord_df)
  
  colnames(metadata)[1] <- "SampleID"
  merged <- dplyr::inner_join(ord_df, metadata, by = "SampleID")
  if (nrow(merged) == 0) stop("Ninguna muestra en común entre la tabla de abundancia y el metadata. Verifica que los nombres coincidan.")
  
  x_lab <- if (!is.null(expl_var)) paste0(names(ord_df)[1], " (", expl_var[1], "%)") else names(ord_df)[1]
  y_lab <- if (!is.null(expl_var)) paste0(names(ord_df)[2], " (", expl_var[2], "%)") else names(ord_df)[2]
  
  p <- ggplot2::ggplot(merged, aes_string(x = names(ord_df)[1],
                                          y = names(ord_df)[2],
                                          color = color_by,
                                          shape = shape_by)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::theme_linedraw() +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::theme(
      axis.text = element_text(colour = "black", size = 12),
      axis.title = element_text(colour = "black", size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12),
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ggplot2::ggtitle(paste(ordination, "-", distance))
  
  if (ordination == "PCA") {
    rot_df <- as.data.frame(ord_res$rotation)
    rot_df$Feature.ID <- rownames(rot_df)
    rot_df$mag <- sqrt(rot_df$PC1^2 + rot_df$PC2^2)
    rot_df <- rot_df[order(rot_df$mag, decreasing = TRUE), ][1:n_taxa, ]
    rot_df$PC1 <- rot_df$PC1 * 5
    rot_df$PC2 <- rot_df$PC2 * 5
    rot_df$Taxon <- taxonomy[match(rot_df$Feature.ID, feature_ids)]
    rot_df$label <- stringr::str_extract(rot_df$Taxon, "[^;]*$")
    
    p <- p +
      ggplot2::geom_segment(data = rot_df,
                            aes(x = 0, y = 0, xend = PC1, yend = PC2),
                            arrow = ggplot2::arrow(length = unit(0.3, "cm")),
                            color = "gray30",
                            inherit.aes = FALSE) +
      ggrepel::geom_label_repel(data = rot_df,
                                aes(x = PC1, y = PC2, label = label),
                                fill = "#EEEEEE", color = "black",
                                fontface = "italic", size = 4,
                                inherit.aes = FALSE)
  }
  
  return(p)
}
