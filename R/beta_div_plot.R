#' Beta diversity plot with multiple distance and ordination methods
#'
#' This function computes beta diversity using several distance metrics and ordination methods
#' (PCA, PCoA, NMDS). It requires an abundance table with taxonomy, metadata, and allows customization
#' of color and shape aesthetics. It also supports compositional transformation via ALDEx2.
#'
#' @param table A data frame with abundances. The last column must contain taxonomy information, regardless of its name.
#' @param metadata A data frame with sample metadata. The first column must contain the sample IDs, regardless of the column name.
#' @param distance Distance method: one of "euclidean", "bray", "jaccard", "sorensen", or "aitchison" (default).
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

beta_div_plot <- function(table, metadata, distance = "aitchison",
                          ordination = "PCA", color_by, shape_by = NULL,
                          n_taxa = 5) {
  
  requireNamespace("vegan")
  requireNamespace("ggplot2")
  requireNamespace("ggrepel")
  requireNamespace("ALDEx2")
  requireNamespace("dplyr")
  requireNamespace("stringr")
  
  # Detect taxonomy column (assume it's the last one)
  tax_col <- names(table)[ncol(table)]
  taxonomy <- table[[tax_col]]
  abund_table <- table[, -ncol(table)]
  
  # Extract feature IDs
  if ("Feature.ID" %in% names(abund_table)) {
    feature_ids <- abund_table$Feature.ID
    abund_table <- abund_table[, !(names(abund_table) == "Feature.ID")]
  } else {
    feature_ids <- rownames(abund_table)
  }
  rownames(abund_table) <- feature_ids
  
  otu_table <- as.data.frame(lapply(abund_table, as.numeric))
  
  # Aitchison transformation
  if (distance == "aitchison") {
    set.seed(123)
    aldex_obj <- ALDEx2::aldex.clr(t(otu_table), mc.samples = 128,
                                   denom = "all", verbose = FALSE, useMC = FALSE)
    otu_trans <- t(ALDEx2::getMonteCarloSample(aldex_obj, 1))
    dist_matrix <- dist(otu_trans, method = "euclidean")
  } else {
    dist_matrix <- vegan::vegdist(t(otu_table), method = distance)
    otu_trans <- NULL
  }
  
  # Ordination
  ord_res <- switch(ordination,
                    "PCA" = prcomp(otu_trans),
                    "PCoA" = cmdscale(dist_matrix, eig = TRUE, k = 2),
                    "NMDS" = vegan::metaMDS(dist_matrix, k = 2, trymax = 100),
                    stop("Invalid ordination method."))
  
  ord_df <- switch(ordination,
                   "PCA" = as.data.frame(ord_res$x),
                   "PCoA" = as.data.frame(ord_res$points),
                   "NMDS" = as.data.frame(ord_res$points))
  ord_df$SampleID <- rownames(ord_df)
  
  # Metadata merge
  colnames(metadata)[1] <- "SampleID"
  merged <- dplyr::inner_join(ord_df, metadata, by = "SampleID")
  
  # Plot
  p <- ggplot2::ggplot(merged, aes_string(x = names(ord_df)[1],
                                          y = names(ord_df)[2],
                                          color = color_by,
                                          shape = shape_by)) +
    ggplot2::geom_point(size = 4) +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle(paste(ordination, "-", distance))
  
  # PCA taxon arrows
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
                            inherit.aes = FALSE) +
      ggrepel::geom_label_repel(data = rot_df,
                                aes(x = PC1, y = PC2, label = label),
                                fill = "#EEEEEE", color = "black",
                                fontface = "italic", size = 4)
  }
  
  return(p)
}