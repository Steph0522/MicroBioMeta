#' CCA/RDA Biplot with ggplot2
#'
#' Performs Canonical Correspondence Analysis (CCA) or Redundancy Analysis (RDA) based on a species abundance table
#' and selected environmental variables, returning a biplot with ggplot2 that visualizes
#' sample scores and environmental vectors.
#'
#' @param table A data frame or matrix of species abundances (samples as rows, species as columns).
#' @param env_data A data frame of environmental variables (rows must match `table`).
#' @param env_vars A character vector with the names of environmental variables to include in the analysis.
#' @param method Transformation method passed to `decostand` (default is `"hell"` for Hellinger).
#' @param metadata Optional data frame with sample metadata for grouping in the plot.
#' @param group_col Optional name of the column in `metadata` used to define sample groups.
#' @param group_colors Optional named vector of colors to use for each group.
#' @param legend_title Optional custom title for the group legend.
#' @param scale_env Logical; whether to scale environmental variables (default is `TRUE`).
#' @param pval_threshold P-value threshold for selecting significant environmental variables (default is `0.05`).
#' @param show_all_env_vectors Logical; if TRUE, plot all environmental vectors regardless of significance.
#' @param analysis Either `"CCA"` or `"RDA"` (default is `"CCA"`).
#' @param seed Random seed for reproducibility (default is `126`).
#' @param scale_arrows Numeric value to scale environmental vectors in the plot.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object displaying the biplot with sample scores and environmental vectors.
#' @export

cca_rda_biplot <- function(table,
                       env_data,
                       env_vars,
                       method = "hell",
                       metadata = NULL,
                       group_col = NULL,
                       group_colors = NULL,
                       legend_title = NULL,
                       scale_env = TRUE,
                       pval_threshold = 0.05,
                       show_all_env_vectors = FALSE,
                       analysis = "CCA",
                       seed = 126,
                       scale_arrows = 1,
                       title = NULL) {
  require(vegan)
  require(ggplot2)
  require(dplyr)
  
  # 1. Process species table (remove last column = taxonomy)
  taxonomy_col <- ncol(table)
  spp_table <- table[, -taxonomy_col, drop = FALSE]
  
  # Transpose table
    spp_table <- t(spp_table)
  
  
  # 2. Process metadata and handle hash IDs
  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata)
    
    # Store original IDs for error messages
    original_table_ids <- rownames(spp_table)
    original_meta_ids <- if(all(rownames(metadata) == as.character(seq_len(nrow(metadata))))) {
      metadata[[1]]
    } else {
      rownames(metadata)
    }
    
    # Case 1: Exact matching possible
    common_samples <- intersect(rownames(spp_table), original_meta_ids)
    
    # Case 2: No exact matches but same number of samples -> match by position
    if (length(common_samples) == 0 && nrow(spp_table) == length(original_meta_ids)) {
      warning("No exact ID matches found. Matching samples by position.", call. = FALSE)
      rownames(spp_table) <- original_meta_ids
      common_samples <- original_meta_ids
    }
    
    if (length(common_samples) == 0) {
      stop("No matching samples found.\n",
           "Table samples (first 6): ", paste(head(original_table_ids), collapse = ", "), "\n",
           "Metadata samples (first 6): ", paste(head(original_meta_ids), collapse = ", "), "\n\n",
           "Solutions:\n",
           "1. Ensure metadata has a column with matching sample IDs\n",
           "2. Provide metadata in the same order as the table\n")
    }
    
    # Apply the matching
    spp_table <- spp_table[common_samples, , drop = FALSE]
    metadata <- metadata[match(common_samples, original_meta_ids), , drop = FALSE]
    rownames(metadata) <- common_samples
  }
  

  # 4. Transform species data
  spp_hell <- vegan::decostand(spp_table, method = method)
  
  # 2. Scale env data
  if (scale_env) {
    env_scaled <- scale(env_data[, env_vars], scale = TRUE, center = FALSE) %>% as.data.frame()
  } else {
    env_scaled <- env_data[, env_vars]
  }
  
  # 3. Verificar correspondencia de filas
  stopifnot(identical(rownames(spp_table), rownames(env_data)))
  
  # 4. Ejecutar CCA o RDA según análisis
  set.seed(seed)
  if (toupper(analysis) == "RDA") {
    ord_result <- vegan::rda(spp_hell ~ ., data = env_scaled)
    axis_names <- c("RDA1", "RDA2")
  } else {
    ord_result <- vegan::cca(spp_hell ~ ., data = env_scaled)
    axis_names <- c("CCA1", "CCA2")
  }
  
  # 5. Ajuste de vectores ambientales
  fit <- vegan::envfit(ord_result, env_scaled)
  
  # 6. Selección de variables a graficar
  if (show_all_env_vectors) {
    vars_to_plot <- rownames(vegan::scores(fit, display = "vectors"))
  } else {
    sig_vars <- names(which(fit$vectors$pvals < pval_threshold))
    if (length(sig_vars) == 0) {
      warning("No hay variables ambientales significativas (p <", pval_threshold, ")")
      return(ggplot() + theme_void() + ggtitle("Sin variables significativas"))
    }
    vars_to_plot <- sig_vars
  }
  
  vectors_scores <- vegan::scores(fit, display = "vectors")[vars_to_plot, , drop = FALSE] %>%
    as.data.frame()
  vectors_scores$Variable <- rownames(vectors_scores)
  
  # 7. Coordenadas de sitios
  site_scores <- vegan::scores(ord_result, display = "sites") %>% as.data.frame()
  colnames(site_scores)[1:2] <- axis_names
  site_scores$SampleID <- rownames(site_scores)
  
  # 8. Agregar metadata
  if (!is.null(metadata) && !is.null(group_col) && group_col %in% colnames(metadata)) {
    if (!"SampleID" %in% colnames(metadata)) {
      metadata$SampleID <- rownames(metadata)
    }
    site_scores <- merge(site_scores, metadata[, c("SampleID", group_col)], by = "SampleID", all.x = TRUE)
    colnames(site_scores)[colnames(site_scores) == group_col] <- "Group"
    
    default_colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
                        "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
    groups_present <- unique(site_scores$Group)
    
    if (is.null(group_colors)) {
      color_values <- rep(default_colors, length.out = length(groups_present))
      names(color_values) <- groups_present
      group_colors <- color_values
    }
    
    legend_name <- ifelse(is.null(legend_title), group_col, legend_title)
    
    plot <- ggplot(site_scores, aes_string(x = axis_names[1], y = axis_names[2], fill = "Group")) +
      geom_point(size = 4, shape=21 ) +
      scale_fill_manual(name = legend_name, values = group_colors)
  } else {
    plot <- ggplot(site_scores, aes_string(x = axis_names[1], y = axis_names[2])) +
      geom_point(size = 4, shape=21)
  }
  
  # 9. Añadir vectores ambientales
  plot <- plot +
    geom_segment(
      data = vectors_scores,
      aes_string(x = 0, y = 0,
                 xend = paste0(axis_names[1], " * scale_arrows"),
                 yend = paste0(axis_names[2], " * scale_arrows")),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "black",
      inherit.aes = FALSE
    ) +
    geom_text(
      data = vectors_scores,
      aes_string(
        x = paste0(axis_names[1], " * scale_arrows"),
        y = paste0(axis_names[2], " * scale_arrows"),
        label = "Variable"
      ),
      color = "black",
      hjust = 0.5,
      vjust = -0.5,
      inherit.aes = FALSE
    ) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5)
    ) +
    geom_hline(yintercept = 0, color = "black") +
    geom_vline(xintercept = 0, color = "black")
  
  # 10. Escalar límites del gráfico
  max_range <- max(abs(c(site_scores[[axis_names[1]]],
                         vectors_scores[[axis_names[1]]] * scale_arrows,
                         site_scores[[axis_names[2]]],
                         vectors_scores[[axis_names[2]]] * scale_arrows)))
  buffer <- 1.1
  plot <- plot +
    scale_x_continuous(limits = c(-max_range, max_range) * buffer) +
    scale_y_continuous(limits = c(-max_range, max_range) * buffer)+
    theme_linedraw() +
    theme(axis.text = element_text(colour = "black", size = 8),
          axis.title = element_text(colour = "black", size = 12),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.position = "right",
          legend.box = "vertical",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    geom_vline(xintercept = 0, linetype = 2, color="#88929b") +   #lines-cross
    geom_hline(yintercept = 0, linetype = 2, color="#88929b") +
    
    guides(
      fill=guide_legend(title=legend_title))+#ylab("CAP2")+xlab("CAP1")+
    theme(axis.text = element_text(size = 12))
  
  # 11. Título del gráfico
  if (!is.null(title)) {
    plot <- plot + ggtitle(title)
  }
  
  return(plot)
}
