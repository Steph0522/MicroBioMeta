#' Beta diversity partition (Jaccard/Sorensen) 
#'
#' This function computes beta diversity partition (Jaccard or Sorensen) 
#' using the betapart package, and plots the ordination (PCoA/NMDS) with gg_ordiplot.
#'
#' @param table Abundance matrix or data frame with taxa/features as rows and
#'   samples as columns (same orientation as the rest of the package). If a
#'   taxonomy column is present it is detected and removed automatically.
#' @param metadata Data frame with sample metadata. First column must be SampleID.
#' @param index Family of dissimilarity: "jaccard" (default) or "sorensen".
#' @param group_col Column in metadata to use as color grouping.
#' @param shape_col Optional column in metadata for point shapes.
#' @param legend_title Optional custom legend title.
#' @param group_colors Optional named vector of colors for groups.
#' @param point_size Numeric. Size of points in ordination plots. Default \code{3}.
#' @param panel_label_case Character. Case of the auto-generated A/B/C panel
#'   tags. One of \code{"upper"} (default, "A", "B", "C") or \code{"lower"}
#'   ("a", "b", "c"). Ignored if \code{panel_labels} is supplied.
#' @param panel_labels Optional character vector of 3 custom panel tags (one
#'   per jaccard/turnover/nestedness panel), used as-is (e.g.
#'   \code{c("(a)", "(b)", "(c)")} or \code{c("a.", "b.", "c.")}) — for
#'   journal styles that \code{panel_label_case} alone can't produce.
#'   Overrides \code{panel_label_case} when provided.
#' @param panel_label_bold Logical. If \code{TRUE} (default), panel tags are
#'   bold. Set to \code{FALSE} for journals that require plain (non-bold)
#'   panel tags.
#' @param save_table Logical. If \code{TRUE}, saves the dissimilarity table to disk. Default \code{FALSE}.
#' @param table_filename Character. Base name for the saved table file. Default \code{"SAMPLE1"}.
#'
#' @return A combined cowplot panel of beta diversity partition plots.
#' @export
#'
#' @examples
#' \dontrun{
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SampleID"
#'
#' beta_partition_plot(
#'   table      = table,
#'   metadata   = metadata,
#'   group_col  = "Location",
#'   point_size = 4
#' )
#' }


beta_partition_plot <- function(table, metadata, 
                                index = "jaccard", 
                                group_col = NULL, 
                                shape_col = NULL, 
                                legend_title = NULL,
                                point_size = 3,
                                group_colors = NULL,
                                panel_label_case = "upper",
                                panel_labels = NULL,
                                panel_label_bold = TRUE,
                                save_table = FALSE,
                                table_filename = "SAMPLE1") {

  if (!requireNamespace("ggordiplots", quietly = TRUE)) {
    stop(
      "Package 'ggordiplots' is required but not installed.\n",
      "Install it with: remotes::install_github(\"jfq3/ggordiplots\")",
      call. = FALSE
    )
  }

  suppressWarnings({

    # --- Aceptar también data.frame como matriz ---
    if (is.data.frame(table)) {
      
      # buscar columna de taxonomía
      tax_cols <- grep("taxonomy|taxon|Taxonomy|Taxa", names(table))
      
      if (length(tax_cols) > 0) {
        table <- table[, -tax_cols[1], drop = FALSE]
      }
      
      # conservar solo columnas numéricas
      table <- table[, sapply(table, is.numeric), drop = FALSE]
      
      # convertir a matriz
      table <- as.matrix(table)
    }
  
  # --- 1. Convert to presence/absence ---
  table_pa <- table
  table_pa[table_pa > 0] <- 1
  table_pa <- table_pa %>% 
    tibble::as_tibble(rownames = "SampleID") %>% 
    dplyr::arrange(SampleID) %>% 
    tibble::column_to_rownames(var = "SampleID") %>% 
    dplyr::select_if(is.numeric) %>% 
    t() %>% as.data.frame()
  
  colnames(metadata)[1] <- "SampleID"
  
  # --- 2. Beta diversity partition ---
  beta <- betapart::beta.pair(table_pa, index.family = index)
  if(index == "jaccard"){
    jac <- beta$beta.jac
    jtu <- beta$beta.jtu
    jne <- beta$beta.jne
  } else if(index == "sorensen"){
    jac <- beta$beta.sor
    jtu <- beta$beta.sim
    jne <- beta$beta.sne
  } else {
    stop("Only 'jaccard' or 'sorensen' are supported for index")
  }
  
  # --- 3. Merge with metadata ---
  
  #if (!"SampleID" %in% colnames(metadata)) {
   # stop("La tabla metadata no contiene una columna llamada 'SampleID'")
  #}
  
  names(metadata)[1] <- "SampleID"
  
  env1 <- table_pa %>% tibble::as_tibble(rownames = "SampleID") %>%
    dplyr::inner_join(metadata, by = "SampleID")
  
  # --- 4. betadisper for ordination ---
  jacs <- vegan::betadisper(jac, factor(env1[[group_col]]))
  jtus <- vegan::betadisper(jtu, factor(env1[[group_col]]))
  jnes <- vegan::betadisper(jne, factor(env1[[group_col]]))
  
  # Guardar tabla si se solicita
  if (save_table) {
    
    # convertir betadisper en data.frame
    betadisper_to_df <- function(bd_obj) {
      data.frame(
        SampleID = names(bd_obj$distances),
        Group = bd_obj$group,
        Distance_to_Centroid = bd_obj$distances,
        bd_obj$vectors,  # coordenadas PCoA
        check.names = FALSE
      )
    }
    
    # Crear los data frames
    jacs_df <- betadisper_to_df(jacs)
    jtus_df <- betadisper_to_df(jtus)
    jnes_df <- betadisper_to_df(jnes)
    
    # Obtener nombre base sin extensión
    base_name <- tools::file_path_sans_ext(table_filename)
    
    # Generar nombres únicos
    file_jacs <- paste0(base_name, "_jacs.txt")
    file_jtus <- paste0(base_name, "_jtus.txt")
    file_jnes <- paste0(base_name, "_jnes.txt")
    
    # Guardar cada tabla
    utils::write.table(jacs_df, file = file_jacs, sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(jtus_df, file = file_jtus, sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(jnes_df, file = file_jnes, sep = "\t", quote = FALSE, row.names = FALSE)
    
    message("Tables saved as:")
    message(file_jacs)
    message(file_jtus)
    message(file_jnes)
  }
  
  mean_jac <- round(mean(as.dist(jac)), 3)
  mean_turn <- round(mean(as.dist(jtu)), 3)
  mean_nes <- round(mean(as.dist(jne)), 3)
  
  # --- 5. Internal plotting function ---
  function_plot_beta <- function(x, env){
    y <- ggordiplots::gg_ordiplot(
      x, groups = env[[group_col]], hull = FALSE, 
      spiders = TRUE, ellipse = FALSE, plot = FALSE, label = TRUE
    )
    
    xlabs <- y$plot$labels$x
    ylabs <- y$plot$labels$y
    
    z <- ggplot2::ggplot() + 
      ggplot2::geom_point(
        data = y$df_ord %>% tibble::rownames_to_column(var = "SampleID") %>% 
          dplyr::inner_join(env, by = "SampleID"),
        ggplot2::aes(
          x = x, y = y,
          color = .data[[group_col]],
          shape = if (!is.null(shape_col)) .data[[shape_col]] else NULL
        ),
        size = point_size
      ) +
      ggplot2::xlab(xlabs) + ggplot2::ylab(ylabs) +
      ggplot2::geom_segment(
        data = y$df_spiders,
        ggplot2::aes(x = cntr.x, xend = x, y = cntr.y, yend = y, color = Group),
        show.legend = FALSE
      )
    
    color_scale <- if(!is.null(group_colors)) {
      ggplot2::scale_color_manual(values = group_colors)
    } else {
      ggplot2::scale_color_manual(values = .mbm_colors)
    }
    
    a <- z +
      ggplot2::geom_label(data = y$df_mean.ord, ggplot2::aes(x = x, y = y, label = Group),
                          fill = "white", color = "black", family = "serif",
                          fontface = "bold", size = 3.5) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      color_scale +
      ggplot2::labs(color = if(!is.null(legend_title)) legend_title else group_col) +
      .mbm_theme(
        legend_position = "right",
        extra = ggplot2::theme(
          legend.box        = "vertical",
          panel.grid.major  = ggplot2::element_blank(),
          plot.margin       = grid::unit(c(0, 0, 0, 0), "cm"),
          aspect.ratio      = 3/10,
          panel.border      = ggplot2::element_blank(),
          axis.line         = ggplot2::element_line(),
          axis.line.y.right = ggplot2::element_blank(),
          axis.line.x.top   = ggplot2::element_blank()
        )
      )
    
    
    return(a)  # <- faltaba este return + cierre
  }
  
  # --- 6. Generate plots ---
  plot_jac  <- function_plot_beta(jacs, env1) + 
    ggplot2::guides(
      colour = ggplot2::guide_legend(nrow = 1, title = if(!is.null(legend_title)) legend_title else group_col),
      shape  = ggplot2::guide_legend(nrow = 1)
    ) + ggplot2::theme(legend.position = "top")
  
  plot_turn <- function_plot_beta(jtus, env1)
  plot_nes  <- function_plot_beta(jnes, env1)
  
  leg <- cowplot::get_legend(plot_jac)
  resolved_labels <- if (!is.null(panel_labels)) {
    panel_labels
  } else if (identical(panel_label_case, "lower")) c("a", "b", "c") else c("A", "B", "C")
  panel_fontface <- if (panel_label_bold) "bold" else "plain"
  panel <- cowplot::plot_grid(
    plot_jac + ggplot2::theme(legend.position = "none") + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, color = "black", family = "serif", face = "bold")) +
      ggplot2::ylab("DIM2") + ggplot2::xlab("DIM1") + ggplot2::theme(aspect.ratio = 10/10) + ggplot2::ggtitle(paste0(index, " dissimilarity (mean = ", mean_jac, ")")),
     plot_turn + ggplot2::theme(legend.position = "none") + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, color = "black", family = "serif", face = "bold")) +
      ggplot2::ylab("") + ggplot2::xlab("DIM1") + ggplot2::theme(aspect.ratio = 10/10) + ggplot2::ggtitle(paste0("Turnover component (mean = ", mean_turn, ")")),
    plot_nes + ggplot2::theme(legend.position = "none") + ggplot2::theme(plot.title = ggplot2::element_text(size = 14, color = "black", family = "serif", face = "bold")) +
      ggplot2::ylab("") + ggplot2::xlab("DIM1") + ggplot2::theme(aspect.ratio = 10/10) + ggplot2::ggtitle(paste0("Nestedness component (mean = ", mean_nes, ")")),
    ncol = 3, align = "hv", labels = resolved_labels, label_fontfamily = "serif",
    label_fontface = panel_fontface
    )
  
  combined_plot <- cowplot::plot_grid(leg, panel, ncol = 1, rel_heights = c(0.1,1))
  
  return(combined_plot)
  })  # <- aquí se cierra suppressWarnings
}
