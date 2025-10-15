#' Beta diversity partition (Jaccard/Sorensen) 
#'
#' This function computes beta diversity partition (Jaccard or Sorensen) 
#' using the betapart package, and plots the ordination (PCoA/NMDS) with gg_ordiplot.
#'
#' @param table Abundance matrix (samples in rows, species/features in columns).
#' @param metadata Data frame with sample metadata. First column must be SampleID.
#' @param index Family of dissimilarity: "jaccard" (default) or "sorensen".
#' @param group_col Column in metadata to use as color grouping.
#' @param shape_col Optional column in metadata for point shapes.
#' @param legend_title Optional custom legend title.
#' @param colors Optional named vector of colors for groups.
#'
#' @return A combined cowplot panel of beta diversity partition plots.
#' @export

beta_partition_plot <- function(table, metadata, 
                                index = "jaccard", 
                                group_col = NULL, 
                                shape_col = NULL, 
                                legend_title = NULL,
                                point_size = 3,
                                colors = NULL,
                                save_table = TRUE,
                                table_filename = "SAMPLE1") {   
  suppressWarnings({
    
  
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
        ggplot2::aes_string(
          x = "x", y = "y", 
          color = group_col, 
          shape = if(!is.null(shape_col)) shape_col else NULL
        ),
        size = point_size
      ) +
      ggplot2::xlab(xlabs) + ggplot2::ylab(ylabs) +
      ggplot2::geom_segment(
        data = y$df_spiders,
        ggplot2::aes(x = cntr.x, xend = x, y = cntr.y, yend = y, color = Group),
        show.legend = FALSE
      )
    
    color_scale <- if(!is.null(colors)) {
      ggplot2::scale_color_manual(values = colors)
    } else {
      ggplot2::scale_color_viridis_d(option = "turbo")
    }
    
    a <- z +
      ggplot2::geom_label(data = y$df_mean.ord, ggplot2::aes(x = x, y = y, label = Group)) +
      ggplot2::theme_linedraw() +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      color_scale +
      ggplot2::labs(color = if(!is.null(legend_title)) legend_title else group_col) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 12, color = "black", family = "Times New Roman"),
        axis.title = ggplot2::element_text(size = 14, color = "black", family = "Times New Roman"),
        legend.text = ggplot2::element_text(size = 12, color = "black", family = "Times New Roman"),
        legend.title = ggplot2::element_text(size = 14, color = "black", family = "Times New Roman", face = "bold"),
        legend.position = "right",
        legend.box = "vertical",
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = grid::unit(c(0, 0, 0, 0), "cm"),
        aspect.ratio = 3/10
      ) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(),
        axis.line.y.right = ggplot2::element_blank(),
        axis.line.x.top = ggplot2::element_blank()
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
  panel <- cowplot::plot_grid(
    plot_jac + ggplot2::theme(legend.position = "none") + ggplot2::theme(plot.title = element_text(size = 12, color = "black", family = "Times New Roman", face = "bold")) +
      ggplot2::ylab("DIM2") + ggplot2::xlab("DIM1") + ggplot2::theme(aspect.ratio = 10/10) + ggplot2::ggtitle(paste0(index, " dissimilarity (mean = ", mean_jac, ")")),
     plot_turn + ggplot2::theme(legend.position = "none") + ggplot2::theme(plot.title = element_text(size = 12, color = "black", family = "Times New Roman", face = "bold")) +
      ggplot2::ylab("") + ggplot2::xlab("DIM1") + ggplot2::theme(aspect.ratio = 10/10) + ggplot2::ggtitle(paste0("Turnover component (mean = ", mean_turn, ")")),
    plot_nes + ggplot2::theme(legend.position = "none") + ggplot2::theme(plot.title = element_text(size = 12, color = "black", family = "Times New Roman", face = "bold")) +
      ggplot2::ylab("") + ggplot2::xlab("DIM1") + ggplot2::theme(aspect.ratio = 10/10) + ggplot2::ggtitle(paste0("Nestedness component (mean = ", mean_nes, ")")),
    ncol = 3, align = "hv", labels = c("A", "B", "C"), label_fontfamily = "Times New Roman"
    )
  
  combined_plot <- cowplot::plot_grid(leg, panel, ncol = 1, rel_heights = c(0.1,1))
  
  return(combined_plot)
  })  # <- aquí se cierra suppressWarnings
}
