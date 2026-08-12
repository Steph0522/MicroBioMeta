#' Beta diversity plot with multiple distance and ordination methods
#'
#' This function computes beta diversity using several distance metrics and ordination methods
#' (PCA, PCoA, NMDS). It requires an abundance table with taxonomy, metadata, and allows customization
#' of color and shape aesthetics. It also supports compositional transformation via ALDEx2.
#'
#' @param table A data frame with abundances. The last column must contain taxonomy information.
#' @param metadata A data frame with sample metadata. The first column must contain the sample IDs.
#' @param distance Distance method: one of "euclidean", "bray", "jaccard", "sorensen",
#'  "compositional" (default), "aitchison", or "robust.aitchison".
#' @param ordination Ordination method: one of "PCA" (default), "PCoA", or "NMDS".
#' @param group_col Column in `metadata` to fill/color points. Its type
#'   decides the scale automatically: numeric columns (e.g. \code{"dist_km"})
#'   get a continuous scale; character/factor columns (e.g. \code{"estado2"})
#'   get a discrete qualitative scale.
#' @param palette Either a palette \strong{name} or a \strong{vector of fixed
#'   colors}; which scale it produces depends on whether \code{group_col} is
#'   discrete or continuous.
#'   \itemize{
#'     \item Named, discrete \code{group_col}: one of \code{"colorb"}
#'       (default; qualitative colorblind-friendly palette), \code{"grey"},
#'       \code{"viridis"}, or \code{"brewer"} (\code{"Set2"}).
#'     \item Named, continuous \code{group_col}: \code{"viridis"} (default;
#'       \code{option = "cividis"}, matching the urban-distance map figure)
#'       or \code{"gradient"} (colorblind-friendly blue-to-orange two-color
#'       gradient).
#'     \item Vector of colors, discrete \code{group_col}: used as-is, one
#'       color per level (\code{scale_*_manual}).
#'     \item Vector of colors, continuous \code{group_col}: used as gradient
#'       stops (\code{scale_*_gradientn}).
#'   }
#' @param shape_col Optional column in `metadata` to shape points.
#' @param legend_title Optional legend title.
#' @param top_n Number of top contributing taxa to display as arrows in PCA.
#' @param arrows_size Numeric. Size/length scaling factor for biplot arrows. Default \code{10}.
#' @param title Plot title. \code{"auto"} (default) generates \code{"Ordination - distance"};
#'   \code{NULL} shows no title; any other string is used as-is.
#' @param save_table Logical. If \code{TRUE}, saves a combined table of sample
#'   ordination scores and (when \code{ordination = "PCA"}) taxon loadings to
#'   disk, distinguished by a \code{type} column (\code{"site"} or
#'   \code{"loading"}). Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"ordination_scores.txt"}.
#'
#' @return A `ggplot2` object.
#' @export
#' @examples
#' \dontrun{
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SampleID"
#'
#' beta_div_plot(
#'   table      = table,
#'   metadata   = metadata,
#'   distance   = "aitchison",
#'   ordination = "NMDS",
#'   group_col  = "Type_of_soil",
#'   top_n      = 5
#' )
#' }


beta_div_plot <- function(table, metadata,
                          distance = "compositional",
                          ordination = "PCA",
                          group_col = NULL,
                          palette = "colorb",
                          shape_col = NULL,
                          legend_title = NULL,
                          arrows_size = 10,
                          top_n = 5,
                          title = "auto",
                          save_table = FALSE,
                          table_filename = "ordination_scores.txt") {
  
  requireNamespace("vegan")
  requireNamespace("ggplot2")
  requireNamespace("ggrepel")
  requireNamespace("ALDEx2")
  requireNamespace("dplyr")
  requireNamespace("stringr")
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  if (ordination == "PCA" && distance != "compositional") {
    stop("PCA is only available with 'compositional' (Aitchison). Use PCoA or NMDS instead.")
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
  otu_table <- as.data.frame(lapply(abund_table, as.numeric), check.names = FALSE)
  rownames(otu_table) <- feature_ids
  
  
  metadata_ids <- trimws(as.character(metadata[[1]]))
  sample_ids <- trimws(colnames(otu_table))
  colnames(otu_table) <- sample_ids
  metadata[[1]] <- metadata_ids
  
  common_samples <- intersect(sample_ids, metadata_ids)
  if (length(common_samples) == 0) stop("No matching samples between table and metadata.")
  otu_table <- otu_table[, common_samples, drop = FALSE]
  metadata <- metadata[metadata_ids %in% common_samples, ]
  
  if (distance == "compositional") {
    set.seed(123)
    aldex_obj <- ALDEx2::aldex.clr(otu_table, mc.samples = 128,
                                   denom = "all", verbose = FALSE, useMC = FALSE)
    otu_trans <- t(ALDEx2::getMonteCarloSample(aldex_obj, 1))
    dist_matrix <- dist(otu_trans, method = "euclidean")
  } else if (distance %in% c("aitchison", "robust.aitchison")) {
    dist_matrix <- vegan::vegdist(t(otu_table), method = distance, pseudocount = 0.5)
    otu_trans <- NULL
  } else {
    dist_matrix <- vegan::vegdist(t(otu_table), method = distance)
    otu_trans <- NULL
  }
  
  expl_var <- NULL
  res <- switch(ordination,
                "PCA" = {
                  pca_input <- if (!is.null(otu_trans)) otu_trans else t(otu_table)
                  pca <- prcomp(pca_input)
                  list(ord = pca, expl_var = round(100 * summary(pca)$importance[2, 1:2], 1))
                },
                "PCoA" = {
                  pcoa <- cmdscale(dist_matrix, eig = TRUE, k = 2)
                  eigs <- pcoa$eig
                  list(ord = pcoa, expl_var = round(100 * eigs[1:2] / sum(eigs[eigs > 0]), 1))
                },
                "NMDS" = list(ord = vegan::metaMDS(dist_matrix, k = 2, trymax = 100), expl_var = NULL),
                stop("Invalid ordination method.")
  )
  
  ord_res <- res$ord
  expl_var <- res$expl_var
  
  
  ord_df <- switch(ordination,
                   "PCA" = as.data.frame(ord_res$x),
                   "PCoA" = {
                     df <- as.data.frame(ord_res$points)
                     colnames(df) <- c("PCoA1", "PCoA2")  
                     df
                   },
                   "NMDS" = as.data.frame(ord_res$points))
  
  ord_df$SampleID <- rownames(ord_df)
  colnames(metadata)[1] <- "SampleID"
  merged <- dplyr::inner_join(ord_df, metadata, by = "SampleID")
  if (nrow(merged) == 0) stop("No common samples between table and metadata.")
  
  x_lab <- if (!is.null(expl_var)) paste0(names(ord_df)[1], " (", expl_var[1], "%)") else names(ord_df)[1]
  y_lab <- if (!is.null(expl_var)) paste0(names(ord_df)[2], " (", expl_var[2], "%)") else names(ord_df)[2]
  
  is_continuous <- !is.null(group_col) && is.numeric(merged[[group_col]])
  legend_name   <- ifelse(is.null(legend_title), group_col, legend_title)

  # `palette` can be a palette *name* (single string) or a *vector of fixed
  # colors*; either way, is_continuous decides whether it becomes a
  # gradient/gradientn scale or a discrete manual/qualitative one.
  is_named_palette <- is.character(palette) && length(palette) == 1

  .group_scale <- function(aesthetic) {
    fill <- aesthetic == "fill"
    if (is_continuous) {
      if (!is_named_palette) {
        if (fill) ggplot2::scale_fill_gradientn(name = legend_name, colours = palette)
        else ggplot2::scale_color_gradientn(name = legend_name, colours = palette)
      } else if (palette == "gradient") {
        if (fill) ggplot2::scale_fill_gradient(name = legend_name, low = "#0072B2", high = "#E69F00")
        else ggplot2::scale_color_gradient(name = legend_name, low = "#0072B2", high = "#E69F00")
      } else {
        # "viridis" (default) or any other name -> viridis cividis
        if (fill) ggplot2::scale_fill_viridis_c(name = legend_name, option = "cividis")
        else ggplot2::scale_color_viridis_c(name = legend_name, option = "cividis")
      }
    } else {
      if (!is_named_palette) {
        if (fill) ggplot2::scale_fill_manual(name = legend_name, values = palette)
        else ggplot2::scale_color_manual(name = legend_name, values = palette)
      } else {
        switch(palette,
          "grey"    = if (fill) ggplot2::scale_fill_grey(name = legend_name, start = 0.9, end = 0.3)
                      else ggplot2::scale_color_grey(name = legend_name, start = 0.9, end = 0.3),
          "viridis" = if (fill) ggplot2::scale_fill_viridis_d(name = legend_name)
                      else ggplot2::scale_color_viridis_d(name = legend_name),
          "brewer"  = if (fill) ggplot2::scale_fill_brewer(name = legend_name, palette = "Set2")
                      else ggplot2::scale_color_brewer(name = legend_name, palette = "Set2"),
          # "colorb" (default) or any other name -> package qualitative palette
          if (fill) ggplot2::scale_fill_manual(name = legend_name, values = .mbm_colors)
          else ggplot2::scale_color_manual(name = legend_name, values = .mbm_colors)
        )
      }
    }
  }

  if (is.null(shape_col)) {
    p <- ggplot2::ggplot(merged, ggplot2::aes(
      x = .data[[names(ord_df)[1]]],
      y = .data[[names(ord_df)[2]]],
      fill = .data[[group_col]]
    )) +
      ggplot2::geom_point(size = 4, shape = 21) +
      .group_scale("fill")

  } else {
    p <- ggplot2::ggplot(merged, ggplot2::aes(
      x = .data[[names(ord_df)[1]]],
      y = .data[[names(ord_df)[2]]],
      color = .data[[group_col]],
      shape = .data[[shape_col]]
    )) +
      ggplot2::geom_point(size = 4) +
      .group_scale("color")
  }

  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::labs(
      x     = x_lab,
      y     = y_lab,
      title = if (identical(title, "auto")) paste(ordination, "-", distance)
              else title   # NULL → no title; custom string → that text
    ) +
    .mbm_theme(
      legend_position = "right",
      extra = ggplot2::theme(
        legend.box        = "vertical",
        panel.grid.major  = ggplot2::element_blank()
      )
    )
  
  rot_df_out <- NULL
  if (ordination == "PCA") {
    rot_df <- as.data.frame(ord_res$rotation)
    rot_df$Feature.ID <- rownames(rot_df)
    rot_df$mag <- sqrt(rot_df$PC1^2 + rot_df$PC2^2)
    rot_df <- rot_df[order(rot_df$mag, decreasing = TRUE), ][1:top_n, ]
    rot_df$PC1 <- rot_df$PC1 * arrows_size
    rot_df$PC2 <- rot_df$PC2 * arrows_size
    rot_df$Taxon <- taxonomy[match(rot_df$Feature.ID, feature_ids)]
    #rot_df$label <- stringr::str_extract(rot_df$Taxon, "(?<=__)[^;]*$") 
    #rot_df$label <- gsub(" ", "\n", rot_df$label)
    
    ###
    extract_clean_label <- function(taxon_string, taxonomy_db = "silva") {
      
      # --- 1) Manejo de NA o vacio ---
      if (is.na(taxon_string) || taxon_string == "" || taxon_string == "Other") {
        return("Other")
      }
      
      # --- 2) Separar por niveles taxonomicos ---
      levels <- unlist(strsplit(taxon_string, ";"))
      levels <- trimws(levels)
      
      # --- 3) Limpieza por base de datos ---
      clean_by_db <- list(
        
        silva = function(x) sub("^[a-zA-Z]__", "", x),
        gg    = function(x) sub("^[a-zA-Z]__", "", x),
        unite = function(x) sub("^[a-zA-Z]__", "", x),
        Kraken2 = function(x) sub("^[a-zA-Z]__", "", x,) 
      )
      
      cleaner <- clean_by_db[[taxonomy_db]]
      
      # Si no existe el limpiador, usar limpieza generica
      if (is.null(cleaner)) {
        cleaner <- function(x) sub(".*__", "", x)
      }
      
      # --- 4) Limpiar niveles ---
      levels_clean <- vapply(levels, cleaner, FUN.VALUE = character(1))
      levels_clean <- trimws(levels_clean)
      
      
      # --- 5) Seleccion del nivel mas especifico valido ---
      invalid_literals <- c("", " ", "NA", "na", "unclassified", "Unassigned",
                            "uncultured", "uncultured_soil", "__")
      
      invalid_regex <- c("bacteriap[0-9]+")
      
      for (i in length(levels_clean):1) {
        
        lvl <- levels_clean[i]
        
        # ESTA ES LA LiNEA CORREGIDA
        if (!(lvl %in% invalid_literals) && 
            !any(grepl(invalid_regex, lvl, ignore.case = TRUE))) {
          
          # > Regla especial: Kraken2 species -> concatenar "Genus species"
          if (taxonomy_db == "Kraken2" && grepl("s__", levels[i])) {
            
            genus_full <- stringr::str_extract(taxon_string, "g__[^;]*")
            
            if (!is.na(genus_full)) {
              genus <- sub("g__", "", genus_full)
              species <- lvl
              
              # Evitar errores por empties
              if (genus != "" && species != "") {
                # opcional: reemplazar underscores
                species <- gsub("_", " ", species)
                return(paste(genus, species))
              }
            }
          }
          
          # Nivel normal
          return(lvl)
        }
      }
      
      return("Unclassified")
    }
    
    
    
    rot_df$label <- sapply(rot_df$Taxon, extract_clean_label)
    rot_df$label <- gsub(" ", "\n", rot_df$label)
    ###
    rot_df_out <- rot_df

    p <- p +
      ggplot2::geom_segment(data = rot_df,
                            ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
                            arrow = ggplot2::arrow(length = grid::unit(0.7, "cm")),
                            color = "gray30",
                            inherit.aes = FALSE) +
      ggrepel::geom_label_repel(data = rot_df,
                                ggplot2::aes(x = PC1, y = PC2, label = label),
                                fill = "white", color = "black",
                                fontface = "italic", size = 4, family= "serif",
                                inherit.aes = FALSE)
  }
  
  # --- Centrar ejes simetricamente ---
  x_limits <- range(merged[[names(ord_df)[1]]], na.rm = TRUE)
  y_limits <- range(merged[[names(ord_df)[2]]], na.rm = TRUE)
  max_range <- max(abs(x_limits), abs(y_limits))
  p <- p + ggplot2::coord_cartesian(xlim = c(-max_range, max_range),
                                    ylim = c(-max_range, max_range)) +
    ggplot2::coord_fixed()  # Mantiene proporcion 1:1

  if (save_table) {
    site_out <- merged
    colnames(site_out)[colnames(site_out) %in% names(ord_df)[1:2]] <- c("Axis1", "Axis2")
    site_out$type <- "site"
    if (!is.null(rot_df_out)) {
      loading_out <- rot_df_out
      colnames(loading_out)[colnames(loading_out) %in% c("PC1", "PC2")] <- c("Axis1", "Axis2")
      loading_out$type <- "loading"
      combined_table <- dplyr::bind_rows(site_out, loading_out)
    } else {
      combined_table <- site_out
    }
    utils::write.table(combined_table, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  return(p)
}
