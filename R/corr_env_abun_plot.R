#' Plot Correlation Between Environmental Variables and Taxonomic Groups
#'
#' This function calculates the relative abundances of taxa at a specified taxonomic level 
#' (phylum, genus, or species) from a count table, computes correlations between these 
#' abundances and environmental variables, and visualizes the results as either a heatmap 
#' (tile) or a bubble plot (circle).
#'
#
#'
#' @param table A data frame containing a taxonomic abundance table with a
#'   taxonomy column and sample columns with numeric counts.
#' @param env_table A data frame or matrix of environmental variables, with
#'   samples as row names.
#' @param metadata A data frame containing sample metadata. The first column
#'   must correspond to sample identifiers.
#' @param cond_vect Character vector of environmental variables to include in
#'   the correlation analysis. If NULL, all available variables are used.
#' @param method Character. Correlation method passed to stats::cor and
#'   stats::cor.test. Supported options include ("spearman", "pearson", "kendall").
#' @param hc.order Logical. If TRUE, applies hierarchical clustering to
#'   reorder taxa and environmental variables in the plot.
#' @param geom Character. Type of visualization to generate: "tile"
#'   (heatmap) or "circle" (bubble plot).
#' @param show_labels Logical. If TRUE, displays correlation values on
#'   the plot.
#' @param col_palette Character vector defining the color palette for correlation
#'   values. If NULL, the palette is chosen via \code{diverging_palette}.
#' @param diverging_palette Character. Name of a built-in colorblind-friendly
#'   diverging palette to use when \code{col_palette} is NULL. One of
#'   \code{"BuOr"} (blue-orange, default), \code{"BuVm"} (blue-vermillion),
#'   \code{"BuPk"} (blue-pink), \code{"GnPk"} (green-pink).
#' @param invert_axes Logical. If TRUE, swaps x and y axes in the plot.
#' @param taxonomy_db Character. Taxonomic database used for annotation and
#'   parsing. Supported options include ("silva", "unite","Kraken2" and "gg2").
#' @param level Character. Taxonomic level to collapse taxa to. One of
#'   ("kingdom", "phylum", "class", "order", "family", "genus", "species").
#' @param pval_threshold Numeric. Optional p-value threshold to retain only taxa
#'   showing significant correlations with at least one environmental variable.
#'   If \code{NULL}, no significance filtering is applied.
#' @param save_table Logical. If TRUE, saves the correlation matrix as a
#'   tab-delimited text file.
#' @param table_filename Character. Name of the output file used when
#'   save_table = TRUE.
#'
#' @return A ggplot2 object.
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
#' # env_table must have rownames matching the sample names in `table`
#' env_data <- metadata
#' rownames(env_data) <- env_data$SampleID
#'
#' corr_env_abund_plot(
#'   table          = table,
#'   env_table      = env_data,
#'   metadata       = metadata,
#'   cond_vect      = c("pH", "TOC", "FW", "Root_FW", "DW", "Root_L", "Stem_L"),
#'   method         = "pearson",
#'   geom           = "tile",
#'   hc.order       = FALSE,
#'   col_palette    = c("pink", "white", "purple"),
#'   invert_axes    = TRUE,
#'   show_labels    = FALSE,
#'   level          = "phylum",
#'   taxonomy_db    = "silva",
#'   pval_threshold = 0.05
#' )
#' }

corr_env_abund_plot <- function(table,
                                env_table,
                                metadata = NULL,
                                cond_vect = NULL,
                                method = "spearman",
                                hc.order = TRUE,
                                geom = c("tile", "circle"),
                                show_labels = TRUE,
                                col_palette = NULL,
                                diverging_palette = "BuOr",
                                invert_axes = TRUE,
                                taxonomy_db = "silva",
                                level = "genus",
                                pval_threshold = NULL,
                                save_table = FALSE,
                                table_filename = "corr.txt") {
  geom <- match.arg(geom)
  rownames(table) <- NULL
  
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  # Remove uninformative taxonomy strings
  table <- table %>%
    dplyr::filter(taxonomy != "d__Bacteria;__;__;__;__;__") %>%
    dplyr::filter(taxonomy != "d__Bacteria") %>%
    dplyr::filter(taxonomy != "d__Archaea;__;__;__;__;__") %>%
    dplyr::filter(taxonomy != "d__Archaea") %>%
    dplyr::filter(taxonomy != "d__Bacteria;p__;c__;o__;f__;g__;s__") %>%
    dplyr::filter(taxonomy != "d__Archaea;p__;c__;o__;f__;g__;s__") %>%
    dplyr::filter(taxonomy != "k__Bacteria;__;__;__;__;__")%>%
    dplyr::filter(taxonomy != "k__Fungi;__;__;__;__;__")%>%
    dplyr::filter(taxonomy != "k__Fungi;p__;c__;o__;f__;g__")%>%
    dplyr::filter(taxonomy != "k__Fungi")%>%
    dplyr::filter(taxonomy != "Unassigned")%>%
    dplyr::filter(taxonomy != "d__Eukaryota")
  
  
  # Align samples across table, env_table, and (optionally) metadata
  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata)
    common_samples <- Reduce(intersect, list(colnames(table), rownames(env_table), metadata[[1]]))
    table <- table[, c(tax_col, match(common_samples, colnames(table))), drop = FALSE]
    env_table <- env_table[common_samples, , drop = FALSE]
    metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
    rownames(metadata) <- metadata[[1]]
    table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
    ordered_samples <- intersect(metadata[, 1], colnames(table)[-1])
    table <- table[, c("taxonomy", ordered_samples)]
  } else {
    common_samples <- intersect(colnames(table), rownames(env_table))
    table <- table[, c(tax_col, match(common_samples, colnames(table))), drop = FALSE]
    env_table <- env_table[common_samples, , drop = FALSE]
    table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
    ordered_samples <- common_samples
    table <- table[, c("taxonomy", ordered_samples)]
  }
  
  # --- Asegurar identificadores únicos ---
  if (!is.null(rownames(table))) {
    table <- tibble::rownames_to_column(table, var = "OTU_ID")
  } else {
    table$OTU_ID <- paste0("OTU_", seq_len(nrow(table)))
  }
  
  # --- Limpiar terminaciones vacías en la taxonomía ---
  table$taxonomy <- gsub("(;__)+$", "", table$taxonomy)
  
  # --- Determinar índice del nivel taxonómico ---
  level_idx <- switch(level,
                      kingdom = 1, phylum = 2, class = 3, order = 4,
                      family = 5, genus = 6, species = 7)
  
  # --- Calcular profundidad taxonómica de cada OTU ---
  get_depth <- function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    sum(grepl("__", levels))
  }
  table$depth <- sapply(table$taxonomy, get_depth)
  
  # --- Separar filas según resolución taxonómica ---
  lowres <- table[table$depth < level_idx, ]    
  highres <- table[table$depth >= level_idx, ]  
  
  # --- Recortar y colapsar taxonomías con suficiente resolución ---
  highres$taxonomy <- sapply(highres$taxonomy, function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    paste(levels[1:level_idx], collapse = ";")
  })
  
  # --- Colapsar correctamente taxones repetidos ---
  highres <- highres %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(
      OTU_ID = paste(unique(OTU_ID), collapse = ";"),
      dplyr::across(dplyr::where(is.numeric), \(x) sum(x, na.rm = TRUE)),
      .groups = "drop"
    )
  
  
  # --- Combinar lowres y highres ---
  table_final <- dplyr::bind_rows(
    lowres[, c("OTU_ID", "taxonomy", ordered_samples)],
    highres[, c("OTU_ID", "taxonomy", ordered_samples)]
  )
  
  # --- Limpiar tabla final ---
  table_final <- table_final[, c("OTU_ID",  ordered_samples, "taxonomy")]
  table_final <- tibble::column_to_rownames(table_final, "OTU_ID")
  
  
  # rownames(table) <- make.unique(table$taxonomy)
  #table$taxonomy <- NULL
  
  table <- table_final
  
  
  #modificar la columna taxonomy para solo conservar el nombre al nivel que colapsamos
  #esto hace que al graficar salga sólo ese nombre y no toda la taxonomía
  
  
  if (taxonomy_db %in% c("unite","silva", "gg2") && level == "species") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*;.*s__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$|s__uncultured|s__$", taxonomy) ~
            paste0(
              stringr::str_extract(taxonomy, "s__[^;]*") %>% sub("s__", "", .)
            ),
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .)),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  if (taxonomy_db == "Kraken2" && level == "species") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*;.*s__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$|s__uncultured|s__$", taxonomy) ~
            paste0(
              stringr::str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .), " ",
              stringr::str_extract(taxonomy, "s__[^;]*") %>% sub("s__", "", .)
            ),
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "g__[^;]*") %>% sub("g__", "", .)),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~
            paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  
  # Simplify taxonomy for SILVA
  if (taxonomy_db %in% c("silva") && level == "genus") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$|g__Incertae_Sedis", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$|f__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$|o__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  
  if (taxonomy_db %in% c("unite", "Kraken2", "gg2") && level == "genus") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  if (taxonomy_db %in% c("unite","silva", "Kraken2","gg2") && level == "family") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$|f__Incertae_Sedis", taxonomy) ~ sub(".*f__([^;]*).*", "\\1", taxonomy),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$|o__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  if (taxonomy_db %in% c("unite","silva", "Kraken2","gg2") && level == "order") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$|o__Incertae_Sedis", taxonomy) ~ sub(".*o__([^;]*).*", "\\1", taxonomy),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  if (taxonomy_db %in% c("unite","silva", "Kraken2","gg2") && level == "class") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ sub(".*c__([^;]*).*", "\\1", taxonomy),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  if (taxonomy_db %in% c("unite","silva", "Kraken2","gg2") && level == "phylum") {
    table <- table %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("p__[^;]*", taxonomy) ~ stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  # 1️⃣ Convertir nombres de fila en columna OTU_ID
  table <- table %>%
    tibble::rownames_to_column(var = "OTU_ID")

  # 2️⃣ Crear una nueva tabla llamada taxon con OTU_ID y taxonomy
  taxon <- table %>%
    dplyr::select(OTU_ID, taxonomy)

  # 3️⃣ Eliminar la columna taxonomy de table
  table <- table %>%
    dplyr::select(-taxonomy)

  # 4️⃣ Convertir OTU_ID nuevamente en nombres de fila
  table <- table %>%
    tibble::column_to_rownames(var = "OTU_ID")
  
  
  
  # Paleta por defecto
  if (is.null(col_palette)) {
    preset <- .mbm_div_palettes[[diverging_palette]]
    if (is.null(preset)) {
      warning("Unknown diverging_palette '", diverging_palette,
              "'. Using 'BuOr'. Valid options: ",
              paste(names(.mbm_div_palettes), collapse = ", "))
      preset <- .mbm_div_palettes[["BuOr"]]
    }
    col_palette <- grDevices::colorRampPalette(preset)(200)
  }
  # --- Filas comunes
  common_samples <- base::intersect(colnames(table), rownames(env_table))
  counts <- table[, common_samples, drop = FALSE]
  env <- env_table[common_samples, , drop = FALSE]
  
  
  # Seleccionar solo las variables ambientales indicadas en cond_vect, verificando coincidencias
  if (!is.null(cond_vect)) {
    cond_vect <- cond_vect[cond_vect %in% colnames(env)]
    if(length(cond_vect) == 0) stop("No matching variables found in env_table")
    env <- env[, cond_vect, drop = FALSE]
  }
  
  
  # Descartar columnas no numericas (p.ej. IDs o variables categoricas
  # mezcladas en env_table); stats::cor() requiere que env sea todo numerico.
  is_num <- vapply(env, is.numeric, logical(1))
  if (!all(is_num)) {
    warning("Dropping non-numeric columns from `env_table`: ",
            paste(names(env)[!is_num], collapse = ", "))
    env <- env[, is_num, drop = FALSE]
  }

  # Filtrar variables constantes
  env <- env[, apply(env, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  abund <- sweep(counts, 2, colSums(counts, na.rm = TRUE), FUN = "/") * 100
  abund <- abund[apply(abund, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
  
  # Matriz de correlación general
  corr_mat <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
  
  
  
  # Guardar tabla si se solicita
  if (save_table) {
    # Convertir a data frame y conservar nombres de filas
    corr_df <- as.data.frame(corr_mat)
    corr_df <- cbind(Taxon = rownames(corr_df), corr_df)
    
    utils::write.table(
      corr_df,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    message(paste("Table saved as:", table_filename))
  }
  
  # --- Calcular p-values si se indica pval_threshold
  if (!is.null(pval_threshold)) {
    pval_mat <- matrix(NA, 
                       nrow = ncol(env), 
                       ncol = nrow(abund),
                       dimnames = list(colnames(env), rownames(abund)))
    
    for (env_var in colnames(env)) {
      for (tax_id in rownames(abund)) {
        test <- suppressWarnings(
          cor.test(env[[env_var]], as.numeric(abund[tax_id, ]), method = method)
        )
        pval_mat[env_var, tax_id] <- test$p.value
      }
    }
    
    # Mantener solo taxones significativos en al menos una variable
    signif_taxa <- rownames(abund)[apply(pval_mat, 2, function(x) any(x < pval_threshold, na.rm = TRUE))]
    abund <- abund[signif_taxa, , drop = FALSE]
    corr_mat <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
  }
  
  
  
  # Clustering jerárquico
  if (hc.order) {
    d_row <- stats::dist(1 - corr_mat)
    hc_row <- stats::hclust(d_row)
    row_ord <- hc_row$labels[hc_row$order]
    d_col <- stats::dist(1 - t(corr_mat))
    hc_col <- stats::hclust(d_col)
    col_ord <- hc_col$labels[hc_col$order]
    corr_mat <- corr_mat[row_ord, col_ord]
  }
  
  # Convertir a formato largo
  corr_df <- reshape2::melt(corr_mat,
                            varnames = c("Environmental", "Group"),
                            value.name = "Correlation")
  # --- Añadir taxonomía a los resultados derretidos ---
  # corr_df$Group son los OTU_ID, vamos a unir con la tabla taxon
  corr_df <- dplyr::left_join(
    corr_df,
    taxon,
    by = c("Group" = "OTU_ID")
  )
  
  # Reemplazar "Group" por el nombre taxonómico
  corr_df$Taxon <- corr_df$taxonomy
  corr_df$taxonomy <- NULL
  # Base del gráfico
  if (invert_axes) {
    p <- ggplot2::ggplot(corr_df,
                         ggplot2::aes(x = Environmental, y = Taxon, fill = Correlation))
  } else {
    p <- ggplot2::ggplot(corr_df,
                         ggplot2::aes(x = Taxon, y = Environmental, fill = Correlation))
  }
  
  # Elegir tipo de gráfico, "tile" es como heatmap y "circle" como bubbleplot
  if (geom == "tile") {
    p <- p + ggplot2::geom_tile(color = "gray80")
    if (show_labels) {
      p <-
        p + ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)),
                               size = 3,
                               color = "black")
    }
  } else if (geom == "circle") {
    p <- p + ggplot2::geom_point(ggplot2::aes(size = abs(Correlation)),
                                 shape = 21,
                                 color = "gray") +
      ggplot2::scale_size(range = c(2, 10))
    if (show_labels) {
      p <-
        p + ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 2)),
                               size = 3,
                               vjust = 0.5)
    }
  }
  
  # Colores y tema
  p <- p +
    ggplot2::scale_fill_gradientn(colours = col_palette,
                                  limits = c(-1, 1),
                                  name = "Correlation") +
    .mbm_theme(
      legend_position = "right",
      extra = ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, vjust = 1,
                                            hjust = 1, size = 12,
                                            color = "black"),
        panel.grid.major = ggplot2::element_blank()
      )
    ) +
    ggplot2::coord_fixed()
  
  return(p)
}
