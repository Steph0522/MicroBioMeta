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
#' @param env_data A data frame or matrix of environmental variables, with
#'   samples as row names.
#' @param metadata A data frame containing sample metadata. The first column
#'   must correspond to sample identifiers.
#' @param env_vars Character vector of environmental variables to include in
#'   the correlation analysis. If NULL, all available variables are used.
#' @param method Character. Correlation method passed to stats::cor and
#'   stats::cor.test. Supported options include ("spearman", "pearson", "kendall").
#' @param hc.order Logical. If TRUE, applies hierarchical clustering to
#'   reorder taxa and environmental variables in the plot.
#' @param geom Character. Type of visualization to generate: `"tile"`
#'   (heatmap, default) or `"circle"` (bubble plot). Case-insensitive.
#' @param show_labels Logical. If TRUE, displays correlation values on
#'   the plot.
#' @param col_palette Character vector defining the color palette for correlation
#'   values. If NULL, the palette is chosen via \code{diverging_palette}.
#' @param diverging_palette Character. Name of a built-in colorblind-friendly
#'   diverging palette to use when \code{col_palette} is NULL. One of
#'   \code{"BuOr"} (default; blue-white-orange, the same Okabe-Ito blue/
#'   orange pairing used for the two-group color convention elsewhere in the
#'   package, e.g. \code{aldex_volcano_plot}'s col_inf/col_sup and the
#'   Rhizosphere/Roots colors in the bundled examples), \code{"BuVm"}
#'   (blue-vermillion), \code{"BuPk"} (blue-pink), \code{"GnPk"}
#'   (green-pink), \code{"PuYl"} (purple-white-yellow, viridis endpoints);
#'   or \code{"viridis"} for the plain sequential purple-to-yellow scale
#'   used in \code{aldex_heatmap_plot} (no neutral midpoint - not
#'   recommended for correlations, where 0 should look distinct from either
#'   extreme). All presets have a true white midpoint at 0 except
#'   \code{"viridis"}.
#' @param invert_axes Logical. If TRUE, swaps x and y axes in the plot.
#' @param taxonomy_db Character. Taxonomic database whose prefix style is used
#'   for parsing/annotation. One of `"silva"` (default), `"gg2"` (also accepts
#'   `"gg"` / `"greengenes2"`), `"unite"`, or `"Kraken2"` (also accepts
#'   `"kraken"`). Case-insensitive.
#' @param level Character. Taxonomic level to collapse taxa to. One of
#'   `"kingdom"`, `"phylum"`, `"class"`, `"order"`, `"family"`, `"genus"`
#'   (default), or `"species"`. Case-insensitive.
#' @param pval_threshold Numeric. Optional p-value threshold to retain only taxa
#'   showing significant correlations with at least one environmental variable.
#'   If \code{NULL}, no significance filtering is applied.
#' @param p_adjust_method Multiple-comparison correction applied to the
#'   p-values of all taxon x variable correlations before filtering with
#'   \code{pval_threshold}; any method of \code{stats::p.adjust()}. Default
#'   \code{"BH"} (false discovery rate); \code{"none"} uses the raw p-values.
#' @param x_label_angle Numeric. Rotation (in degrees) of the x-axis labels,
#'   so long variable or taxon names don't overlap. Default \code{45};
#'   \code{0} for horizontal labels.
#' @param save_table Logical. If TRUE, saves the plotted correlations (one row
#'   per taxon and variable, with the taxon name, raw p-value and p-value
#'   adjusted with \code{p_adjust_method}) as a
#'   tab-delimited text file.
#' @param table_filename Character. Name of the output file used when
#'   save_table = TRUE.
#'
#' @param ... Old names of renamed arguments (\code{env_table}, \code{cond_vect}), still accepted
#'   with a warning. Any other extra argument is an error.
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#' colnames(metadata)[1] <- "SampleID"
#'
#' # env_data must have rownames matching the sample names in `table`
#' env_data <- metadata
#' rownames(env_data) <- env_data$SampleID
#'
#' # Uses the default colorblind-friendly diverging palette (purple-white-yellow)
#' corr_env_abund_plot(
#'   table          = table,
#'   env_data      = env_data,
#'   metadata       = metadata,
#'   env_vars      = c("pH", "TOC", "FW", "Root_FW", "DW", "Root_L", "Stem_L"),
#'   method         = "pearson",
#'   geom           = "tile",
#'   hc.order       = FALSE,
#'   invert_axes    = TRUE,
#'   show_labels    = FALSE,
#'   level          = "phylum",
#'   taxonomy_db    = "silva"
#' )
#' # pval_threshold = 0.05 would keep only taxa with a significant correlation
#' # after the p_adjust_method correction (BH by default).

corr_env_abund_plot <- function(table,
                                env_data,
                                metadata = NULL,
                                env_vars = NULL,
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
                                p_adjust_method = "BH",
                                x_label_angle = 45,
                                save_table = FALSE,
                                table_filename = "corr.txt",
                                ...) {
  # Old argument names still work, with a warning (see .mbm_renamed_args)
  renamed <- .mbm_renamed_args(list(...), c(env_table = "env_data", cond_vect = "env_vars"), "corr_env_abund_plot")
  for (nm in names(renamed)) assign(nm, renamed[[nm]])

  geom <- match.arg(tolower(geom), c("tile", "circle"))

  # Accept taxonomy_db / level case-insensitively (mapping to the exact value
  # the rest of the function expects), so users don't have to remember the
  # internal capitalisation (e.g. Kraken2).
  taxonomy_db <- switch(
    tolower(taxonomy_db),
    "silva"       = "silva",
    "gg"          = ,
    "gg2"         = ,
    "greengenes2" = "gg2",
    "unite"       = "unite",
    "kraken"      = ,
    "kraken2"     = "Kraken2",
    stop("Invalid `taxonomy_db`: '", taxonomy_db,
         "'. Choose one of: \"silva\", \"gg2\", \"unite\", \"Kraken2\" ",
         "(case-insensitive).", call. = FALSE)
  )
  level <- tolower(level)
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
  
  
  # Align samples across table, env_data, and (optionally) metadata
  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata)
    common_samples <- Reduce(intersect, list(colnames(table), rownames(env_data), metadata[[1]]))
    table <- table[, c(tax_col, match(common_samples, colnames(table))), drop = FALSE]
    env_data <- env_data[common_samples, , drop = FALSE]
    metadata <- metadata[metadata[[1]] %in% common_samples, , drop = FALSE]
    rownames(metadata) <- metadata[[1]]
    table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
    ordered_samples <- intersect(metadata[, 1], colnames(table)[-1])
    table <- table[, c("taxonomy", ordered_samples)]
  } else {
    common_samples <- intersect(colnames(table), rownames(env_data))
    table <- table[, c(tax_col, match(common_samples, colnames(table))), drop = FALSE]
    env_data <- env_data[common_samples, , drop = FALSE]
    table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
    ordered_samples <- common_samples
    table <- table[, c("taxonomy", ordered_samples)]
  }
  
  # --- Ensure unique identifiers ---
  if (!is.null(rownames(table))) {
    table <- tibble::rownames_to_column(table, var = "OTU_ID")
  } else {
    table$OTU_ID <- paste0("OTU_", seq_len(nrow(table)))
  }
  
  # --- Strip trailing empty taxonomy levels ---
  table$taxonomy <- gsub("(;__)+$", "", table$taxonomy)
  
  # --- Determine the index of the taxonomic level ---
  level_idx <- switch(level,
                      kingdom = 1, phylum = 2, class = 3, order = 4,
                      family = 5, genus = 6, species = 7)
  
  # --- Compute the taxonomic depth of each OTU ---
  get_depth <- function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    sum(grepl("__", levels))
  }
  table$depth <- vapply(table$taxonomy, get_depth, integer(1))

  # --- Split rows by taxonomic resolution ---
  lowres <- table[table$depth < level_idx, ]    
  highres <- table[table$depth >= level_idx, ]  
  
  # --- Trim and collapse taxonomies with sufficient resolution ---
  highres$taxonomy <- vapply(highres$taxonomy, function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    paste(levels[seq_len(level_idx)], collapse = ";")
  }, character(1))
  
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
  #this way the plot shows only that name instead of the full taxonomy string
  
  
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
  
  # 1. Convert row names into an OTU_ID column
  table <- table %>%
    tibble::rownames_to_column(var = "OTU_ID")

  # 2. Create a new table called taxon with OTU_ID and taxonomy
  taxon <- table %>%
    dplyr::select(OTU_ID, taxonomy)

  # 3. Remove the taxonomy column from table
  table <- table %>%
    dplyr::select(-taxonomy)

  # 4. Convert OTU_ID back into row names
  table <- table %>%
    tibble::column_to_rownames(var = "OTU_ID")
  
  
  
  # Default palette
  if (is.null(col_palette)) {
    if (identical(diverging_palette, "viridis")) {
      # Sequential, colorblind-friendly viridis scale (same family used for
      # the heatmap body in aldex_heatmap_plot's "Median clr value")
      col_palette <- viridis::viridis(200, option = "C", direction = -1)
    } else {
      preset <- .mbm_div_palettes[[diverging_palette]]
      if (is.null(preset)) {
        warning("Unknown diverging_palette '", diverging_palette,
                "'. Using 'PuYl'. Valid options: 'viridis', ",
                paste(names(.mbm_div_palettes), collapse = ", "))
        preset <- .mbm_div_palettes[["PuYl"]]
      }
      col_palette <- grDevices::colorRampPalette(preset)(200)
    }
  }
  # --- Filas comunes
  common_samples <- base::intersect(colnames(table), rownames(env_data))
  counts <- table[, common_samples, drop = FALSE]
  env <- env_data[common_samples, , drop = FALSE]
  
  
  # Seleccionar solo las variables ambientales indicadas en env_vars, verificando coincidencias
  if (!is.null(env_vars)) {
    env_vars <- env_vars[env_vars %in% colnames(env)]
    if(length(env_vars) == 0) stop("No matching variables found in env_data")
    env <- env[, env_vars, drop = FALSE]
  }
  
  
  # Descartar columnas no numericas (p.ej. IDs o variables categoricas
  # mezcladas en env_data); stats::cor() requiere que env sea todo numerico.
  is_num <- vapply(env, is.numeric, logical(1))
  if (!all(is_num)) {
    warning("Dropping non-numeric columns from `env_data`: ",
            paste(names(env)[!is_num], collapse = ", "))
    env <- env[, is_num, drop = FALSE]
  }

  # Filtrar variables constantes
  env <- env[, apply(env, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  abund <- sweep(counts, 2, colSums(counts, na.rm = TRUE), FUN = "/") * 100
  abund <- abund[apply(abund, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]
  
  # Overall correlation matrix
  corr_mat <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
  
  
  
  
  # --- Calcular p-values si se indica pval_threshold
  pval_mat <- NULL
  if (!is.null(pval_threshold) || save_table) {
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
    
    # Every taxon x variable pair is a separate test, so correct for
    # multiple comparisons over the whole matrix before filtering.
    pval_raw <- pval_mat
    pval_mat[] <- stats::p.adjust(pval_mat, method = p_adjust_method)
  }
  if (!is.null(pval_threshold)) {

    # Mantener solo taxones significativos en al menos una variable
    signif_taxa <- rownames(abund)[apply(pval_mat, 2, function(x) any(x < pval_threshold, na.rm = TRUE))]
    if (length(signif_taxa) == 0) {
      stop("No taxa have a significant correlation (p < ", pval_threshold,
           " after p_adjust_method = '", p_adjust_method, "').\n",
           "Try a higher pval_threshold, a coarser level, or p_adjust_method = 'none'.",
           call. = FALSE)
    }
    abund <- abund[signif_taxa, , drop = FALSE]
    corr_mat <- stats::cor(env, t(abund), method = method, use = "pairwise.complete.obs")
  }
  
  
  
  # Hierarchical clustering
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
  # --- Add taxonomy to the melted results ---
  # corr_df$Group holds the OTU_IDs; join with the taxon table
  corr_df <- dplyr::left_join(
    corr_df,
    taxon,
    by = c("Group" = "OTU_ID")
  )
  
  # Replace "Group" with the taxonomic name
  corr_df$Taxon <- corr_df$taxonomy
  corr_df$taxonomy <- NULL

  # Save the plotted correlations (long format, with taxon names and p-values)
  if (save_table) {
    ids <- as.character(corr_df$Group)
    vars <- as.character(corr_df$Environmental)
    out <- data.frame(
      Taxon       = corr_df$Taxon,
      Variable    = vars,
      Correlation = corr_df$Correlation,
      p_value     = pval_raw[cbind(vars, ids)],
      p_adj       = pval_mat[cbind(vars, ids)],
      check.names = FALSE
    )
    utils::write.table(out, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message("Table saved as: ", table_filename)
  }
  # Base plot
  if (invert_axes) {
    p <- ggplot2::ggplot(corr_df,
                         ggplot2::aes(x = Environmental, y = Taxon, fill = Correlation))
  } else {
    p <- ggplot2::ggplot(corr_df,
                         ggplot2::aes(x = Taxon, y = Environmental, fill = Correlation))
  }
  
  # Choose plot type: "tile" behaves like a heatmap, "circle" like a bubble plot
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
    # Taxon names are italicized, matching every other function that labels
    # an axis with taxa (ancombc_plot, abundance_heatmap_plot,
    # ratios_bubble_plot, random_forest_lollipop_plot) - Taxon sits on
    # whichever axis invert_axes puts it on, so the 45deg/italic styling
    # (meant for the often-long taxon names) follows it there instead of
    # always defaulting to axis.text.x.
    .mbm_theme(
      legend_position = "right",
      extra = ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = x_label_angle,
          vjust = if (x_label_angle == 0) 0.5 else 1,
          hjust = if (x_label_angle == 0) 0.5 else 1,
          size  = 12, color = "black",
          face  = if (invert_axes) "plain" else "italic"
        ),
        axis.text.y = ggplot2::element_text(
          size = 12, color = "black",
          face = if (invert_axes) "italic" else "plain"
        ),
        panel.grid.major = ggplot2::element_blank()
      )
    ) +
    ggplot2::coord_fixed()
  
  return(p)
}
