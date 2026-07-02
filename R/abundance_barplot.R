#' Relative abundance barplot
#'
#' Generates a bar plot of relative abundance (%) for the most abundant taxa groups across samples or sample groups.
#'
#' @param table A data frame with taxa in rows and samples in columns. The first column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. Must include a `SAMPLEID` column matching sample names in `table`.
#' @param taxonomy_db Character. Reference taxonomy database: `"silva"` (default). Affects how taxonomic names are simplified in the plot.
#' @param level Character. Taxonomic level to collapse: `"genus"` (default) or `"phylum"`.
#' @param x_col Character. Column name in `metadata` to use for the x-axis (e.g., environment, condition).
#' @param facet_col Optional. Character. Column name in `metadata` to facet the plot by (e.g., treatment group). Default is `NULL`.
#' @param label Character. Legend title for the taxa groups. Default is `"taxonomy"`.
#' @param top_n Integer. Number of most abundant taxa groups to display. Default is `15`.
#' @param x_axis_title Character. The title for the x-axis (default = "Samples")
#' @param add_remained Logical indicating whether to include an "Other" category to sum remaining groups; default is FALSE.
#' @param width_equal Logical. If \code{TRUE}, all bars have equal width regardless of sample count per group. Default \code{FALSE}.
#' @param save_table Logical. If \code{TRUE}, saves the relative-abundance table to disk. Default \code{TRUE}.
#' @param table_filename Character. File path/name for the saved table (used when \code{save_table = TRUE}). Default \code{"relative_abundance.txt"}.
#' @return A `ggplot2` object showing a stacked barplot of relative abundances.
#'
#' @details
#' - Relative abundances are calculated per sample (%).
#' - Taxa names are collapsed to the specified taxonomic `level` ("genus" or "phylum").
#' - Only the top `top_n` taxa are shown; others are filtered out.
#' - Samples are grouped and ordered according to `x_col`.
#' - Optional faceting by `facet_col` if provided.
#' - Taxonomic strings matching `"d__Bacteria;__;__;__;__;__"` are automatically removed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' abundance_barplot(
#'   table      = table_taxa,
#'   metadata   = metadata,
#'   taxonomy_db = "silva",
#'   level      = "genus",
#'   x_col      = "MUESTRA",
#'   label = "Genus",
#'   facet_col  = "SITIO",
#'   top_n      = 30,
#'   add_remained = TRUE
#' )
#' }



abundance_barplot <- function(table,
                              metadata,
                              taxonomy_db = "silva",
                              level = "genus",
                              x_col,
                              facet_col = NULL,
                              width_equal = FALSE,
                              label = "taxonomy",
                              top_n = 15,
                              x_axis_title = "Samples",
                              add_remained = FALSE,
                              save_table = TRUE,
                              table_filename = "relative_abundance.txt") {
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  names(table)[ncol(table)] <- "taxonomy"
  
  table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
  
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
  
  table <- table %>%
  # quitar taxonomías que solo llegan al dominio/reino
  dplyr::filter(!grepl("^(d__|k__)[^;]*;[ _;]*$", taxonomy))
  
  
  if (taxonomy_db == "gg2") {
    table <- table %>%
      # quitar taxonomías con todos los niveles vacíos (__)
      dplyr::filter(!grepl("(__;?)+$", taxonomy))
  }
  
  
  # Reorder columns based on SAMPLEID order in metadata
  ordered_samples <- metadata[[1]]
  sample_columns <- colnames(table)[-1]
  ordered_samples <- intersect(ordered_samples, sample_columns)
  table <- table[, c("taxonomy", ordered_samples)]
  if (taxonomy_db %in% c("silva", "Kraken2", "gg2")) {
    if (level == "kingdom") {
      table$taxonomy <- sub(";.*", "", table$taxonomy)
    }
    if (level == "phylum") {
      table$taxonomy <- sub(";\\s?c__.*", "", table$taxonomy)
    }
    if (level == "class") {
      table$taxonomy <- sub(";\\s?o__.*", "", table$taxonomy)
    }
    if (level == "order") {
      table$taxonomy <- sub(";\\s?f__.*", "", table$taxonomy)
    }
    if (level == "family") {
      table$taxonomy <- sub(";\\s?g__.*", "", table$taxonomy)
    }
    if (level == "genus") {
      table$taxonomy <- sub(";\\s?s__.*", "", table$taxonomy)
    }
    if (level == "species") {
      table$taxonomy <- table$taxonomy
    }
  }
  
  # ---- Colapsar taxonomía según nivel ----
  if (taxonomy_db == "unite") {
    if (level == "kingdom") {
      table$taxonomy <- sub(";.*", "", table$taxonomy)
    }
    if (level == "phylum") {
      table$taxonomy <- sub(";\\s?c__.*", "", table$taxonomy)
    }
    if (level == "class") {
      table$taxonomy <- sub(";\\s?o__.*", "", table$taxonomy)
    }
    if (level == "order") {
      table$taxonomy <- sub(";\\s?f__.*", "", table$taxonomy)
    }
    if (level == "family") {
      table$taxonomy <- sub(";\\s?g__.*", "", table$taxonomy)
    }
    if (level == "genus") {
      table$taxonomy <- sub(";\\s?s__.*", "", table$taxonomy)
    }
    if (level == "species") {
      table$taxonomy <- sub(";\\s?sh__.*", "", table$taxonomy)
    }
  }
  
  table <- table %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE))
  
  # Calculate relative abundance
  table[,-1] <- sweep(table[,-1], 2, colSums(table[,-1], na.rm = TRUE), FUN = "/") * 100
  
  # Guardar tabla si se solicita
  if (save_table) {
    utils::write.table(
      table,
      file = table_filename,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    ) 
    message(paste("Table saved as:", table_filename))
  }
  
  # Convert to long format
  table_long <- table %>%
    tidyr::pivot_longer(cols = -taxonomy, names_to = "SAMPLEID", values_to = "RelativeAbundance")
  
  # Join with metadata
  columns_to_join <- c("SAMPLEID", x_col, facet_col)
  columns_to_join <- columns_to_join[!is.na(columns_to_join) & columns_to_join != "NULL"]
  table_long <- dplyr::left_join(
    table_long,
    metadata %>% dplyr::select(dplyr::all_of(columns_to_join)),
    by = "SAMPLEID"
  )
  
  # Grouping
  grouping_vars <- c(x_col, "taxonomy")
  if (!is.null(facet_col)) {
    grouping_vars <- c(grouping_vars, facet_col)
  }
  
  avg_by_group <- table_long %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
    dplyr::summarise(MeanAbundance = mean(RelativeAbundance, na.rm = TRUE), .groups = "drop")
  
  overall_means <- avg_by_group %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(MeanAbundance = mean(MeanAbundance, na.rm = TRUE))
  
  top_groups <- overall_means %>%
    dplyr::arrange(dplyr::desc(MeanAbundance)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(taxonomy)
  
  if (add_remained) {
    top_avg <- table_long %>%
      dplyr::filter(taxonomy %in% top_groups) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
      dplyr::summarise(MeanAbundance = mean(RelativeAbundance, na.rm = TRUE), .groups = "drop")
    
    grouping_vars_no_tax <- setdiff(grouping_vars, "taxonomy")
    summed <- top_avg %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars_no_tax))) %>%
      dplyr::summarise(SumAbundance = sum(MeanAbundance), .groups = "drop")
    
    other_rows <- summed %>%
      dplyr::mutate(MeanAbundance = pmax(0, 100 - SumAbundance)) %>%
      dplyr::mutate(taxonomy = "Other") %>%
      dplyr::select(all_of(grouping_vars), MeanAbundance)
    
    avg_by_group <- dplyr::bind_rows(top_avg, other_rows)
    
    all_combinations <- tidyr::expand_grid(
      taxonomy = unique(avg_by_group$taxonomy),
      !!!setNames(
        lapply(grouping_vars_no_tax, function(v) unique(avg_by_group[[v]])),
        grouping_vars_no_tax
      )
    )
    
    avg_by_group <- dplyr::right_join(all_combinations, avg_by_group, by = grouping_vars)
    avg_by_group$MeanAbundance[is.na(avg_by_group$MeanAbundance)] <- 0
  } else {
    avg_by_group <- avg_by_group %>%
      dplyr::filter(taxonomy %in% top_groups)
  }
  
  if (taxonomy_db %in% c("unite","silva", "gg2") && level == "species") {
    avg_by_group <- avg_by_group %>%
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
    avg_by_group <- avg_by_group %>%
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
    avg_by_group <- avg_by_group %>%
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
    avg_by_group <- avg_by_group %>%
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
    avg_by_group <- avg_by_group %>%
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
    avg_by_group <- avg_by_group %>%
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
    avg_by_group <- avg_by_group %>%
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
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          taxonomy == "Other" ~ "Other",
          grepl("p__[^;]*", taxonomy) ~ stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  if (!is.null(facet_col)) {
    facet_levels <- unique(avg_by_group[[facet_col]])
    avg_by_group[[facet_col]] <- factor(avg_by_group[[facet_col]], levels = facet_levels)
  }
  
  x_levels <- unique(metadata[[x_col]])
  avg_by_group[[x_col]] <- factor(avg_by_group[[x_col]], levels = x_levels)
  
  taxonomy_order <- avg_by_group %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(max_abund = max(MeanAbundance, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(max_abund)) %>%
    dplyr::pull(taxonomy)
  
  taxonomy_order <- unique(c(setdiff(taxonomy_order, c("Other", "Unclassified")), "Unclassified", "Other"))
  avg_by_group$taxonomy <- factor(avg_by_group$taxonomy, levels = rev(taxonomy_order))
  
  # ==== CORRECCIÓN PALLETA ====
  tax_levels <- levels(avg_by_group$taxonomy)
  tax_levels_no_other <- setdiff(tax_levels, c("Other", "Unclassified"))
  
  cbPalette <- grDevices::colorRampPalette(
    c("#99ff10", "#0099CC", "#ff6600", "#FF0066", "#99FF33",
      "#CC00cc", "#009E73", "#F0E442", "#0072B2", "#ff9900",
      "#56B4E9", "#FFFFFF", "#99ff90", "#ffff00", "#FF0000")
  )(length(tax_levels_no_other))
  
  names(cbPalette) <- tax_levels_no_other
  cbPalette["Other"] <- "#D3D3D3"
  cbPalette["Unclassified"] <- "#666666"
  
  # Plot
  p <- ggplot2::ggplot(avg_by_group,
                       ggplot2::aes(x = !!rlang::sym(x_col),
                                    y = MeanAbundance,
                                    fill = taxonomy)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5, color = "#000000") +
    ggplot2::scale_fill_manual(
      name = label,
      values = cbPalette,
      labels = function(taxa) {
        if (level %in% c("genus", "species")) {
          sapply(taxa, function(x) {
            if (x %in% c("Other", "Unclassified") || grepl("^other ", x)) {
              x   # texto plano
            } else {
              bquote(italic(.(x)))  # en cursivas
            }
          })
        } else {
          taxa  # siempre texto plano
        }
      }
      
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 14, color = "black", family = "serif", face = "bold"),
      strip.text = ggplot2::element_text(size = 12, color = "black", family = "serif", face = "bold"),
      axis.title.x = ggplot2::element_text(size = 14, color = "black", family = "serif"),
      axis.title.y = ggplot2::element_text(size = 14, color = "black", family = "serif"),
      axis.text.x = ggplot2::element_text(size = 12, color = "black", family = "serif"),
      axis.text.y = ggplot2::element_text(size = 12, color = "black", family = "serif"),
      legend.text = ggplot2::element_text(size = 12, family = "serif")
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 100)) +
    ggplot2::ylab("Relative abundance (%)") +
    ggplot2::xlab(x_axis_title)
  
  if (!is.null(facet_col)) {
    if (width_equal) {
      p <- p + ggplot2::facet_grid(rows = NULL, cols = vars(!!rlang::sym(facet_col)),
                                   scales = "free_x", space = "free")
    } else {
      p <- p + ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet_col)), scales = "free_x")
    }
  }
  
  return(p)
}
