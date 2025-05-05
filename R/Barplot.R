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
#' @param top_n_groups Integer. Number of most abundant taxa groups to display. Default is `15`.
#' @param x_axis_title Character. The tittle that should be in the x-axis (deault = "Samples")
#' @return A `ggplot2` object showing a stacked barplot of relative abundances.
#'
#' @details
#' - Relative abundances are calculated per sample (%).
#' - Taxa names are collapsed to the specified taxonomic `level` ("genus" or "phylum").
#' - Only the top `top_n_groups` taxa are shown; others are filtered out.
#' - Samples are grouped and ordered according to `x_col`.
#' - Optional faceting by `facet_col` if provided.
#' - Taxonomic strings matching `"d__Bacteria;__;__;__;__;__"` are automatically removed.
#'
#' @export
#'
#' @examples
#' # relative_abundance_plot(table = your_table, metadata = your_metadata, ...)
relative_abundance_plot <- function(table,
                                    metadata,
                                    taxonomy_db = "silva",
                                    level = "genus",
                                    x_col,
                                    facet_col = NULL,
                                    label = "taxonomy",
                                    top_n_groups = 15,
                                    x_axis_title="Samples") {

  table<- table[,-1]
  table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
  # Remove uninformative taxonomy strings
  table <- table %>%
    dplyr::filter(taxonomy != "d__Bacteria;__;__;__;__;__")
  
  # Reorder columns based on SAMPLEID order in metadata
  ordered_samples <- metadata$SAMPLEID
  sample_columns <- colnames(table)[-1]
  ordered_samples <- intersect(ordered_samples, sample_columns)
  table <- table[, c("taxonomy", ordered_samples)]
  
  # Collapse to genus level if specified
  if (level == "genus") {
    table$taxonomy <- gsub(";\\s?s__.*", "", table$taxonomy)
    table <- table %>%
      dplyr::group_by(taxonomy) %>%
      dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE))
  }
  
  # Collapse to phylum level if specified
  if (level == "phylum") {
    table$taxonomy <- gsub(";\\s?c__.*", "", table$taxonomy)
    table <- table %>%
      dplyr::group_by(taxonomy) %>%
      dplyr::summarise(dplyr::across(where(is.numeric), sum, na.rm = TRUE))
  }
  
  # Calculate relative abundance (%)
  table[,-1] <- sweep(table[,-1], 2, colSums(table[,-1], na.rm = TRUE), FUN = "/") * 100
  
  # Convert to long format
  table_long <- table %>%
    tidyr::pivot_longer(cols = -taxonomy,
                        names_to = "SAMPLEID",
                        values_to = "RelativeAbundance")
  
  # Join with metadata
  columns_to_join <- c("SAMPLEID", x_col, facet_col, x_col)
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
    dplyr::slice_head(n = top_n_groups) %>%
    dplyr::pull(taxonomy)
  
  avg_by_group <- avg_by_group %>%
    dplyr::filter(taxonomy %in% top_groups)
  
  # Taxonomic name simplification for SILVA
  if (taxonomy_db == "silva" && level == "genus") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(
        taxonomy = dplyr::case_when(
          grepl("g__[^;]*", taxonomy) & !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
          grepl("f__[^;]*", taxonomy) & !grepl("f__uncultured|f__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)),
          grepl("o__[^;]*", taxonomy) & !grepl("o__uncultured|o__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)),
          grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)),
          grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  if (taxonomy_db == "silva" && level == "phylum") {
    avg_by_group <- avg_by_group %>%
      dplyr::mutate(taxonomy = dplyr::case_when(
        grepl("p__[^;]*", taxonomy) ~ stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
        TRUE ~ "Unclassified"
      ))
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
  
  avg_by_group$taxonomy <- factor(avg_by_group$taxonomy, levels = rev(taxonomy_order))
  
  cbPalette <- grDevices::colorRampPalette(
    c(
      "#999999", "#0099CC", "#ff6600", "#FF0066", "#99FF33",
      "#CC00cc", "#009E73", "#F0E442", "#0072B2", "#ff9900",
      "#56B4E9", "#FFFFFF", "#99ff90", "#ffff00", "#FF0000"
    )
  )(length(unique(avg_by_group$taxonomy)))
  
  p <- ggplot2::ggplot(avg_by_group,
                       ggplot2::aes(
                         x = !!rlang::sym(x_col),
                         y = MeanAbundance,
                         fill = taxonomy
                       )) +
    ggplot2::geom_bar(
      position = "stack",
      stat = "identity",
      width = 0.5,
      color = "#000000"
    ) +
    ggplot2::scale_fill_manual(name = label, values = cbPalette) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 14, color = "black"),
      axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
      axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
      legend.text = ggplot2::element_text(size = 10, face = if (level == "genus") "italic" else "plain")
    ) +
    ggplot2::ylim(0, 100) +
    ggplot2::ylab("Relative abundance (%)") +
    ggplot2::xlab(x_axis_title)
  
  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet_col)), scales = "free_x")
  }
  
  return(p)
}
