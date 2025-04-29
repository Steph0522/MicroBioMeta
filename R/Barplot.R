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
#' @param group_var Character. Column name in `metadata` used to group samples for plotting (e.g., replicate, subject).
#' @param label Character. Legend title for the taxa groups. Default is `"taxonomy"`.
#' @param top_n_groups Integer. Number of most abundant taxa groups to display. Default is `15`.
#' 
#' @return A `ggplot2` object showing a stacked barplot of relative abundances.
#'
#' @details
#' - Relative abundances are calculated per sample (%).
#' - Taxa names are collapsed to the specified taxonomic `level` ("genus" or "phylum").
#' - Only the top `top_n_groups` taxa are shown; others are filtered out.
#' - Samples are grouped and ordered according to `group_var`.
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
                                    group_var,
                                    label = "taxonomy",
                                    top_n_groups = 15,
                                    name_vector = NULL) {
  library(ggplot2)
  library(ggtext)
  library(dplyr)
  library(tidyr)
  library(forcats)
  
  # Remove uninformative taxonomy strings
  table <- table %>%
    filter(taxonomy != "d__Bacteria;__;__;__;__;__")
  
  # Collapse to genus level if specified
  if (level == "genus") {
    table$taxonomy <-
      gsub(";\\s?s__.*", "", table$taxonomy)  # Remove species-level annotations
    table <- table %>%
      group_by(taxonomy) %>%
      summarise(across(where(is.numeric), sum, na.rm = TRUE))
  }
  
  # Collapse to phylum level if specified
  if (level == "phylum") {
    table$taxonomy <-
      gsub(";\\s?c__.*", "", table$taxonomy)  # Remove class-level annotations
    table <- table %>%
      group_by(taxonomy) %>%
      summarise(across(where(is.numeric), sum, na.rm = TRUE))
  }
  
  # Calculate relative abundance (%)
  table[,-1] <-
    sweep(table[,-1], 2, colSums(table[,-1], na.rm = TRUE), FUN = "/") * 100
  
  # Convert table to long format
  table_long <- table %>%
    pivot_longer(cols = -taxonomy,
                 names_to = "SAMPLEID",
                 values_to = "RelativeAbundance")
  
  # Join with metadata
  columns_to_join <- c("SAMPLEID", x_col, facet_col, group_var)
  columns_to_join <-
    columns_to_join[!is.na(columns_to_join) &
                      columns_to_join != "NULL"]
  
  table_long <- table_long %>%
    left_join(metadata %>% select(all_of(columns_to_join)), by = "SAMPLEID")
  
  # Group by required variables
  grouping_vars <- c(group_var, "taxonomy")
  if (!is.null(facet_col)) {
    grouping_vars <- c(grouping_vars, facet_col)
  }
  
  avg_by_group <- table_long %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(MeanAbundance = mean(RelativeAbundance, na.rm = TRUE),
              .groups = "drop")
  
  # Identify top N abundant taxa
  overall_means <- avg_by_group %>%
    group_by(taxonomy) %>%
    summarise(MeanAbundance = mean(MeanAbundance, na.rm = TRUE))
  
  top_groups <- overall_means %>%
    arrange(desc(MeanAbundance)) %>%
    slice_head(n = top_n_groups) %>%
    pull(taxonomy)
  
  avg_by_group <- avg_by_group %>%
    filter(taxonomy %in% top_groups)
  
  # Simplify taxonomic names for "silva" database at genus level
  if (taxonomy_db == "silva" && level == "genus") {
    avg_by_group <- avg_by_group %>%
      mutate(
        taxonomy = case_when(
          grepl("g__[^;]*", taxonomy) &
            !grepl("g__uncultured|g__$", taxonomy) ~ sub(".*g__([^;]*).*", "\\1", taxonomy),
          grepl("f__[^;]*", taxonomy) &
            !grepl("f__uncultured|f__$", taxonomy) ~ paste0(
              "other ",
              str_extract(taxonomy, "f__[^;]*") %>% sub("f__", "", .)
            ),
          grepl("o__[^;]*", taxonomy) &
            !grepl("o__uncultured|o__$", taxonomy) ~ paste0(
              "other ",
              str_extract(taxonomy, "o__[^;]*") %>% sub("o__", "", .)
            ),
          grepl("c__[^;]*", taxonomy) &
            !grepl("c__uncultured|c__$", taxonomy) ~ paste0(
              "other ",
              str_extract(taxonomy, "c__[^;]*") %>% sub("c__", "", .)
            ),
          grepl("p__[^;]*", taxonomy) &
            !grepl("p__uncultured|p__$", taxonomy) ~ paste0(
              "other ",
              str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)
            ),
          TRUE ~ "Unclassified"
        )
      )
  }
  
  # Simplify taxonomic names for "silva" database at phylum level
  if (taxonomy_db == "silva" && level == "phylum") {
    avg_by_group <- avg_by_group %>%
      mutate(taxonomy = case_when(
        grepl("p__[^;]*", taxonomy) ~ str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
        TRUE ~ "Unclassified"
      ))
  }
  
  # Set facet levels if applicable
  if (!is.null(facet_col)) {
    facet_levels <- unique(avg_by_group[[facet_col]])
    avg_by_group[[facet_col]] <-
      factor(avg_by_group[[facet_col]], levels = facet_levels)
  }
  
  # Set x-axis order based on metadata
  x_levels <- unique(metadata[[group_var]])
  avg_by_group[[group_var]] <-
    factor(avg_by_group[[group_var]], levels = x_levels)
  
  # Order taxonomy by maximum abundance
  taxonomy_order <- avg_by_group %>%
    group_by(taxonomy) %>%
    summarise(max_abund = max(MeanAbundance, na.rm = TRUE)) %>%
    arrange(desc(max_abund)) %>%
    pull(taxonomy)
  
  avg_by_group$taxonomy <-
    factor(avg_by_group$taxonomy, levels = rev(taxonomy_order))
  
  # Define color palette
  cbPalette <- colorRampPalette(
    c(
      "#999999",
      "#0099CC",
      "#ff6600",
      "#FF0066",
      "#99FF33",
      "#CC00cc",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#ff9900",
      "#56B4E9",
      "#FFFFFF",
      "#99ff90",
      "#ffff00",
      "#FF0000"
    )
  )(length(unique(avg_by_group$taxonomy)))
  
  # Generate plot
  p <- ggplot(avg_by_group,
              aes(
                fill = taxonomy,
                y = MeanAbundance,
                x = !!sym(group_var)
              )) +
    geom_bar(
      position = "stack",
      stat = "identity",
      width = 0.5,
      color = "#000000"
    ) +
    scale_fill_manual(name = label, values = cbPalette) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_text(size = 12),
      legend.text = element_markdown(size = 10)
    ) +
    ylim(0, 100) +
    ylab("Relative abundance (%)") +
    xlab("Sample")
  
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(vars(!!sym(facet_col)), scales = "free_x")
  }
  
  return(p)
}
