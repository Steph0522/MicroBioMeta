#' Relative abundance barplot
#'
#' Generates a bar plot of relative abundance (%) for the most abundant taxa groups across samples or sample groups.
#'
#' @param table A data frame with taxa in rows and samples in columns.
#' The last column must be named `taxonomy`, containing full taxonomic strings.
#' @param metadata A data frame containing sample metadata. Its first column
#'   must hold the sample IDs (the column names of `table`). Optional: if
#'   `NULL` (default), one bar is drawn per sample.
#' @param taxonomy_db Character. Database the taxonomy strings come from:
#'   \code{"silva"} (default), \code{"gg2"} (Greengenes2, also \code{"gg"}),
#'   \code{"unite"} or \code{"Kraken2"} (also \code{"kraken"}). Case-insensitive.
#' @param level Character. Taxonomic level: \code{"kingdom"}, \code{"phylum"},
#'   \code{"class"}, \code{"order"}, \code{"family"}, \code{"genus"} (default) or
#'   \code{"species"}. Case-insensitive.
#' @param x_col Character. Name of the column in \code{metadata} for the
#'   x-axis. If \code{NULL} (default), one bar per sample.
#' @param facet_by Character. Name of the column in \code{metadata} to facet
#'   the plot by. Optional; \code{NULL} (default) for no facets.
#' @param label Character. Title of the taxa legend. Default \code{"taxonomy"}.
#' @param top_n Integer. Number of most abundant taxa to show. Default
#'   \code{15}.
#' @param x_axis_title Character. Title of the x-axis. Default
#'   \code{"Samples"}.
#' @param y_axis_title Character. Title of the y-axis. Default \code{"Relative
#'   abundance (\%)"}.
#' @param x_label_angle Numeric. Rotation (degrees) of the x-axis labels.
#'   Default \code{0} (horizontal); use \code{45} or \code{90} if they overlap.
#' @param strip_text_bold Logical. If \code{TRUE}, the facet strip labels are
#'   bold. Default \code{FALSE}.
#' @param strip_color Character. Background color of the facet strips (used
#'   when \code{facet_by} is set). Default \code{"grey"}.
#' @param aspect_ratio Numeric. Aspect ratio (height/width) of each panel.
#'   Default \code{NULL} (automatic).
#' @param add_remained Logical. If \code{TRUE}, adds an "Other" category with
#'   the sum of the remaining taxa. Default \code{FALSE}.
#' @param width_equal Logical. If \code{TRUE}, all bars have equal width regardless of sample count per group. Default \code{FALSE}.
#' @param save_table Logical. If \code{TRUE}, saves the relative-abundance
#'   table as a tab-delimited file. Default \code{FALSE}.
#' @param table_filename Character. Name or path of the saved file (used when
#'   \code{save_table = TRUE}). Default \code{"relative_abundance.txt"}.
#' @return A ggplot object with the stacked bar plot of relative abundances.
#'
#' @details
#' - Relative abundances are calculated per sample (%).
#' - Taxa are collapsed to the taxonomic `level`; taxa not resolved to that level are shown as "other <rank>".
#' - Only the top `top_n` taxa are shown; others are filtered out.
#' - Samples are grouped and ordered according to `x_col`.
#' - Optional faceting by `facet_by` if provided.
#' - Features assigned only to a kingdom/domain (e.g. `"d__Bacteria;__;__;__;__;__"`) or `"Unassigned"` are removed.
#'
#' @export
#'
#' @examples
#' table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
#' table <- read.delim(table_path, row.names = 1, check.names = FALSE)
#'
#' metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
#' metadata <- read.delim(metadata_path, check.names = FALSE)
#'
#' abundance_bar_plot(
#'     table = table,
#'     metadata = metadata,
#'     taxonomy_db = "silva",
#'     level = "genus",
#'     x_col = "Location",
#'     label = "Genus",
#'     facet_by = "Treatment",
#'     top_n = 30,
#'     add_remained = TRUE
#' )
abundance_bar_plot <- function(table,
                               metadata = NULL,
                               taxonomy_db = "silva",
                               level = "genus",
                               x_col = NULL,
                               facet_by = NULL,
                               width_equal = FALSE,
                               label = "taxonomy",
                               top_n = 15,
                               x_axis_title = "Samples",
                               y_axis_title = "Relative abundance (%)",
                               x_label_angle = 0,
                               strip_text_bold = FALSE,
                               strip_color = "grey",
                               aspect_ratio = NULL,
                               add_remained = FALSE,
                               save_table = FALSE,
                               table_filename = "relative_abundance.txt") {
    taxonomy_db <- .mbm_taxonomy_db(taxonomy_db)
    level <- .mbm_check_level(level)

    .mbm_check_table(table)
    tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
    if (length(tax_col) != 1) stop("There is no taxonomy column in the table")

    if (is.null(metadata)) {
        if (!is.null(facet_by)) stop("`facet_by` needs `metadata`.", call. = FALSE)
        metadata <- data.frame(SAMPLEID = names(table)[-tax_col])
    }
    colnames(metadata)[1] <- "SAMPLEID"
    if (is.null(x_col)) x_col <- "SAMPLEID"
    if (!x_col %in% names(metadata)) {
        stop("`x_col` '", x_col, "' is not a column of metadata.", call. = FALSE)
    }

    names(table)[tax_col] <- "taxonomy"

    table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]

    table <- table %>%
        dplyr::filter(taxonomy != "d__Bacteria;__;__;__;__;__") %>%
        dplyr::filter(taxonomy != "d__Bacteria") %>%
        dplyr::filter(taxonomy != "d__Archaea;__;__;__;__;__") %>%
        dplyr::filter(taxonomy != "d__Archaea") %>%
        dplyr::filter(taxonomy != "d__Bacteria;p__;c__;o__;f__;g__;s__") %>%
        dplyr::filter(taxonomy != "d__Archaea;p__;c__;o__;f__;g__;s__") %>%
        dplyr::filter(taxonomy != "k__Bacteria;__;__;__;__;__") %>%
        dplyr::filter(taxonomy != "k__Fungi;__;__;__;__;__") %>%
        dplyr::filter(taxonomy != "k__Fungi;p__;c__;o__;f__;g__") %>%
        dplyr::filter(taxonomy != "k__Fungi") %>%
        dplyr::filter(taxonomy != "Unassigned") %>%
        dplyr::filter(taxonomy != "d__Eukaryota")

    table <- table %>%
        dplyr::filter(!grepl("^(d__|k__)[^;]*;[ _;]*$", taxonomy))

    if (taxonomy_db == "gg2") {
        table <- table %>%
            dplyr::filter(!grepl("(__;?)+$", taxonomy))
    }

    ordered_samples <- metadata[[1]]
    sample_columns <- colnames(table)[-1]
    ordered_samples <- intersect(ordered_samples, sample_columns)
    table <- table[, c("taxonomy", ordered_samples)]
    table <- .mbm_truncate_taxonomy(table, level, taxonomy_db)

    table <- table %>%
        dplyr::group_by(taxonomy) %>%
        dplyr::summarise(dplyr::across(where(is.numeric), \(x) sum(x, na.rm = TRUE)))

    table[, -1] <- sweep(table[, -1], 2, colSums(table[, -1], na.rm = TRUE), FUN = "/") * 100

    if (save_table) {
        utils::write.table(
            table,
            file = table_filename,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        message("Table saved as: ", table_filename)
    }

    table_long <- table %>%
        tidyr::pivot_longer(cols = -taxonomy, names_to = "SAMPLEID", values_to = "RelativeAbundance")

    columns_to_join <- unique(c("SAMPLEID", x_col, facet_by))
    columns_to_join <- columns_to_join[!is.na(columns_to_join) & columns_to_join != "NULL"]
    table_long <- dplyr::left_join(
        table_long,
        metadata %>% dplyr::select(dplyr::all_of(columns_to_join)),
        by = "SAMPLEID"
    )

    grouping_vars <- c(x_col, "taxonomy")
    if (!is.null(facet_by)) {
        grouping_vars <- c(grouping_vars, facet_by)
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

    avg_by_group <- .mbm_label_taxa(avg_by_group, level, taxonomy_db)

    avg_by_group$taxonomy <- .mbm_composite_genus(avg_by_group$taxonomy)

    avg_by_group <- avg_by_group %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
        dplyr::summarise(
            MeanAbundance = sum(MeanAbundance, na.rm = TRUE),
            .groups = "drop"
        )

    if (!is.null(facet_by)) {
        facet_levels <- unique(avg_by_group[[facet_by]])
        avg_by_group[[facet_by]] <- factor(avg_by_group[[facet_by]], levels = facet_levels)
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

    tax_levels <- levels(avg_by_group$taxonomy)
    tax_levels_no_other <- setdiff(tax_levels, c("Other", "Unclassified"))

    cbPalette <- grDevices::colorRampPalette(
        c(
            "#99ff10", "#0099CC", "#ff6600", "#FF0066", "#99FF33",
            "#CC00cc", "#009E73", "#F0E442", "#0072B2", "#ff9900",
            "#56B4E9", "#FFFFFF", "#99ff90", "#ffff00", "#FF0000"
        )
    )(length(tax_levels_no_other))

    names(cbPalette) <- tax_levels_no_other
    cbPalette["Other"] <- "#D3D3D3"
    cbPalette["Unclassified"] <- "#666666"

    p <- ggplot2::ggplot(
        avg_by_group,
        ggplot2::aes(
            x = !!rlang::sym(x_col),
            y = MeanAbundance,
            fill = taxonomy
        )
    ) +
        ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.5, color = "#000000") +
        ggplot2::scale_fill_manual(
            name = label,
            values = cbPalette,
            labels = function(taxa) {
                if (level %in% c("genus", "species")) {
                    lapply(taxa, function(x) {
                        if (x %in% c("Other", "Unclassified") || grepl("^other ", x)) {
                            x
                        } else {
                            bquote(italic(.(x)))
                        }
                    })
                } else {
                    taxa
                }
            }
        ) +
        .mbm_theme(
            legend_position = "right",
            extra = ggplot2::theme(
                panel.grid        = ggplot2::element_blank(),
                axis.text.x       = .mbm_x_text(x_label_angle),
                strip.text        = .mbm_strip_text(strip_text_bold),
                strip.background  = ggplot2::element_rect(fill = strip_color, color = "black")
            )
        ) +
        ggplot2::coord_cartesian(ylim = c(0, 100)) +
        ggplot2::ylab(y_axis_title) +
        ggplot2::xlab(x_axis_title)

    if (!is.null(facet_by)) {
        strip_labeller <- ggplot2::labeller(
            .default = ggplot2::label_wrap_gen(width = 10)
        )
        if (width_equal) {
            p <- p + ggplot2::facet_grid(
                rows = NULL, cols = vars(!!rlang::sym(facet_by)),
                scales = "free_x", space = "free",
                labeller = strip_labeller
            )
        } else {
            p <- p + ggplot2::facet_wrap(ggplot2::vars(!!rlang::sym(facet_by)),
                scales = "free_x",
                labeller = strip_labeller
            )
        }
    }

    if (!is.null(aspect_ratio)) {
        p <- p + ggplot2::theme(aspect.ratio = aspect_ratio)
    }

    return(p)
}
