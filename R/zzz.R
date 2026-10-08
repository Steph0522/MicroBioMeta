#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details. Re-exported
#' here so every function's \code{@examples} (most of which chain steps with
#' \code{\%>\%}) work with just \code{library(MicroBioMeta)}, without also
#' requiring \code{library(dplyr)}.
#'
#' @param lhs A value.
#' @param rhs A function call using the magrittr semantics.
#' @return \code{rhs(lhs)}, i.e. \code{lhs} piped into \code{rhs}.
#' @importFrom dplyr %>%
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @usage lhs \%>\% rhs
#' @examples
#' c(1, 2, 3) %>% sum()
NULL

utils::globalVariables(c(
    ".grp_idx", "Index", "Variable", "geo_km", "group", "hill", "hjust",
    "row_id", "side", "similarity", "slot", "vjust", "x_pos", "y_pos"
))

.mbm_parse_taxonomy <- function(taxonomy, tax_sep = "; |;", trim_extra = TRUE) {
    if (sum(colnames(taxonomy) %in% c("Feature.ID", "Taxon")) != 2) {
        stop("Table does not match expected format, i.e. does not have columns Feature.ID and Taxon.")
    }

    taxonomy <- taxonomy[, c("Feature.ID", "Taxon")]
    if (trim_extra) {
        taxonomy$Taxon <- gsub("[kdpcofgs]__", "", taxonomy$Taxon)
        taxonomy$Taxon <- gsub("D_\\d__", "", taxonomy$Taxon)
    }
    taxonomy <- suppressWarnings(
        tidyr::separate(taxonomy, Taxon,
            c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
            sep = tax_sep, fill = "right", extra = "merge"
        )
    )
    taxonomy <- apply(taxonomy, 2, function(x) ifelse(x == "", NA_character_, x))
    taxonomy <- as.data.frame(taxonomy)
    rownames(taxonomy) <- taxonomy$Feature.ID
    taxonomy$Feature.ID <- NULL
    taxonomy
}

.mbm_check_table <- function(table, arg = "table") {
    if (!is.data.frame(table) && !is.matrix(table)) {
        what <- if (is.function(table)) {
            "a function"
        } else {
            sprintf("of class '%s'", class(table)[1])
        }
        stop("`", arg, "` must be a data frame with taxa in rows and samples in ",
            "columns, but it is ", what,
            ". Check the name of the object you passed.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.mbm_aldex_clr <- function(counts, mc_samples = 1) {
    if (!is.numeric(mc_samples) || length(mc_samples) != 1 || mc_samples < 1) {
        stop("`mc_samples` must be a single number >= 1 (e.g. 1 or 128).", call. = FALSE)
    }
    aldex_obj <- withCallingHandlers(
        ALDEx2::aldex.clr(counts,
            mc.samples = mc_samples,
            denom = "all", verbose = FALSE, useMC = FALSE
        ),
        warning = function(w) {
            if (grepl("so few MC smps", conditionMessage(w))) invokeRestart("muffleWarning")
        }
    )
    if (mc_samples == 1) {
        return(t(ALDEx2::getMonteCarloSample(aldex_obj, 1)))
    }
    instances <- ALDEx2::getMonteCarloInstances(aldex_obj)
    clr_mean <- vapply(instances, rowMeans, numeric(nrow(instances[[1]])))
    rownames(clr_mean) <- rownames(instances[[1]])
    t(clr_mean)
}

.mbm_betadisper_spider_df <- function(betadisper_obj, groups) {
    groups <- factor(groups)

    n_axes <- length(vegan::eigenvals(betadisper_obj))
    if (n_axes < 2) {
        stop(
            "This dissimilarity partition only has ", n_axes, " usable ordination axis, ",
            "so a 2D spider plot can't be drawn for it. This happens when one component ",
            "of the beta-diversity partition (often nestedness) is close to zero across ",
            "the whole dataset - i.e. it's a property of this data, not a fixable bug.",
            call. = FALSE
        )
    }

    df_ord <- vegan::scores(betadisper_obj, display = "sites", choices = c(1, 2))
    df_ord <- data.frame(x = df_ord[, 1], y = df_ord[, 2], Group = groups)

    df_mean.ord <- stats::aggregate(df_ord[, c("x", "y")], by = list(Group = df_ord$Group), mean)

    df_spiders <- df_ord
    df_spiders$cntr.x <- df_mean.ord$x[match(df_spiders$Group, df_mean.ord$Group)]
    df_spiders$cntr.y <- df_mean.ord$y[match(df_spiders$Group, df_mean.ord$Group)]

    eig <- vegan::eigenvals(betadisper_obj)
    pos_sum <- sum(eig[eig > 0])
    axis_label <- function(i) {
        if (is.na(pos_sum) || pos_sum <= 0 || is.na(eig[i]) || eig[i] <= 0) {
            paste0("PCoA", i)
        } else {
            sprintf("PCoA%d (%.1f%%)", i, 100 * eig[i] / pos_sum)
        }
    }

    list(
        df_ord      = df_ord,
        df_mean.ord = df_mean.ord,
        df_spiders  = df_spiders,
        xlab        = axis_label(1),
        ylab        = axis_label(2)
    )
}

.mbm_format_pval <- function(p, accuracy = 0.001) {
    scales::label_pvalue(accuracy = accuracy)(p)
}

.mbm_p_label <- function(p, accuracy = 0.001) {
    scales::pvalue(p, accuracy = accuracy, add_p = TRUE)
}

.mbm_stat_layer <- function(stat, p_adjust_method = "holm",
                            data = NULL, label.x = NULL, label.y = NULL) {
    if (stat %in% c("wilcox.test", "t.test")) {
        p_col <- if (p_adjust_method == "none") "p" else "p.adj"
        ggpubr::stat_pwc(
            data            = data,
            method          = stat,
            p.adjust.method = p_adjust_method,
            p.adjust.by     = "panel",
            label           = paste0("{scales::pvalue(", p_col, ", 0.001, add_p = TRUE)}"),
            hide.ns         = p_col,
            size            = 0.4,
            label.size      = 3,
            family          = "serif",
            tip.length      = 0.01,
            step.increase   = 0.14,
            vjust           = -0.2
        )
    } else {
        ggpubr::stat_compare_means(
            method = stat,
            mapping = ggplot2::aes(
                label = scales::pvalue(after_stat(p), 0.001, add_p = TRUE)
            ),
            size = 3.5, family = "serif", hide.ns = TRUE,
            data = data, label.x = label.x, label.y = label.y
        )
    }
}

.mbm_patchwork_grid <- function(plots, ncol, tags = NULL, title = NULL,
                                show_legend = TRUE, legend_position = "bottom",
                                tag_bold = TRUE, nrow = NULL, widths = NULL,
                                heights = NULL) {
    p <- patchwork::wrap_plots(plots,
        ncol = ncol, nrow = nrow,
        widths = widths, heights = heights
    ) +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(
            title = title,
            tag_levels = if (!is.null(tags)) list(tags[seq_along(plots)]) else NULL,
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(
                    face = "bold", family = "serif",
                    size = 14, hjust = 0.5
                )
            )
        )
    p & ggplot2::theme(
        legend.position = if (show_legend) legend_position else "none",
        plot.tag = ggplot2::element_text(
            family = "serif", size = 14,
            face = if (tag_bold) "bold" else "plain"
        )
    )
}

.mbm_check_color_names <- function(group_colors, groups, arg = "group_colors") {
    if (!is.null(names(group_colors)) && !any(names(group_colors) %in% groups)) {
        warning("None of the names of `", arg, "` (", paste(names(group_colors), collapse = ", "),
            ") match the groups in the plot (", paste(groups, collapse = ", "),
            "), so they are drawn grey. Name the colors after these groups.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

.mbm_row_axis_title <- function(title, n_rows, width = 18) {
    if (n_rows > 1 && is.character(title) && length(title) == 1 && nchar(title) > width) {
        paste(strwrap(title, width = width), collapse = "\n")
    } else {
        title
    }
}

.mbm_facet_tags <- function(p, tags, bold = TRUE) {
    layout <- ggplot2::ggplot_build(p)$layout$layout
    layout <- layout[order(layout$ROW, layout$COL), , drop = FALSE]
    facet_vars <- setdiff(names(layout), c(
        "PANEL", "ROW", "COL", "SCALE_X", "SCALE_Y",
        "AXIS_X", "AXIS_Y", "COORD"
    ))
    tag_df <- layout[, facet_vars, drop = FALSE]
    tag_df$.tag <- tags[seq_len(nrow(tag_df))]
    p + ggplot2::geom_text(
        data = tag_df,
        mapping = ggplot2::aes(x = -Inf, y = Inf, label = .tag),
        inherit.aes = FALSE, hjust = -0.4, vjust = 1.4,
        family = "serif", size = 4.5, fontface = if (bold) "bold" else "plain"
    )
}

.mbm_align_metadata <- function(sample_ids, metadata) {
    sample_ids <- trimws(as.character(sample_ids))
    meta_ids <- trimws(as.character(metadata[[1]]))

    if (anyDuplicated(meta_ids)) {
        stop("Duplicated sample IDs in the first column of metadata: ",
            paste(head(unique(meta_ids[duplicated(meta_ids)])), collapse = ", "),
            call. = FALSE
        )
    }
    keep <- sample_ids[sample_ids %in% meta_ids]
    if (length(keep) == 0) {
        stop("None of the samples of the table are in the first column of metadata, ",
            "which must hold the sample IDs.\n",
            "Table samples (first 6): ", paste(head(sample_ids), collapse = ", "), "\n",
            "Metadata IDs (first 6): ", paste(head(meta_ids), collapse = ", "),
            call. = FALSE
        )
    }
    missing_ids <- setdiff(sample_ids, keep)
    if (length(missing_ids) > 0) {
        message(
            length(missing_ids), " sample(s) of the table are not in metadata ",
            "and were left out: ", paste(head(missing_ids), collapse = ", "),
            if (length(missing_ids) > 6) ", ..." else ""
        )
    }

    if (!identical(meta_ids[meta_ids %in% keep], keep)) {
        message("metadata rows were reordered to match the sample order of the table.")
    }
    aligned <- metadata[match(keep, meta_ids), , drop = FALSE]
    aligned[[1]] <- keep

    if (!identical(as.character(aligned[[1]]), keep)) {
        stop("Metadata could not be aligned with the table (internal problem).", call. = FALSE)
    }
    aligned
}

.mbm_env_data <- function(env_data, metadata, env_vars) {
    if (!is.null(env_data)) {
        return(env_data)
    }

    if (is.null(metadata)) {
        stop("Give `metadata` (with the environmental variables as columns) ",
            "or a separate `env_data` table.",
            call. = FALSE
        )
    }
    if (is.null(env_vars) || length(env_vars) == 0) {
        stop("`env_vars` is required when the environmental variables are taken ",
            "from `metadata` (env_data = NULL): name the columns to use.",
            call. = FALSE
        )
    }
    metadata <- as.data.frame(metadata)
    missing_vars <- setdiff(env_vars, colnames(metadata))
    if (length(missing_vars) > 0) {
        stop("Variables of `env_vars` not found in metadata: ",
            paste(missing_vars, collapse = ", "),
            call. = FALSE
        )
    }
    ids <- trimws(as.character(metadata[[1]]))
    if (anyDuplicated(ids)) {
        stop("Duplicated sample IDs in the first column of metadata: ",
            paste(head(unique(ids[duplicated(ids)])), collapse = ", "),
            call. = FALSE
        )
    }
    env <- metadata[, env_vars, drop = FALSE]
    rownames(env) <- ids
    env
}

.mbm_theme <- function(legend_position = "bottom", extra = NULL) {
    base <- .mbm_base_theme() + ggplot2::theme(legend.position = legend_position)
    if (!is.null(extra)) base <- base + extra
    base
}

.mbm_cache <- new.env(parent = emptyenv())
.mbm_base_theme <- function() {
    if (is.null(.mbm_cache$theme)) .mbm_cache$theme <- .mbm_build_theme()
    .mbm_cache$theme
}
.mbm_build_theme <- function() {
    ggplot2::theme_bw(base_family = "serif") +
        ggplot2::theme(
            plot.title = ggplot2::element_text(
                hjust = 0.5, face = "bold",
                size = 14, color = "black"
            ),
            plot.subtitle = ggplot2::element_text(
                hjust = 0.5, size = 12,
                color = "black"
            ),
            axis.text.x = ggplot2::element_text(size = 12, color = "black"),
            axis.text.y = ggplot2::element_text(size = 12, color = "black"),
            axis.title.x = ggplot2::element_text(size = 14, color = "black"),
            axis.title.y = ggplot2::element_text(size = 14, color = "black"),
            legend.title = ggplot2::element_text(
                size = 14, face = "bold",
                color = "black"
            ),
            legend.text = ggplot2::element_text(size = 12, color = "black"),
            strip.text = ggplot2::element_text(
                size = 12, face = "bold",
                color = "black"
            ),
            strip.background = ggplot2::element_rect(
                fill = "white",
                color = "black"
            ),
            panel.grid.minor = ggplot2::element_blank()
        )
}

.mbm_x_text <- function(angle = 0, colour = "black", size = 12) {
    ggplot2::element_text(
        angle = angle,
        hjust = if (angle == 0) 0.5 else 1,
        vjust = if (angle == 0) 0.5 else 1,
        size  = size,
        color = colour
    )
}

.mbm_strip_text <- function(bold = FALSE, size = 12, angle = NULL, colour = "black") {
    ggplot2::element_text(
        size   = size,
        face   = if (bold) "bold" else "plain",
        family = "serif",
        color  = colour,
        angle  = angle
    )
}

.mbm_div_palettes <- list(
    "BuOr" = c("#0072B2", "white", "#E69F00"),
    "BuVm" = c("#0072B2", "white", "#D55E00"),
    "BuPk" = c("#0072B2", "white", "#CC79A7"),
    "GnPk" = c("#009E73", "white", "#CC79A7"),
    "PuYl" = c("#440154", "white", "#FDE725")
)

.mbm_q_labels <- c(
    q0 = 'bolditalic(q)*bold("=0")',
    q1 = 'bolditalic(q)*bold("=1")',
    q2 = 'bolditalic(q)*bold("=2")'
)

.mbm_q_labels_plain <- c(
    q0 = 'italic(q)*"=0"',
    q1 = 'italic(q)*"=1"',
    q2 = 'italic(q)*"=2"'
)

.mbm_q_labeller <- function(bold = FALSE) {
    ggplot2::as_labeller(
        if (bold) .mbm_q_labels else .mbm_q_labels_plain,
        default = ggplot2::label_parsed
    )
}

.mbm_normalize_pair <- function(x) {
    parts <- strsplit(x, "_vs_", fixed = TRUE)
    vapply(parts, function(p) paste0(pmin(p[1], p[2]), "_vs_", pmax(p[1], p[2])), character(1))
}

.mbm_colors <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
)

.mbm_colors_2group <- c("#E69F00", "#56B4E9")

.mbm_colors_safe <- c(
    "#CC6677", "#332288", "#999933", "#882255",
    "#44AA99", "#AA4499", "#661100", "#88CCEE",
    "#DDCC77", "#117733", "#6699CC", "#888888"
)

.mbm_match_colors <- function(colors, values, arg = "colors") {
    if (is.null(names(colors))) {
        return(stats::setNames(colors[seq_along(values)], values))
    }
    .mbm_check_color_names(colors, values, arg)
    out <- stats::setNames(unname(colors[values]), values)
    out[is.na(out)] <- "grey70"
    out
}

.mbm_composite_genus <- function(x) {
    gsub("(?<![A-Za-z-])(?:[A-Z][a-z]+-){2,}([A-Z][a-z]+)(?![A-Za-z-])",
        "\\1 group", x,
        perl = TRUE
    )
}

.mbm_shorten_labels <- function(x, max_length = 35, composite = TRUE) {
    if (composite) x <- .mbm_composite_genus(x)
    if (!is.null(max_length)) {
        long <- nchar(x) > max_length
        x[long] <- paste0(substr(x[long], 1, max_length - 1), "\u2026")
    }
    make.unique(x, sep = "_")
}

utils::globalVariables(c(
    "%>%",
    ".",
    ".data",
    ":=",
    "taxonomy",
    "taxonomy2",
    "taxonomy_original",
    "OTUID",
    "OTU_ID",
    "SampleID",
    "Taxon",
    "Phylum",
    "taxon",
    "taxRank",
    "Abundance",
    "RelativeAbundance",
    "MeanAbundance",
    "SumAbundance",
    "MeanAbund",
    "max_abund",
    "abun",
    "abund",
    "PC1",
    "PC2",
    "x",
    "y",
    "cntr.x",
    "cntr.y",
    "label",
    "ids",
    "TD_beta",
    "site1",
    "site2",
    "value",
    "name",
    "comparison",
    "grupo",
    "diff.btw",
    "effect",
    "wi.eBH",
    "wi.ep",
    "log_pvalue",
    ".p",
    ".tag",
    "significant",
    "after_stat",
    "p",
    "direction",
    "lfc",
    "direct",
    "Ratio",
    "Dominant",
    "Condition",
    "MeanDecreaseGini",
    "Correlation",
    "Environmental",
    "compar_condition1",
    "compar_condition2",
    "k", "p", "o", "f", "g", "s",
    "asv",
    "phylum",
    "Group",
    "condition1_group",
    "across",
    "all_of",
    "where",
    "everything",
    "desc",
    "vars",
    "unit",
    "setNames",
    "reorder",
    "head",
    "as.dist",
    "as.formula",
    "cmdscale",
    "cor.test",
    "dist",
    "prcomp",
    "sd",
    "var",
    ".grp_idx",
    "side",
    "slot",
    "hill",
    "x_pos",
    "y_pos",
    "hjust",
    "vjust",
    "geo_km",
    "similarity",
    "group",
    "Index",
    "Variable",
    "row_id"
))

.mbm_hill_pairwise <- function(comm, q) {
    comm <- as.matrix(comm)
    sites <- rownames(comm)
    n <- nrow(comm)
    P <- comm / rowSums(comm)

    if (q == 0) {
        a_site <- rowSums(P > 0)
    } else if (q == 1) {
        a_site <- rowSums((P / 2) * log(P / 2), na.rm = TRUE)
    } else {
        a_site <- rowSums((P / 2)^q)
    }

    idx <- which(upper.tri(diag(n)), arr.ind = TRUE)
    idx <- idx[order(idx[, "col"], idx[, "row"]), , drop = FALSE]
    i <- idx[, "row"]
    j <- idx[, "col"]

    gamma <- numeric(nrow(idx))
    for (s1 in unique(i)) {
        k <- which(i == s1)
        G <- (P[s1, ] + t(P[j[k], , drop = FALSE])) / 2
        G <- sweep(G, 2, colSums(G), "/")
        gamma[k] <- if (q == 0) {
            colSums(G > 0)
        } else if (q == 1) {
            exp(colSums(-G * log(G), na.rm = TRUE))
        } else {
            colSums(G^q)^(1 / (1 - q))
        }
    }

    alpha <- if (q == 0) {
        (a_site[i] + a_site[j]) / 2
    } else if (q == 1) {
        exp(-(a_site[i] + a_site[j])) / 2
    } else {
        (1 / 2) * (a_site[i] + a_site[j])^(1 / (1 - q))
    }
    beta <- gamma / alpha
    if (q == 1) {
        local <- (log(2) - log(gamma) + log(alpha)) / log(2)
        region <- local
    } else {
        local <- (2^(1 - q) - beta^(1 - q)) / (2^(1 - q) - 1)
        region <- ((1 / beta)^(1 - q) - (1 / 2)^(1 - q)) / (1 - (1 / 2)^(1 - q))
    }
    tibble::tibble(
        q = q, site1 = sites[i], site2 = sites[j],
        TD_gamma = gamma, TD_alpha = unname(alpha), TD_beta = unname(beta),
        local_similarity = unname(local), region_similarity = unname(region)
    )
}

.mbm_from_parts <- function(counts, tax, samples) {
    table <- as.data.frame(counts, check.names = FALSE)
    table$taxonomy <- .mbm_taxonomy_string(tax, rownames(counts))

    if (is.null(samples)) {
        metadata <- data.frame(SAMPLEID = colnames(counts))
    } else {
        samples <- samples[colnames(counts), , drop = FALSE]
        samples$SAMPLEID <- NULL
        metadata <- data.frame(
            SAMPLEID = colnames(counts), samples,
            check.names = FALSE, row.names = NULL
        )
    }
    list(table = table, metadata = metadata)
}

.mbm_taxonomy_string <- function(tax, feature_ids) {
    if (is.null(tax)) {
        warning("No taxonomy found: the `taxonomy` column is left empty.", call. = FALSE)
        return(rep("", length(feature_ids)))
    }
    tax <- as.data.frame(tax, check.names = FALSE)
    tax_col <- grep("^taxonomy$", names(tax), ignore.case = TRUE)
    if (length(tax_col) == 1) {
        return(as.character(tax[[tax_col]]))
    }

    ranks <- c(
        "kingdom|domain|superkingdom", "phylum", "class", "order",
        "family", "genus", "species"
    )
    prefixes <- c("k", "p", "c", "o", "f", "g", "s")
    cols <- vapply(ranks, function(r) {
        hit <- grep(paste0("^(", r, ")$"), names(tax), ignore.case = TRUE)
        if (length(hit) == 0) NA_integer_ else hit[1]
    }, integer(1))
    if (all(is.na(cols))) {
        stop("No taxonomic rank columns (Kingdom/Domain, Phylum, ..., Species) ",
            "found in the taxonomy.",
            call. = FALSE
        )
    }
    prefixes <- prefixes[!is.na(cols)]
    tax <- tax[, cols[!is.na(cols)], drop = FALSE]

    vapply(seq_len(nrow(tax)), function(i) {
        values <- trimws(sub("^([kdpcofgs]__|D_[0-9]+__)", "", as.character(unlist(tax[i, ]))))
        missing <- is.na(values) | values == ""
        n <- if (any(missing)) which(missing)[1] - 1 else length(values)
        if (n == 0) {
            return("Unassigned")
        }
        paste0(prefixes[seq_len(n)], "__", values[seq_len(n)], collapse = "; ")
    }, character(1))
}

.mbm_taxonomy_db <- function(taxonomy_db) {
    switch(tolower(taxonomy_db),
        "silva" = "silva",
        "gg" = ,
        "gg2" = ,
        "greengenes2" = "gg2",
        "unite" = "unite",
        "kraken" = ,
        "kraken2" = "Kraken2",
        stop("Invalid `taxonomy_db`: '", taxonomy_db,
            "'. Choose one of: \"silva\", \"gg2\", \"unite\", \"Kraken2\" ",
            "(case-insensitive).",
            call. = FALSE
        )
    )
}

.mbm_check_level <- function(level) {
    level <- tolower(level)
    valid <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    if (length(level) != 1 || !level %in% valid) {
        stop("Invalid `level`: '", paste(level, collapse = ", "), "'. Choose one of: ",
            paste0('"', valid, '"', collapse = ", "), " (case-insensitive).",
            call. = FALSE
        )
    }
    level
}

.mbm_truncate_taxonomy <- function(table, level, taxonomy_db) {
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
    table
}

.mbm_label_taxa <- function(df, level, taxonomy_db) {
    if (taxonomy_db %in% c("unite", "silva", "gg2") && level == "species") {
        df <- df %>%
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
        df <- df %>%
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

    if (taxonomy_db %in% c("silva") && level == "genus") {
        df <- df %>%
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
        df <- df %>%
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
    if (taxonomy_db %in% c("unite", "silva", "Kraken2", "gg2") && level == "family") {
        df <- df %>%
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

    if (taxonomy_db %in% c("unite", "silva", "Kraken2", "gg2") && level == "order") {
        df <- df %>%
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

    if (taxonomy_db %in% c("unite", "silva", "Kraken2", "gg2") && level == "class") {
        df <- df %>%
            dplyr::mutate(
                taxonomy = dplyr::case_when(
                    taxonomy == "Other" ~ "Other",
                    grepl("c__[^;]*", taxonomy) & !grepl("c__uncultured|c__$|c__Incertae_Sedis", taxonomy) ~ sub(".*c__([^;]*).*", "\\1", taxonomy),
                    grepl("p__[^;]*", taxonomy) & !grepl("p__uncultured|p__$|p__Incertae_Sedis", taxonomy) ~ paste0("other ", stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .)),
                    TRUE ~ "Unclassified"
                )
            )
    }
    if (taxonomy_db %in% c("unite", "silva", "Kraken2", "gg2") && level == "phylum") {
        df <- df %>%
            dplyr::mutate(
                taxonomy = dplyr::case_when(
                    taxonomy == "Other" ~ "Other",
                    grepl("p__[^;]*", taxonomy) ~ stringr::str_extract(taxonomy, "p__[^;]*") %>% sub("p__", "", .),
                    TRUE ~ "Unclassified"
                )
            )
    }
    df
}
