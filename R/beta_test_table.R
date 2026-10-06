#' Table of Permanova or Betadisper 
#' 
#' This function create a table with the results of permanova or betadisper
#' 
#' 
#' @param table A precomputed distance matrix (\code{matrix} or \code{dist}
#'   object, e.g. from \code{vegan::vegdist}), or a data frame with taxonomy,
#'   where the columns are the samples and rows are ASVs or taxa.
#' @param metadata Data frame of characteristics or important information of the samples
#' @param formula_str Model formula
#' @param distance Method for calculating pairwise distances. Same convention
#'   as \code{beta_div_plot}'s \code{distance} argument: \code{"compositional"}
#'   runs ALDEx2's CLR transform (\code{ALDEx2::aldex.clr}) on the raw counts
#'   and then a Euclidean distance on the CLR values (requires \code{table} to
#'   be a raw abundance data frame with a taxonomy column, not a precomputed
#'   distance matrix); \code{"aitchison"}/\code{"robust.aitchison"} use
#'   vegan's built-in Aitchison distance (\code{vegan::vegdist} with a
#'   pseudocount); any other value (e.g. \code{"euclidean"}, \code{"bray"})
#'   is passed straight to \code{vegan::vegdist} on the raw values
#'   (\code{"euclidean"} default; case-insensitive).
#' @param test Statistical test to run: one of \code{"permanova"} (default) or
#'   \code{"betadisper"}. Case-insensitive.
#' @param permutations Number of permutations required
#' @param strata_var Group or variable within which permutations are restricted
#' @param digits Number of decimal places for the numeric columns (except the p-value).
#' @param save_table Logical. If \code{TRUE}, saves the results table to disk.
#'   Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"beta_test_results.txt"}.
#' @param mc_samples Number of ALDEx2 Monte Carlo instances used when
#'   \code{distance = "compositional"}. With \code{1} (default) the clr values
#'   of one random instance are used: fast, but the result changes a little
#'   between runs (use \code{set.seed()}). With more, the clr values are
#'   averaged across instances, which gives an almost identical result in every
#'   run; \code{128} (ALDEx2's default) is suggested for final analyses, and
#'   takes longer. Ignored for other distances.
#'
#' @details The first column of \code{metadata} must hold the sample IDs;
#'   metadata rows are matched to the samples by ID, so their order doesn't
#'   matter. A precomputed distance is used as-is (\code{distance} is ignored).
#'   With \code{distance = "compositional"} and \code{mc_samples = 1}, the
#'   clr values come from one random Monte Carlo instance of
#'   \code{ALDEx2::aldex.clr()}; call \code{set.seed()} before the function
#'   to make the result reproducible, or use \code{mc_samples = 128}.
#'
#' @return A data frame (class \code{mbm_test_table}) with the test results:
#'   one row per term and the columns returned by \code{vegan::adonis2()}
#'   (\code{Df}, \code{SumOfSqs}, \code{R2}, \code{F}, \code{Pr(>F)}) or
#'   \code{vegan::permutest()}, plus \code{Term}. Printing it (e.g. typing its
#'   name) draws the formatted table figure; \code{ggplot2::autoplot()}
#'   returns that figure as a \code{ggplot} object, e.g. to combine it with
#'   other plots (\code{cowplot::plot_grid()}, \code{patchwork}) or save it
#'   with \code{ggplot2::ggsave()}.
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
#' # Example using a data frame
#' beta_test_table(
#'   table       = table,
#'   metadata    = metadata,
#'   formula_str = "Location*Treatment",
#'   distance    = "bray",
#'   test        = "permanova",
#'   permutations = 999,
#'   strata_var  = "Plot"
#' )
#'
#' # Example using a distance matrix
#' dist_matrix <- vegan::vegdist(
#'   t(table[, setdiff(colnames(table), "taxonomy")]),
#'   method = "bray"
#' )
#' beta_test_table(
#'   table       = dist_matrix,
#'   metadata    = metadata,
#'   formula_str = "Location",
#'   test        = "betadisper"
#' )
#'
#' # Compositional PERMANOVA (CLR via ALDEx2, then Euclidean)
#' beta_test_table(
#'   table       = table,
#'   metadata    = metadata,
#'   formula_str = "Location",
#'   distance    = "compositional",
#'   test        = "permanova",
#'   permutations = 999
#' )
beta_test_table <- function(table,
                            metadata,
                            formula_str,
                            distance = "euclidean", 
                            test = c("permanova", "betadisper"),
                            permutations = 999,
                            mc_samples = 1,
                            strata_var = NULL,
                            digits = 3,
                            save_table = FALSE,
                            table_filename = "beta_test_results.txt") {


  # Accept `test` and `distance` case-insensitively (distance names are matched
  # case-sensitively by vegan::vegdist, so normalising here avoids a cryptic
  # "invalid distance method" error from e.g. distance = "Bray").
  test   <- match.arg(tolower(test), c("permanova", "betadisper"))
  distance <- tolower(distance)
  # A precomputed distance (a dist object, or a square symmetric matrix with
  # a zero diagonal) is used as-is: no distance is computed and `distance` is
  # ignored, instead of running vegdist() on the distances themselves.
  is_square_dist <- is.matrix(table) && nrow(table) == ncol(table) &&
    isSymmetric(unname(table)) && all(diag(table) == 0)
  precomputed <- inherits(table, "dist") || is_square_dist
  raw_input <- is.data.frame(table)

  if (precomputed) {
    if (distance == "compositional") {
      stop("distance = 'compositional' requires a raw abundance table ",
           "(data frame with a taxonomy column), not a precomputed distance matrix.")
    }
    dist_matrix <- stats::as.dist(table)
    sample_ids  <- labels(dist_matrix)
  } else if (is.data.frame(table)) {
    # Si la ultima columna parece taxonomia, eliminarla
    tax_cols <- grep("taxonomy|taxon|Taxonomy|Taxa", names(table))
    if (length(tax_cols) > 0) {
      table <- table[, -tax_cols[1], drop = FALSE]
    }

    # Convertir solo columnas numericas
    num_cols <- vapply(table, is.numeric, logical(1))
    if (!all(num_cols)) {
    }
    table <- as.matrix(table[, num_cols, drop = FALSE])
  } else if (!is.matrix(table)) {
    stop("'table' must be a matrix, dataframe, or dist object")
  }

  # --- Detectar orientacion ---
  # Samples go in rows: transpose when the sample IDs are the column names.
  meta_ids <- trimws(as.character(metadata[[1]]))
  if (!precomputed) {
    n_cols_in_meta <- sum(colnames(table) %in% meta_ids)
    n_rows_in_meta <- sum(rownames(table) %in% meta_ids)
    if (n_cols_in_meta == 0 && n_rows_in_meta == 0) {
      stop("None of the sample names in table are in the first column of metadata.")
    }
    if (n_cols_in_meta >= n_rows_in_meta) table <- t(table)
    sample_ids <- rownames(table)
  }

  # --- Alinear metadata con las muestras (por ID, nunca por posicion) ---
  if (is.null(sample_ids)) {
    stop("table has no sample names; they are needed to match it with ",
         "the sample IDs in the first column of metadata.")
  }
  metadata <- .mbm_align_metadata(sample_ids, metadata)
  keep <- metadata[[1]]
  if (precomputed) {
    if (length(keep) < length(sample_ids)) {
      dist_matrix <- stats::as.dist(as.matrix(dist_matrix)[keep, keep])
    }
  } else {
    table <- table[keep, , drop = FALSE]
  }

  # --- Verificaciones ---
  if (test == "permanova") {
    vars <- all.vars(as.formula(paste("~", formula_str)))
    vars_in_metadata <- vars %in% colnames(metadata)
    if (!all(vars_in_metadata)) {
      stop("Variables ", paste(vars[!vars_in_metadata], collapse = ", "), " are not in metadata")
    }
  }
  
  strata <- NULL
  if (!is.null(strata_var)) {
    if (!strata_var %in% colnames(metadata)) {
      stop("Strata variable ", strata_var, " is not in metadata")
    }
    strata <- metadata[[strata_var]]
  }
  
  # --- Distancia (mismo criterio que beta_div_plot) ---
  if (precomputed) {
    # dist_matrix already set above
  } else if (distance == "compositional") {
    if (!raw_input) {
      stop("distance = 'compositional' requires a raw abundance table ",
           "(data frame with a taxonomy column), not a precomputed distance matrix.")
    }
    clr_samples <- .mbm_aldex_clr(t(table), mc_samples)
    dist_matrix <- stats::dist(clr_samples, method = "euclidean")
  } else if (distance %in% c("aitchison", "robust.aitchison")) {
    dist_matrix <- vegan::vegdist(table, method = distance, pseudocount = 0.5)
  } else {
    dist_matrix <- vegan::vegdist(table, method = distance)
  }

  # --- PERMANOVA ---
  if (test == "permanova") {
    resultado <- vegan::adonis2(as.formula(paste("dist_matrix ~", formula_str)),
                                data = metadata,
                                permutations = permutations,
                                strata = strata,
                                by = "terms")
    tabla <- as.data.frame(resultado)
    tabla$Term <- rownames(tabla)
    rownames(tabla) <- NULL
  }

  # --- BETADISPER / PERMDISP ---
  if (test == "betadisper") {
    var_group <- all.vars(as.formula(paste("~", formula_str)))[1]
    disp <- vegan::betadisper(dist_matrix, metadata[[var_group]])
    perm <- vegan::permutest(disp, permutations = permutations)
    
    tabla <- as.data.frame(perm$tab)
    # vegan::permutest() always labels the between-group row "Groups"; show
    # the actual variable name being tested instead, for clarity.
    tabla$Term <- ifelse(rownames(tabla) == "Groups", var_group, rownames(tabla))
    rownames(tabla) <- NULL
  }
  
  # --- Formato numerico ---
  # The p-value column uses significant-figure formatting (consistent with
  # every other figure in the package) instead of the fixed `digits`
  # rounding applied to the other numeric columns, since fixed decimals can
  # round small p-values (e.g. 0.0004) down to "0".
  results <- tabla
  col_p_name <- grep("Pr", names(tabla), ignore.case = TRUE, value = TRUE)
  tabla <- tabla %>%
    dplyr::mutate(across(where(is.numeric) & !dplyr::any_of(col_p_name), ~ round(., digits))) %>%
    dplyr::mutate(across(dplyr::any_of(col_p_name), ~ .mbm_format_pval(.))) %>%
    dplyr::mutate(across(everything(), as.character)) %>%
    dplyr::mutate(across(everything(), ~ ifelse(is.na(.), "-", .)))

  if (save_table) {
    utils::write.table(tabla, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message("Table saved as: ", table_filename)
  }

  # --- Crear tabla visual ---
  # vegan::adonis2()'s own column is literally named "R2"; only the
  # displayed table gets the proper R-squared (superscript 2) - the saved
  # table (save_table = TRUE) and the `tabla` data itself keep the plain
  # "R2" name for compatibility with scripts that read it back in.
  tabla_display <- tabla
  names(tabla_display)[names(tabla_display) == "R2"] <- "R\u00b2"

  tab <- ggpubr::ggtexttable(tabla_display,
                             rows = NULL,
                             theme = ggpubr::ttheme(
                               colnames = ggpubr::colnames_style(
                                 fill = "gray", 
                                 color = "black",
                                 face = "bold",
                                 size = 12,
                                 fontface = "plain",
                                 fontfamily = "serif"
                               ),
                               tbody.style = ggpubr::tbody_style(
                                 fill = "white",
                                 color = "black",
                                 size = 12,
                                 fontfamily = "serif"
                               )
                             )
  )
  
  # --- Lineas bajo encabezado ---
  tab <- tab %>%
    ggpubr::tab_add_hline(at.row = 1, row.side = "top", linewidth = 4) %>%
    ggpubr::tab_add_hline(at.row = 2, row.side = "top", linewidth = 4)
  
  # --- Resaltar p-valores ---
  col_p <- grep("Pr", names(tabla), ignore.case = TRUE)
  if (length(col_p) > 0) {
    p_text <- tabla[[col_p]]
    # Values formatted as "<0.001" are significant by definition but would
    # become NA under as.numeric(), so flag them separately.
    is_below_threshold <- startsWith(p_text, "<")
    p_values <- suppressWarnings(as.numeric(p_text))
    filas_signif <- which(is_below_threshold | (!is.na(p_values) & p_values < 0.05))
    if (length(filas_signif) > 0) {
      for (fila in filas_signif) {
        tab <- ggpubr::table_cell_font(tab, row = fila + 1, column = col_p, face = "bold")
      }
    }
  }
  
  class(results) <- c("mbm_test_table", class(results))
  attr(results, "plot") <- tab
  results
}

#' @export
print.mbm_test_table <- function(x, ...) {
  print(attr(x, "plot"))
  invisible(x)
}

# a subset of the table (e.g. x[1:2, ]) is a plain data frame, without the figure
#' @export
`[.mbm_test_table` <- function(x, ...) {
  out <- NextMethod()
  if (is.data.frame(out)) {
    attr(out, "plot") <- NULL
    class(out) <- setdiff(class(out), "mbm_test_table")
  }
  out
}

#' @exportS3Method ggplot2::autoplot
autoplot.mbm_test_table <- function(object, ...) {
  attr(object, "plot")
}

