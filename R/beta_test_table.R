#' Table of Permanova or Betadisper 
#' 
#' This function create a table with the results of permanova or betadisper
#' 
#' 
#' @param table Distance matrix or data frame with taxonomy, where the columns are the samples and rows are ASVs or taxa.
#' @param metadata Data frame of characteristics or important information of the samples
#' @param formula_str Model formula
#' @param method Method for calculating pairwise distances. Same convention
#'   as \code{beta_div_plot}'s \code{distance} argument: \code{"compositional"}
#'   runs ALDEx2's CLR transform (\code{ALDEx2::aldex.clr}) on the raw counts
#'   and then a Euclidean distance on the CLR values (requires \code{table} to
#'   be a raw abundance data frame with a taxonomy column, not a precomputed
#'   distance matrix); \code{"aitchison"}/\code{"robust.aitchison"} use
#'   vegan's built-in Aitchison distance (\code{vegan::vegdist} with a
#'   pseudocount); any other value (e.g. \code{"euclidean"}, \code{"bray"})
#'   is passed straight to \code{vegan::vegdist} on the raw values.
#' @param test Permanova or betadisper 
#' @param permutations Number of permutations required
#' @param strata_var Group or variable within which permutations are restricted
#' @param decimales Number of decimales required
#' @param save_table Logical. If \code{TRUE}, saves the results table to disk.
#'   Default \code{FALSE}.
#' @param table_filename Character. File path/name for the saved table (used
#'   when \code{save_table = TRUE}). Default \code{"beta_test_results.txt"}.
#'
#' @return A table with the results of R2, F and p value 
#' @export
#'
#' @examples
#' \dontrun{
#' # Example using a data frame
#' beta_test_table(
#'   table       = table,
#'   metadata    = metadata,
#'   formula_str = "metodo*edad",
#'   method      = "euclidean",
#'   test        = "permanova",
#'   permutations = 999,
#'   strata_var  = "Individuo"
#' )
#'
#' # Example using a distance matrix
#' beta_test_table(
#'   table       = dist_matrix,
#'   metadata    = metadata,
#'   formula_str = "Origen",
#'   test        = "betadisper"
#' )
#'
#' # Compositional PERMANOVA (CLR via ALDEx2, then Euclidean)
#' beta_test_table(
#'   table       = table,
#'   metadata    = metadata,
#'   formula_str = "dist_km",
#'   method      = "compositional",
#'   test        = "permanova",
#'   permutations = 999
#' )
#' }
beta_test_table <- function(table,
                            metadata,
                            formula_str,
                            method = "euclidean", 
                            test = c("permanova", "betadisper"),
                            permutations = 999,
                            strata_var = NULL,
                            decimales = 3,
                            save_table = FALSE,
                            table_filename = "beta_test_results.txt") {
  
  test <- match.arg(test)
  raw_input <- is.data.frame(table)

  # --- Aceptar tambien data.frame como matriz ---
  if (is.data.frame(table)) {
    # Si la ultima columna parece taxonomia, eliminarla
    tax_cols <- grep("taxonomy|taxon|Taxonomy|Taxa", names(table))
    if (length(tax_cols) > 0) {
      table <- table[, -tax_cols[1], drop = FALSE]
    }
    
    # Convertir solo columnas numericas
    num_cols <- sapply(table, is.numeric)
    if (!all(num_cols)) {
    }
    table <- as.matrix(table[, num_cols, drop = FALSE])
  } else if (!is.matrix(table)) {
    stop("'table' must be a matrix or dataframe")
  }
  
  # --- Detectar orientacion ---
  if (ncol(table) == nrow(metadata)) {
    table <- t(table)
  } else if (nrow(table) != nrow(metadata)) {
    stop("Dimensons of table and metadata don't match.")
  }
  
  # --- Verificaciones ---
  if (test == "permanova") {
    vars <- all.vars(as.formula(paste("~", formula_str)))
    vars_in_metadata <- vars %in% colnames(metadata)
    if (!all(vars_in_metadata)) {
      stop(paste("Variables", paste(vars[!vars_in_metadata], collapse=", "), "are not in metadata"))
    }
  }
  
  strata <- NULL
  if (!is.null(strata_var)) {
    if (!strata_var %in% colnames(metadata)) {
      stop(paste("Strata variable", strata_var, "is not in metadata"))
    }
    strata <- metadata[[strata_var]]
  }
  
  # --- Distancia (mismo criterio que beta_div_plot) ---
  if (method == "compositional") {
    if (!raw_input) {
      stop("method = 'compositional' requires a raw abundance table ",
           "(data frame with a taxonomy column), not a precomputed distance matrix.")
    }
    set.seed(123)
    aldex_obj   <- ALDEx2::aldex.clr(t(table), mc.samples = 128,
                                     denom = "all", verbose = FALSE, useMC = FALSE)
    clr_samples <- t(ALDEx2::getMonteCarloSample(aldex_obj, 1))
    dist_matrix <- stats::dist(clr_samples, method = "euclidean")
  } else if (method %in% c("aitchison", "robust.aitchison")) {
    dist_matrix <- vegan::vegdist(table, method = method, pseudocount = 0.5)
  } else {
    dist_matrix <- vegan::vegdist(table, method = method)
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
    tabla$Term <- rownames(tabla)
    rownames(tabla) <- NULL
  }
  
  # --- Formato numerico ---
  tabla <- tabla %>%
    dplyr::mutate(across(where(is.numeric), ~ round(., decimales))) %>%
    dplyr::mutate(across(everything(), as.character)) %>%
    dplyr::mutate(across(everything(), ~ ifelse(is.na(.), "-", .)))

  if (save_table) {
    utils::write.table(tabla, file = table_filename, sep = "\t",
                       quote = FALSE, row.names = FALSE)
    message(paste("Table saved as:", table_filename))
  }

  # --- Crear tabla visual ---
  tab <- ggpubr::ggtexttable(tabla,
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
    p_values <- suppressWarnings(as.numeric(tabla[[col_p]]))
    filas_signif <- which(!is.na(p_values) & p_values < 0.05)
    if (length(filas_signif) > 0) {
      for (fila in filas_signif) {
        tab <- ggpubr::table_cell_font(tab, row = fila + 1, column = col_p, face = "bold")
      }
    }
  }
  

  return(tab)
}

