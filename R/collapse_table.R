#' Collapse table
#' Collapse an OTU/ASV table by taxonomic level
#'
#' Collapses an abundance table by a specified taxonomic rank (e.g. genus,
#' family, phylum), summing counts across features that share the same taxonomy.
#' Features with lower taxonomic resolution than the selected level are retained
#' unchanged. Optionally converts counts to relative abundance and exports the
#' collapsed table to a tab-delimited file.
#'
#' @param table A data frame containing an OTU/ASV abundance table. Must include
#'   a column with OTUID, one with taxonomy and sample columns with numeric counts.
#' @param metadata A data frame containing sample metadata. The first column is
#'   used to define the order of samples in the output table.
#' @param level Character string specifying the taxonomic level to collapse to.
#'   ("kingdom","phylum","class","order","family","genus" or"species".
#' @param rel_abun Logical. If TRUE, converts counts to relative abundance
#'   (%) per sample.
#' @param export_txt Logical. If TRUE, exports the collapsed table as a
#'   tab-delimited text file.
#' @param file_name Character. Name of the output file when
#'   export_txt = TRUE.
#' @return a table collapsed by taxonomic level
#' @export
#' @examples
#' \dontrun{
#' collapse_table(
#'   table      = table,
#'   metadata   = metadata,
#'   level      = "genus",
#'   rel_abun   = FALSE,
#'   export_txt = FALSE
#' )
#' }

collapse_table <- function(table,
                          metadata,
                          level = "genus",
                          rel_abun = FALSE,
                          export_txt = FALSE,
                          file_name = "collapsed_table.txt") {
  # --- Ordenar columnas ---
  table <- table[, c("taxonomy", setdiff(names(table), "taxonomy"))]
  ordered_samples <- intersect(metadata[, 1], colnames(table)[-1])
  table <- table[, c("taxonomy", ordered_samples)]
  
  # Si tiene rownames, guardarlos como columna temporal
  if (!is.null(rownames(table))) {
    table <- tibble::rownames_to_column(table, var = "OTU_ID")
  } else {
    table$OTU_ID <- paste0("OTU_", seq_len(nrow(table)))
  }
  
  # --- Limpiar terminaciones vacías ---
  table$taxonomy <- gsub("(;__)+$", "", table$taxonomy)
  
  # --- Mapear nivel deseado a índice ---
  level_idx <- switch(
    level,
    kingdom = 1,
    phylum = 2,
    class = 3,
    order = 4,
    family = 5,
    genus = 6,
    species = 7
  )
  
  # --- Calcular profundidad de cada fila ---
  get_depth <- function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    sum(grepl("__", levels))
  }
  table$depth <- sapply(table$taxonomy, get_depth)
  
  # --- Separar filas según profundidad ---
  lowres <-
    table[table$depth < level_idx,]    # menos resueltas, se conservan
  highres <-
    table[table$depth >= level_idx,]  # suficientemente resueltas para colapsar
  
  # --- Colapsar highres hasta el nivel deseado ---
  highres$taxonomy <- sapply(highres$taxonomy, function(tax) {
    levels <- unlist(strsplit(tax, ";"))
    paste(levels[1:level_idx], collapse = ";")
  })
  
  # --- Agrupar y sumar, conservando un ID original ---
  highres <- highres %>%
    dplyr::group_by(taxonomy) %>%
    dplyr::summarise(
      OTU_ID = dplyr::first(OTU_ID),
      # conservar uno de los IDs originales
      dplyr::across(where(is.numeric), sum, na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- Combinar lowres y highres ---
  table_final <- dplyr::bind_rows(lowres[, c("OTU_ID", "taxonomy", ordered_samples)],
                                  highres[, c("OTU_ID", "taxonomy", ordered_samples)])
  
  # --- Eliminar columna depth ---
  table_final <-
    table_final[, c("OTU_ID", "taxonomy", ordered_samples)]
  table_final <- tibble::column_to_rownames(table_final, "OTU_ID")
  
  # --- Convertir a abundancia relativa ---
  if (rel_abun) {
    table_final[, ordered_samples] <- sweep(table_final[, ordered_samples],
                                            2,
                                            colSums(table_final[, ordered_samples], na.rm = TRUE),
                                            FUN = "/") * 100
  }
  
  # --- Formato largo ---
  table_long <- table_final %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(ordered_samples),
      names_to = "SAMPLEID",
      values_to = ifelse(rel_abun, "RelativeAbundance", "Counts")
    )
  
  # --- Exportar a TXT si se solicita ---
  if (export_txt) {
    utils::write.table(
      table_final,
      file = file_name,
      sep = "\t",
      quote = FALSE,
      col.names = NA
    )
  }
  
  # --- Asignar al ambiente global ---
  assign("collapsed_table", table_final, envir = .GlobalEnv)
  assign("collapsed_long_table", table_long, envir = .GlobalEnv)
  
  return(list(collapsed_table = table_final,
              long_format = table_long))
}
