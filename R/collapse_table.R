#' Title
#'
#' @param table 
#' @param metadata 
#' @param level 
#' @param rel_abun 
#' @param export_txt 
#' @param file_name 
#'
#' @return
#' @export
#'
#' @examples
#' collapse_table(
#               table,
#               metadata,
#               level = "genus",
#               rel_abun = FALSE,
#               export_txt = FALSE
#               )

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
  table_final <- column_to_rownames(table_final, "OTU_ID")
  
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
      cols = all_of(ordered_samples),
      names_to = "SAMPLEID",
      values_to = ifelse(rel_abun, "RelativeAbundance", "Counts")
    )
  
  # --- Exportar a TXT si se solicita ---
  if (export_txt) {
    write.table(
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
