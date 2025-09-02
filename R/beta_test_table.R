
#Table of Permanova or Betadisper 

#This function create a table with the results of permanova or betadisper 

beta_test_table <- function(matriz,
                            metadata,
                            formula_str,
                            method = "euclidean", 
                            test = c("permanova", "betadisper"),
                            permutations = 999,
                            strata_var = NULL,
                            decimales = 3) {
  
  test <- match.arg(test)
  
  # --- Verificaciones ---
  if (test == "permanova") {
    vars <- all.vars(as.formula(paste("~", formula_str)))
    vars_in_metadata <- vars %in% colnames(metadata)
    if (!all(vars_in_metadata)) {
      stop(paste("Las variables", paste(vars[!vars_in_metadata], collapse=", "), "no están en metadata"))
    }
  }
  
  strata <- NULL
  if (!is.null(strata_var)) {
    if (!strata_var %in% colnames(metadata)) {
      stop(paste("La variable strata", strata_var, "no está en metadata"))
    }
    strata <- metadata[[strata_var]]
  }
  
  
  # --- PERMANOVA ---
  if (test == "permanova") {
    dist_matrix <- vegan::vegdist(matriz, method = method) 
    resultado <- vegan::adonis2(as.formula(paste("dist_matrix ~", formula_str)),
                                data = metadata,
                                method = method,
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
    dist_matrix <- vegan::vegdist(matriz, method = method)
    disp <- vegan::betadisper(dist_matrix, metadata[[var_group]])
    perm <- vegan::permutest(disp, permutations = permutations)
    
    tabla <- as.data.frame(perm$tab)
    tabla$Term <- rownames(tabla)
    rownames(tabla) <- NULL
  }
  
  
  # --- Formato numérico ---
  tabla <- tabla %>%
    dplyr::mutate(across(where(is.numeric), ~ round(., decimales)))
  
  # --- Crear ggtexttable ---
  tab <- ggpubr::ggtexttable(tabla, rows = NULL, theme = ggpubr::ttheme(colnames = ggpubr::colnames_style(fill = "gray", color = "black"), # encabezado gris
                                                                        tbody.style = ggpubr::tbody_style(fill = "white")                  # filas en blanco
  )
  )
  
  # --- Agregar línea gruesa bajo encabezado ---
  tab <- tab %>% ggpubr::tab_add_hline(at.row = 1, row.side = "top", linewidth = 4) %>%
    ggpubr::tab_add_hline(at.row = 2, row.side = "top", linewidth = 4)
  
  # --- Resaltar p-valores automáticamente ---
  col_p <- grep("Pr", names(tabla), ignore.case = TRUE)
  if (length(col_p) > 0) {
    # convertir a num para comparación
    p_values <- suppressWarnings(as.numeric(tabla[[col_p]]))
    # filas con p < 0.05
    filas_signif <- which(!is.na(p_values) & p_values < 0.05)
    if (length(filas_signif) > 0) {
      for (fila in filas_signif) {
        tab <- ggpubr::table_cell_font(tab, row = fila + 1, column = col_p, face = "bold")
      }
    }
  }
  
  return(tab)
}