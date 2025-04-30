#' Hill Number Plot
#'
#' This function generates a boxplot to visualize Hill numbers (q = 0, 1, 2) for a given dataset,
#' stratified by a specified categorical variable (e.g., sample type or treatment). It allows
#' customization of various plot aesthetics, such as color palettes, facet orientation, and legend visibility.
#'
#' @param table A data frame or matrix with samples as columns and taxa as rows. The values represent the
#'        abundance of each taxon in each sample.
#' @param metadata A data frame with metadata, where the first column should contain the sample names.
#' @param x_col The column name in `metadata` that defines the x-axis grouping variable.
#' @param fill_col The column name in `metadata` that defines the fill color for the boxplots.
#' @param palette_type A character string specifying the color palette type. Options are "grey", "viridis", and "brewer". Default is "grey".
#' @param custom_palette A vector of custom colors to use for the `fill_col` variable if `palette_type` is set to `NULL`.
#' @param facet_orientation A character string indicating the facet orientation. Options are "horizontal" (default) or "vertical".
#' @param n_cols The number of columns to display in the facet grid. If `NULL`, the function automatically decides based on `facet_orientation`.
#' @param strip_color The color of the strip background in the facets. Default is "grey".
#' @param show_legend Logical indicating whether to display the legend. Default is `TRUE`.
#'
#' @return A `ggplot` object representing the boxplot of Hill numbers for each q-value, stratified by `x_col` and `fill_col`.
#' @export
#'
#' @examples a
hill_plot <-
    function(table,
    metadata,
    x_col,
    fill_col,
    palette_type = "grey",
    custom_palette = NULL,
    facet_orientation = "horizontal",
    n_cols = NULL,
    strip_color = "grey",
    show_legend = TRUE) {
        library(hillR)
        library(ggplot2)
        library(tidyr)
        library(patchwork)
        library(ggpubr)


        # Extraer los nombres de las muestras desde la primera columna de metadata
        sample_order <- metadata[[1]]

        # Verificar si los nombres de las columnas en table coinciden con el orden de las muestras en metadata
        common_samples <- intersect(colnames(table), sample_order)

        if (length(common_samples) == 0) {
            stop("No matching sample names between table and metadata.")
        }

        # Filtrar table para asegurarnos de que solo tenemos las columnas correspondientes a las muestras en metadata
        table <- table[, common_samples, drop = FALSE]

        # Verificar si los nombres de las columnas están en el mismo orden que sample_order
        table <- table[, match(sample_order, colnames(table))]

        # Transponer la table
        table <- data.frame(t(table))

        # Calcular los valores para q= 0,1 y 2 y unirlos en un dataframe con los metadatos
        resultado <- data.frame(
            q0 = hill_taxa(comm = table, q = 0),
            q1 = hill_taxa(comm = table, q = 1),
            q2 = hill_taxa(comm = table, q = 2),
            metadata
        )

        # Asegurar que los factores x_col y fill_col sean correctos
        resultado[[fill_col]] <-
            factor(resultado[[fill_col]], levels = unique(resultado[[fill_col]]))
        resultado[[x_col]] <-
            factor(resultado[[x_col]], levels = unique(resultado[[x_col]]))

        # Reestructurar el dataframe a formato largo
        resultado_largo <- resultado %>%
            pivot_longer(
                cols = starts_with("q"),
                names_to = "q",
                values_to = "value"
            )

        # Seleccionar la paleta de color según el argumento
        fill_scale <- if (!is.null(custom_palette)) {
            scale_fill_manual(values = custom_palette)
        } else {
            switch(palette_type,
                "grey" = scale_fill_grey(start = 0.9, end = 0.3),
                "viridis" = scale_fill_viridis_d(option = "plasma"),
                "brewer" = scale_fill_brewer(palette = "Set2"),
                scale_fill_grey(start = 0.9, end = 0.3) # Opción por defecto si no se reconoce el argumento
            )
        }

        # Determinar la configuración de facet_wrap()
        facet_wrap_config <- if (is.null(n_cols)) {
            if (facet_orientation == "vertical") {
                facet_wrap(~q, ncol = 1, scales = "free_y") # Una columna (vertical)
            } else {
                facet_wrap(~q, nrow = 1, scales = "free_y") # Una fila (horizontal, por defecto)
            }
        } else {
            facet_wrap(~q, ncol = n_cols, scales = "free_y") # Especificar número de columnas
        }

        # Gráfico con facet_wrap
        p <-
            ggplot(resultado_largo, aes(x = .data[[x_col]], y = value, fill = .data[[fill_col]])) +
            geom_boxplot(width = 0.4, position = position_dodge(width = 0.9)) +
            fill_scale + # Aplicar la paleta de colores seleccionada
            facet_wrap_config +
            labs(title = "", x = "", y = "Effective number of ASVs") +
            theme_bw() +
            theme(
                panel.grid = element_blank(),
                legend.position = if (show_legend) {
                    "right"
                } else {
                    "none"
                },
                panel.spacing = unit(1, "lines"),
                aspect.ratio = 0.4,
                strip.text = element_text(
                    face = "bold",
                    color = "black",
                    size = 15
                ),
                strip.background = element_rect(fill = strip_color)
            )

        return(p)
    }
