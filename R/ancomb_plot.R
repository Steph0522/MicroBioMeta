#' barplot ancombc
#'
#' This function generates a heatmap to visualize alpha diversity Hill numbers (q = 0, 1, 2)
#' for a given dataset, faceted by one or two categorical variables (e.g., sample type or treatment).
#' It supports palette customization, faceting, and statistical comparison.
#'
#' @param table, Data frame with taxonomy, where, the columns are the samples and rows are ASV's or taxa.
#' @param conditions, Vector that defines the categories or classes to compare.
#' @param effect, Effect size (default: >= 0.8)
#' @param pvalue_BH, Value p-ajusted (optional)
#' @param formula, Formula of the model 
#'
#' @return A plot with the deferentially abundant taxonomic groups between two categories or groups of samples.
#' @export
#'
#' @examples ancomb_plot(table = table,
#'                      conditions = conditions,
#'                      effect = 0.8,
#'                      formula = formula, 
#'                      pvalue_BH = NULL)
#'

ancombc_plot <- function(table,
                         metadata,
                         col_cond,
                         prv_cut = 0.1,
                         formula = NULL,
                         rand_formula = NULL) {
  #check ANCOMBC package
  
  if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
    message("El paquete 'ANCOMBC' is not installed. Installing from Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("ANCOMBC")
  }
  
#check condition in metadata
    if (!col_cond %in% colnames(metadata)) {
    stop(paste("Column", col_cond, "not found in metadata."))
  }
  
  
  tax_col <- grep("taxonomy|Taxonomy|taxon|Taxa|taxa|Taxon", names(table), ignore.case = TRUE)
  if(length(tax_col) != 1) stop("There is no taxonomy column in the table")
  
  names(table)[ncol(table)] <- "taxonomy"
  
 # counts table
  table <- table 
  table_counts <- table %>%
    dplyr::select(-taxonomy) 
  
  
  # create phyloseq object 
  
  otumat= table_counts %>% as.matrix()
  taxa= table %>% dplyr::select(Taxon=taxonomy) %>% rownames_to_column(var = "Feature.ID")
  taxmat = qiime2R::parse_taxonomy(taxa) %>% as.matrix()
  library("phyloseq")
  OTU = otu_table(otumat, taxa_are_rows = TRUE)
  TAX = tax_table(taxmat)
  sampledata = sample_data(metadata %>% column_to_rownames(var = "SAMPLEID"))
  physeq = phyloseq(OTU, TAX, sampledata )
  physeq_filt <- prune_taxa(apply(otu_table(physeq), 1, var) > 0, physeq)
  dat <- mia::makeTreeSummarizedExperimentFromPhyloseq(physeq_filt)



  # run ancombc
  ancomb_results <- ANCOMBC::ancombc2(
    data = dat, assay_name = "counts",
    rank = "Genus",
    fix_formula= formula, p_adj_method = "holm",
    pseudo_sens = TRUE,
    prv_cut = prv_cut,lib_cut = 1000,s0_perc = 0.05,
    group = col_cond,
    struc_zero = TRUE,neg_lb = TRUE)
  
  res_prim = ancomb_results$res
  
  
 
  conditions <- unique(metadata[[col_cond]])


  # plot if has 2 conditions

  df_cond = res_prim %>%
    dplyr::select(taxon, contains(col_cond)) 
  pattern <- paste0("^", col_cond, ".*")
  colnames(res_prim) <- gsub(pattern, col_cond, colnames(res_prim))
  lfc_col <- grep("^lfc_", names(df_cond), value = TRUE)[1]
  diff_col <- grep("^diff_", names(df_cond), value = TRUE)[1]
  se_col <- grep("^se_", names(df_cond), value = TRUE)[1]


  
  df_condition <- df_cond %>%
    filter(if_any(matches(diff_col), ~ . == TRUE)) %>%
    arrange(desc(.data[[lfc_col]])) %>%
    mutate(direct = ifelse(.data[[lfc_col]] > 0,
                           "Positive LFC", "Negative LFC"))
  
  df_condition$taxon = factor( df_condition$taxon, levels =  df_condition$taxon)
  df_condition$direct = factor(df_condition$direct, 
                             levels = c("Positive LFC", "Negative LFC"))
  cond2 <- unique(metadata[[col_cond]])[2]
  
  fig = df_condition %>%
    ggplot(aes(x = taxon, y = .data[[lfc_col]], fill = direct)) + 
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = .data[[lfc_col]] - .data[[se_col]], 
                      ymax = .data[[lfc_col]] + .data[[se_col]]), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = paste("Log fold changes as one unit increase", "in",conditions )) + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 1,
                                     color = df_condition$color))
  
  #plot it has more than 2 conditions
  
  df_cond = res_prim %>%
    dplyr::select(taxon, contains(col_cond)) 
  conditions <- unique(metadata[[col_cond]])

  diff_cols <- paste0("diff_", col_cond,conditions[-1])
  lfc_cols  <- paste0("lfc_", col_cond,conditions[-1])
  
  dfs_list <- list()
  
  for(i in seq_along(conditions[-1])) {
    df_temp <- df_cond %>%
      filter(.data[[diff_cols[i]]] == TRUE) %>%
      mutate(
        value = round(.data[[lfc_cols[i]]], 2)
      ) %>%
      select(taxon, value)
    
    # Guardar en la lista con nombre dinámico
    dfs_list[[conditions[i]]] <- df_temp
  }
  
  df_fig1 <- df_cond %>%
    dplyr::filter(diff_cols[[3]] == TRUE | 
                    diff_cols[[3]] == TRUE)  %>%
    mutate(
      lfc1 = ifelse(.data[[diff_cols[2]]] == TRUE, 
                    round(.data[[lfc_cols[2]]], 2), 0),
      lfc2 = ifelse(data[[diff_cols[3]]] == TRUE, 
                    round(.data[[lfc_cols[3]]], 2), 0)
    ) %>%
    pivot_longer(cols = lfc1:lfc2, names_to = "group", values_to = "value") %>%
    arrange(taxon)
  
  
  df_fig2 <- df_cond %>%
    filter(if_any(all_of(diff_cols), ~ . == TRUE)) %>%   
    mutate(
      lfc1 = ifelse(.data[[lfc_cols[1]]] == 1, round(.data[[lfc_cols[1]]], 2), 0),
      lfc2 = ifelse(.data[[lfc_cols[3]]] == 1, round(.data[[lfc_cols[2]]], 2), 0)
    ) %>%
    pivot_longer(cols = lfc1:lfc2, names_to = "group", values_to = "value") %>%
    arrange(taxon)

  df_fig_bmi = df_fig_bmi1 %>%
    dplyr::left_join(df_fig_bmi2, by = c("taxon", "group"))
  
  df_fig_bmi$group = recode(df_fig_bmi$group, 
                            `lfc1` = "Overweight - Obese",
                            `lfc2` = "Lean - Obese")
  df_fig_bmi$group = factor(df_fig_bmi$group, 
                            levels = c("Overweight - Obese",
                                       "Lean - Obese"))
  
  lo = floor(min(df_fig_bmi$value))
  up = ceiling(max(df_fig_bmi$value))
  mid = (lo + up)/2
  fig_bmi = df_fig_bmi %>%
    ggplot(aes(x = group, y = taxon, fill = value)) + 
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         na.value = "white", midpoint = mid, limit = c(lo, up),
                         name = NULL) +
    geom_text(aes(group, taxon, label = value, color = color), size = 4) +
    scale_color_identity(guide = "none") +
    labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  fig_bmi
  
  
  return(invisible(fig))
}
