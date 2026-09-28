# prototype_abundance_bar_plot.R
# abundance_bar_plot ya es rapido (~0.2 s construir + ~0.3 s renderizar; mas rapido que phyloseq y rbiom).
# Hotspots (Rprof): R/abundance_bar_plot.R#447 construccion ggplot (~19%), R/zzz.R#152 .mbm_theme/theme_bw
# (~15%), R/abundance_bar_plot.R#114-126 cadena de 12 dplyr::filter() (~12%), #194 group_by+summarise (~9%).
# Variante:
#   B1: un solo filtro vectorizado (!taxonomy %in% c(...)) en lugar de 12 dplyr::filter encadenados.
#   B2: B1 + colapso por taxonomia con rowsum() (base) en lugar de group_by/summarise(across()).
# Verificacion: p$data identico y render PNG pixel-identico.
# Uso: cd optimization/profiling; Rscript prototype_abundance_bar_plot.R
source("00_setup_data.R")
N <- as.integer(Sys.getenv("NREPS", "11"))

src_lines <- readLines(file.path(PKG_ROOT, "R", "abundance_bar_plot.R"))
i <- grep("table <- table %>%", src_lines, fixed = TRUE)[1]
j <- grep("dplyr::filter(taxonomy != \"d__Eukaryota\")", src_lines, fixed = TRUE)[1]
old_filter <- paste(src_lines[i:j], collapse = "\n")
new_filter <- paste(
  "  table <- table[!(table$taxonomy %in% c(\"d__Bacteria;__;__;__;__;__\", \"d__Bacteria\", \"d__Archaea;__;__;__;__;__\",",
  "    \"d__Archaea\", \"d__Bacteria;p__;c__;o__;f__;g__;s__\", \"d__Archaea;p__;c__;o__;f__;g__;s__\",",
  "    \"k__Bacteria;__;__;__;__;__\", \"k__Fungi;__;__;__;__;__\", \"k__Fungi;p__;c__;o__;f__;g__\", \"k__Fungi\",",
  "    \"Unassigned\", \"d__Eukaryota\")) & !is.na(table$taxonomy), , drop = FALSE]",
  sep = "\n")
old_collapse <- paste(
  "  table <- table %>%",
  "    dplyr::group_by(taxonomy) %>%",
  "    dplyr::summarise(dplyr::across(where(is.numeric), \\(x) sum(x, na.rm = TRUE)))", sep = "\n")
new_collapse <- paste(
  "  .num <- vapply(table, is.numeric, logical(1))",
  "  .m <- as.matrix(table[, .num, drop = FALSE]); .m[is.na(.m)] <- 0",
  "  .g <- factor(table$taxonomy, levels = sort(unique(table$taxonomy), method = \"radix\"))",
  "  .rs <- rowsum(.m, .g, reorder = TRUE)",
  "  table <- tibble::as_tibble(data.frame(taxonomy = rownames(.rs), .rs, check.names = FALSE, row.names = NULL))",
  sep = "\n")

f_B1 <- patch_fun("abundance_bar_plot.R", "abundance_bar_plot", list(c(old_filter, new_filter)))
f_B2 <- patch_fun("abundance_bar_plot.R", "abundance_bar_plot", list(c(old_filter, new_filter), c(old_collapse, new_collapse)))
args <- list(table = table_bac, metadata = metadata_bacteria, taxonomy_db = "silva", level = "phylum", top_n = 15,
             x_col = "Type_of_soil", facet_col = "Treatment", add_remained = TRUE, label = "Phylum",
             x_label_angle = 45, save_table = FALSE)
call_f <- function(f, ...) suppressWarnings(suppressMessages(do.call(f, utils::modifyList(args, list(...)))))

p0 <- call_f(abundance_bar_plot); p1 <- call_f(f_B1); p2 <- call_f(f_B2)
px0 <- render_png_pixels(p0)
cat("B1: p$data identico:", identical(p0$data, p1$data), "| PNG identico:", identical(px0, render_png_pixels(p1)), "\n")
cat("B2: p$data identico:", identical(p0$data, p2$data), "| all.equal:", isTRUE(all.equal(p0$data, p2$data)),
    "| PNG identico:", identical(px0, render_png_pixels(p2)), "\n")
for (lv in c("genus", "family", "species")) {
  a <- call_f(abundance_bar_plot, level = lv); b <- call_f(f_B2, level = lv)
  cat("  level =", lv, "B2 all.equal p$data:", isTRUE(all.equal(a$data, b$data)), "\n")
}

bm <- function(f, n = N) { t <- tr <- numeric(n)
  for (k in seq_len(n)) { gc(FALSE); t[k] <- system.time(p <- call_f(f))[["elapsed"]]
    tr[k] <- system.time(render_null(p))[["elapsed"]] }
  c(first = t[1], build_med = median(t[-1]), render_med = median(tr[-1])) }
res <- rbind("original (paquete)" = bm(abundance_bar_plot),
             "B1 filtro unico"    = bm(f_B1),
             "B2 B1 + rowsum"     = bm(f_B2))
print(round(res, 3))
write.csv(data.frame(variant = rownames(res), round(res, 3)), "results_prototype_abundance_bar_plot.csv", row.names = FALSE)
