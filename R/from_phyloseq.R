#' Convert a phyloseq object to MicroBioMeta tables
#'
#' Extracts the count table, taxonomy and sample metadata of a
#' \code{phyloseq} object into the plain data frames that every
#' MicroBioMeta function takes.
#'
#' @param physeq A \code{phyloseq} object with an OTU table and a
#'   taxonomy table and sample data.
#' @return A list with two data frames:
#'   \describe{
#'     \item{\code{table}}{Taxa in rows, samples in columns, and a last
#'       column \code{taxonomy} with the full taxonomic string
#'       (e.g. \code{"k__Bacteria; p__Firmicutes; c__Bacilli"}).}
#'     \item{\code{metadata}}{Sample metadata whose first column,
#'       \code{SAMPLEID}, holds the sample names of \code{table}.}
#'   }
#' @seealso \code{\link{from_tse}} for \code{(Tree)SummarizedExperiment}
#'   objects.
#' @export
#' @examples
#' if (requireNamespace("phyloseq", quietly = TRUE)) {
#'     counts <- matrix(c(10, 0, 5, 3, 8, 1),
#'         nrow = 3,
#'         dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2"))
#'     )
#'     tax <- matrix(
#'         c(
#'             "Bacteria", "Bacteria", "Bacteria",
#'             "Firmicutes", "Proteobacteria", NA
#'         ),
#'         nrow = 3, dimnames = list(rownames(counts), c("Kingdom", "Phylum"))
#'     )
#'     samples <- data.frame(
#'         Soil = c("Rhizosphere", "Bulk soil"),
#'         row.names = c("S1", "S2")
#'     )
#'     ps <- phyloseq::phyloseq(
#'         phyloseq::otu_table(counts, taxa_are_rows = TRUE),
#'         phyloseq::tax_table(tax),
#'         phyloseq::sample_data(samples)
#'     )
#'     mbm <- from_phyloseq(ps)
#'     mbm$table
#'     mbm$metadata
#' }
from_phyloseq <- function(physeq) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
        stop("Package 'phyloseq' is required but not installed.\n",
            "Install it with: BiocManager::install(\"phyloseq\")",
            call. = FALSE
        )
    }
    if (!methods::is(physeq, "phyloseq")) stop("`physeq` must be a phyloseq object.")

    counts <- methods::as(phyloseq::otu_table(physeq), "matrix")
    if (!phyloseq::taxa_are_rows(physeq)) counts <- t(counts)

    tax <- phyloseq::tax_table(physeq, errorIfNULL = FALSE)
    if (!is.null(tax)) tax <- methods::as(tax, "matrix")

    samples <- phyloseq::sample_data(physeq, errorIfNULL = FALSE)
    samples <- if (is.null(samples)) NULL else data.frame(samples, check.names = FALSE)

    .mbm_from_parts(counts, tax, samples)
}
