## Builds the small 16S example dataset used in the function examples and
## tests (subset_table_bacteria.txt, subset_metadata_bacteria.txt) from the
## full 16S dataset used in the vignettes (table_bacteria.txt,
## taxonomy_bacteria.txt, metadata_bacteria.txt).
##
## Subset: samples of month 2 from Rhizosphere and Roots (46 samples) and the
## ASVs with at least one read in them, with the taxonomy as the last column.
## In the metadata, Type_of_soil is named Location and the treatments are
## renamed (TC = Control, TD = Moderate_drought, TED = Severe_drought).
##
## Run from the package root: Rscript inst/script/make_subset_bacteria.R

extdata <- file.path("inst", "extdata")

table <- read.delim(file.path(extdata, "table_bacteria.txt"),
    row.names = 1, check.names = FALSE
)
taxonomy <- read.delim(file.path(extdata, "taxonomy_bacteria.txt"),
    check.names = FALSE
)
metadata <- read.delim(file.path(extdata, "metadata_bacteria.txt"),
    check.names = FALSE
)

metadata <- metadata[metadata$Month == 2 &
    metadata$Type_of_soil %in% c("Rhizosphere", "Roots") &
    metadata$SAMPLEID %in% colnames(table), ]

counts <- table[, metadata$SAMPLEID]
counts <- counts[rowSums(counts) > 0, ]
counts$taxonomy <- taxonomy$Taxon[match(rownames(counts), taxonomy$Feature.ID)]

colnames(metadata)[colnames(metadata) == "SAMPLEID"] <- "#SampleID"
colnames(metadata)[colnames(metadata) == "Type_of_soil"] <- "Location"
metadata$Treatment <- c(
    TC = "Control", TD = "Moderate_drought", TED = "Severe_drought"
)[metadata$Treatment]

write.table(cbind(OTUID = rownames(counts), counts),
    file.path(extdata, "subset_table_bacteria.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(metadata, file.path(extdata, "subset_metadata_bacteria.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
)
