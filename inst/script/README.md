# Example data in `inst/extdata`

All the example data come from published studies of the package authors.
Tables have taxa in rows and samples in columns; metadata tables have the
sample IDs in the first column.

## 16S rRNA (bacteria)

Source: Hereira-Pacheco, S.E., Navarro-Noya, Y.E. & Dendooven, L. (2021). The
root endophytic bacterial community of *Ricinus communis* L. resembles the seeds
community more than the rhizosphere bacteria independent of soil water content.
*Scientific Reports*, 11, 2173. https://doi.org/10.1038/s41598-021-81551-7

Sequences were processed in QIIME 2 (ASVs, SILVA taxonomy). The feature table
and the taxonomy were exported from the QIIME 2 artifacts (`.qza`) to
tab-delimited text.

- `table_bacteria.txt`: ASV counts (6045 ASVs x 255 samples). Used in the
  vignettes.
- `taxonomy_bacteria.txt`: taxonomy of each ASV (`Feature.ID`, `Taxon`,
  `Confidence`).
- `metadata_bacteria.txt`: sample metadata (294 samples), including soil
  compartment (`Type_of_soil`), treatment (`Treatment`), sampling month
  (`Month`) and soil physicochemical variables.
- `subset_table_bacteria.txt` and `subset_metadata_bacteria.txt`: a small
  subset used in the function examples and tests (46 samples of month 2 from
  Rhizosphere and Roots, and the 1704 ASVs with reads in them, with the
  taxonomy as the last column). In this subset `Type_of_soil` is named
  `Location` and the treatments are renamed (TC = Control,
  TD = Moderate_drought, TED = Severe_drought). It is built from the full
  files with `make_subset_bacteria.R`.

## 18S rRNA (fungi)

Source: Hereira-Pacheco, S.E., Estrada-Torres, A., Dendooven, L. &
Navarro-Noya, Y.E. (2023). Shifts in root-associated fungal communities under
drought conditions in *Ricinus communis*. *Fungal Ecology*, 63, 101225.

Sequences were processed in QIIME 2 (ASVs). The feature table and the
taxonomy were exported from the QIIME 2 artifacts (`.qza`) to tab-delimited
text.

- `table_fungi.txt`: ASV counts (211 ASVs x 65 samples).
- `taxonomy_fungi.txt`: taxonomy of each ASV (`Feature.ID`, `Taxon`,
  `Consensus`).
- `metadata_fungi.txt`: sample metadata (93 samples), with the same variables
  as the bacterial metadata.

## Shotgun metagenomes (fungi, Kraken2/Bracken)

Source: Hereira-Pacheco, S.E. et al. (2025). Metagenomic analysis of the
diversity and composition of the fungal communities of forest soils of the
transect Iztaccíhuatl-Popocatépetl and La Malinche. *PeerJ*, 13, e18323.
Raw data, metadata and analysis scripts:
https://github.com/Steph0522/Fungal_Metagenomes_PNIP_PNML

Reads were classified with Kraken2 (PlusPFP database) and abundances were
re-estimated at the species level with Bracken.

- `table_kraken.txt`: the `bracken_species` reports of the 36 samples merged
  into one table (87 fungal species; one column per sample, named
  `kraken_pluspfp_<id_metagenome>_report_bracken_species`), with the
  taxonomy as the last column.
- `metadata_metagenomic.txt`: sample metadata (36 samples): polygon
  (`Poligono`), site (`Sitio`), transect (`Transecto`), metagenome ID
  (`id_metagenome`) and soil physicochemical variables.
- `coord_metagenomic.csv`: latitude and longitude of the sampling points.
