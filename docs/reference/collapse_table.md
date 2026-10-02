# Collapse table Collapse an OTU/ASV table by taxonomic level

Collapses an abundance table by a specified taxonomic rank (e.g. genus,
family, phylum), summing counts across features that share the same
taxonomy. Features with lower taxonomic resolution than the selected
level are retained unchanged. Optionally converts counts to relative
abundance and exports the collapsed table to a tab-delimited file.

## Usage

``` r
collapse_table(
  table,
  metadata,
  level = "genus",
  rel_abun = FALSE,
  save_table = FALSE,
  table_filename = "collapsed_table.txt"
)
```

## Arguments

- table:

  A data frame containing an OTU/ASV abundance table. Must include a
  column with OTUID, one with taxonomy and sample columns with numeric
  counts.

- metadata:

  A data frame containing sample metadata. The first column is used to
  define the order of samples in the output table.

- level:

  Character. Taxonomic level to collapse to. One of `"kingdom"`,
  `"phylum"`, `"class"`, `"order"`, `"family"`, `"genus"` (default), or
  `"species"`. Case-insensitive.

- rel_abun:

  Logical. If TRUE, converts counts to relative abundance (%) per
  sample.

- save_table:

  Logical. If `TRUE`, saves the collapsed table to disk. Default
  `FALSE`.

- table_filename:

  Character. File path/name for the saved table (used when
  `save_table = TRUE`). Default `"collapsed_table.txt"`.

## Value

A list with two elements: `collapsed_table` (wide format, one row per
collapsed taxon) and `long_format` (the same data pivoted to one row per
taxon/sample combination).

## Examples

``` r
table_path <- system.file("extdata", "tabla_bacteria.txt", package = "MicroBioMeta")
table <- read.delim(table_path, row.names = 1, check.names = FALSE)

metadata_path <- system.file("extdata", "metadata_bacteria.txt", package = "MicroBioMeta")
metadata <- read.delim(metadata_path, check.names = FALSE)
colnames(metadata)[1] <- "SampleID"

collapse_table(
  table      = table,
  metadata   = metadata,
  level      = "genus",
  rel_abun   = FALSE,
  save_table = FALSE
)
#> $collapsed_table
#>                                                                                                                                                                         taxonomy
#> 374a5eb6496d14e8e9d6f3d8c7e34d9d                                                                                          d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae
#> 4c47d5cf81df1ac6e2a0560029a81578                                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae
#> 8fb6741a6685fad0374f85ad3e16156c                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae
#> 310f3de009c95de5c938b8e811dfc96f                                                    d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae
#> 7d3b39e8c0f0fdcb81b0c898a98c38f9                                                                                                                                     d__Bacteria
#> fc211549300b0954dbb4a4bf57e8a605                                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae
#> 53fad3d8c52f4b5a023ff92aee9c0f1a                                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae
#> 5e37342d650cde6cfc5bbe2726b12e3e                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Rhodocyclaceae
#> 56025f74cc8e54d4a054442d786ef096                                                                                                                d__Bacteria; p__Actinobacteriota
#> d9e51d98f1c86a1c426552d252c944e4                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae
#> bc905dcc609da83d3620f1986e63a9f7                                                                                          d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria
#> 90009df763e86e6eca7baa54f01551e6                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae
#> 6bbd12b392461134369d835715cd2765                                                   d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae
#> 83c88a72c2dab171c38aadd98aa2d389                                                                          d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales
#> 2b93708de9f1c24ff6140814a599a3a8                                                                                   d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales
#> 671a53a2acbf0c8c7e240e446daf0000                                                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales
#> bfe35644221421300c0c17e00e465c05                                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae
#> d12bdf20f24b92aa8ec94c6c3f60497f                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae
#> 3bdb72ce060bd89f107c5c2bc95399ed                                                                          d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales
#> e6a82be4b75aac8b928c1dfbe9ee9287                                                   d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Ilumatobacteraceae
#> e3a63ba0a0ec40fbf270a57c180e0d3c                                                     d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae
#> 2877cf981b128a18f788ba5437c1afa6                                                                                              d__Bacteria; p__Acidobacteriota; c__Acidobacteriae
#> c015bc50696fe3e3d94ace01d9dcd691                                                               d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Caldilineales; f__Caldilineaceae
#> eb2eeddf829080f8d9f25c6be5d507cd                                                           d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Opitutales; f__Opitutaceae
#> b1fcda6d8df4c25c8e6d7124c076b1fe                                                   d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae
#> de8bb5c39fb121fd048dfe0d639e6e8a                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae
#> a234ac03223c275785dedf6df3b2aeff                                                                                          d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria
#> 9ea1641eeed1f51fc228c653076d1d08                                                d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micromonosporales; f__Micromonosporaceae
#> 694b926113eb8291fea166446890c8fe                                                               d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Microscillaceae
#> b2af163540d17007833dc819704afb28                                                                              d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Frankiales
#> 6fc0ba1f8ff8f9259c9d492e59771755                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae
#> 7b8c1dba059eac7ec7f3cefa203f6b2d                                                      d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Frankiales; f__Geodermatophilaceae
#> 0b06028193ef1ce19f5339a8ba7d4448                                                                           d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micrococcales
#> 56b4ec05ef3bbded515ec3dd2d573197                                                               d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Gemmatales; f__Gemmataceae
#> fc9c84f8767611251f0eb93af4ca10db                                                d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micromonosporales; f__Micromonosporaceae
#> d3616bc0d8925a452274c484c73406c3                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> 7986b8b93f3afce920bd19a5d2229c09                                                     d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Reyranellales; f__Reyranellaceae
#> bc78ee7297b6912452be0865d113f7ad                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae
#> c472668b46b0d0f575108aa4170dd81f                                                                          d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales
#> 2658a78c5ff06760389ce2e3a548b4a2                                                          d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Devosiaceae
#> 600c9e81fd1d1ddd5a37a5aa0f9f796d                                                                    d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__Polyangiaceae
#> f8d0bb97f093d8e8630ea62ac8a1bf13                                                           d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Pirellulales; f__Pirellulaceae
#> 6c685b83d63d6d56a9dcaba34215ae6b                                                                      d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales
#> e0e7e6208280a5846a8b66f00667f8ca                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> 6d1aa1c1f93c423a856a5dd3b62e5c70                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> e553e766ea0e08f1f9fa2debe7946af6                                                               d__Bacteria; p__Firmicutes; c__Clostridia; o__Oscillospirales; f__Ruminococcaceae
#> 0f1ef43f6b23a567aba80599766d70bc                                            d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Propionibacteriales; f__Propionibacteriaceae
#> 3c4ccaa45e8d148c9a17b5d0f3720cd2                                                 d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Propionibacteriales; f__Nocardioidaceae
#> d9df6265d046a558d87d397b8a177d54                                                                  d__Bacteria; p__Myxococcota; c__Polyangia; o__Nannocystales; f__Nannocystaceae
#> 4bbad82a39655a26a2de61b88c5af7fa                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> f5a6865451a3b138f4915c5ab445a0a7                                                    d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae
#> 4a86fc550fa1ac37aa480d37edea875f                                                         d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Isosphaerales; f__Isosphaeraceae
#> 02fc21e72f66adea757b355edea4a369                                                                                             d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia
#> 25cf10f4f15dab00925cef52b9f08442                                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae
#> fa8c8b98d197c252368c160b43cefb98                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> c0b07d29aed1a875641aede53bc0a20c                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> 341433c1c8cc65ce47dcf32179d08fd8                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> 8bf175f329d16e4aa732cf2b32279df3                                                                     d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales
#> 062f7552b747dccab9586dff7275b5d9                                                   d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Ilumatobacteraceae
#> 059e124ad6ba8f98c3d21db28a075ea8                                                         d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae
#> 3d7a81faea3f783158352c9f739c8d17                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> 90be2fb2c82a98002d33ed7531e877e3                                                                                                                                     d__Bacteria
#> b0e263cdd1c189fc1757d8189db14be2                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> b962ce332fea338b05293d1a5a29cbca                                                   d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae
#> 7a310a8f8e4bcbfb717263bd5376994d                                                                  d__Bacteria; p__Myxococcota; c__Polyangia; o__Nannocystales; f__Nannocystaceae
#> 57136cbdae8fe066c2bb4a661003719b                                                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales
#> 69fe46cab1cc9407470a2397bfc84103                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Rhodocyclaceae
#> b88efba97658a6936bb053985d026ace                                                                                          d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria
#> 4d691fbbe2dbfb2815649fb72c757977                                           d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Steroidobacterales; f__Steroidobacteraceae
#> 49192cc8a222d97b965b839ff5d0f333                                                     d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micrococcales; f__Microbacteriaceae
#> 859873fcd088c294a92a13340744eb6f                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Methylophilaceae
#> 0a8c44ca4d744e313e6dc5b239fca562                                                                                                                      d__Bacteria; p__Firmicutes
#> 09f94345824c82e8fe0d824cc4478fb5                                                                      d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales
#> 9dfbdb741fb0dc25d2206b8680bbf579                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> dc2bea463c720874cccd189ccdcf0e3e                                                                      d__Bacteria; p__Desulfobacterota; c__Desulfuromonadia; o__Desulfuromonadia
#> 6273eb52991d528cf5281215f79b40a3                                                                        d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Planococcaceae
#> dc4bcbe74986e44ec3a040b866dfda9e                                                                   d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae
#> 96ec4781ee61f4083792f252f1b6cd3a                                                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales
#> 9657d52acca51e701e7a3b7cedcfcc12                                                   d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Azospirillales; f__Azospirillaceae
#> f65f183293d4710041057ab98f2c94ed                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae
#> 8f6b3783b31a2dd8640de3df92dbb9c4                                                                                                                                     d__Bacteria
#> 0d90eb4842c9ad9126ed33c538c768ac                                                        d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micrococcales; f__Micrococcaceae
#> 8d3ee46aaaf7728594014be61c649347                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Methylophilaceae
#> e011a95fda1f1f2bafc4e5b83b52d989                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> ef8ea65ba4c2e2da3ba9556a87d2385d                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> dbbb92eb3a8a3feaf96092d46f2b6a45                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> 06cb03cddd8a12dbfaf27e8984b395c4                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> cdf134bb501c54d997945dabd44e35ff                                                           d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Opitutales; f__Opitutaceae
#> f3c4561cc01a3b45ee3f134d709c6b94                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> b121442fd1eefe78ce4f4602aac1842e                                                        d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micrococcales; f__Micrococcaceae
#> ba8707aa88a2a62939d9338189b2a733                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae
#> ecf7289628d49f9244469337ee1622cb                                                    d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae
#> 9f06b400b93f859b11874b0bd7f09295                                                                      d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales
#> f775be5cf7c0dd2177006a44913503d5                                                   d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae
#> e87ceb10f2becb2a6c4cf0692b1923ee                                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae
#> 036574bb8e05f3b9a11ab7a4090d6950                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> 0267ba9db97e00ae01fc2cccbb6263e0                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> be4f701c1dbfbda16c5150af29430141                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> 655a66e6566798ade91f73b2942c0b63                                                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales
#> cd79f2696509d0ad906f5addba264079                                                                           d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae
#> 07ea3284bda2bf7daf1251ad2744e4c4                                                                              d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Frankiales
#> e3870da9de4a2b18a7ae9e9cae208e3e                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> e00bd59c1ea0ca9b1cc5f976801d4e1b                                                                   d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae
#> b39b381d0f533c9923c3b898ffab5614                                                                                                                                     d__Bacteria
#> e2959f0552fa6df5707d9ec8b57816d6                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> 53157acd1642f973842d01e1cd287421                                                d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Nitrosomonadaceae
#> 5991790d1c0aa8b73b2d7f906449a1f9                                                   d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Azospirillales; f__Azospirillaceae
#> adee3552c211d9e33aaa28e82532bdea                                                d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micromonosporales; f__Micromonosporaceae
#> 0089040d041888e3ae4aa4e0010a0e36                                                                      d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales
#> 8e876e13d051992c4c8616cf2db7d73a                                               d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Rhodanobacteraceae
#> 179abbc5c1f1549ffd0d7ae6963e3783                                                                       d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Corynebacteriales
#> 340d095f9a9fb445582b60e1282d985d                                                                                  d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Kallotenuales
#> 3d7e5432d8eaa3fa56ab7138605ea43d                                                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales
#> d4b233740db907d69035915e8657fa8d                                                                          d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales
#> b22119df3b1de76745abaac585b646f8                                                    d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae
#> 64761a24fc66e7c93145217dd153b1ee                                                                    d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__Polyangiaceae
#> eb6d7f9e0601207fdcfc4a88947d513c                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> 0d20271e3e2dbc68e1de868599a9a1ae                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae
#> 369ea22772cb774950c830d0af7e8235                                                                                                                d__Bacteria; p__Actinobacteriota
#> cae63066cfe088a49fcac50f119c9e2d                                                   d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Ilumatobacteraceae
#> d9c57ec20ed90e2cebf6cacc4f9a5511                                                    d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Hyphomicrobiaceae
#> f5828cc3078b4c219a44c67e96e5d455                                                                                                                      d__Bacteria; p__Firmicutes
#> 73847ac2778d4cba27fe9758ca129d6c                                                                                                                                     d__Bacteria
#> d78003f13712da48d0bb36fd45116b91                                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae
#> eb40853d986ad72c82658d939f291ff2                                                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales
#> 5d675a3518222fc99c8cba34feaeac81                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> 3ba390ae11c4e76375af082eefe7e41f                                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae
#> 3aa66daa4e560cf62c1833697bcd3837                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> 6f4388510be5ea70356f2257d8c57ddf                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae
#> b33390d034148edf07d8a742778ddcbd                                                               d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Gemmatales; f__Gemmataceae
#> 09225ec7044394bbd3b6811a6799e481                                                   d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae
#> fe6fbf9b86e0f487911a856009fa4a1e                                                         d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae
#> 0aee6329e0b5b81449980151bf40205a                                                      d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Microtrichaceae
#> 0ef66a78236f90916afc2305f3ff6c80                                                         d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae
#> fb1eb41359005c2e3e308fb668f07220                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> cb59c57800227006aead644df42b8986                                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae
#> 138da9dc244c3292536cacb8df79467b                                                                   d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae
#> 72f92361ec0696cb8df985b93d277a25                                                   d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae
#> 4d5ab121763879fd59e8f00d35473540                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> 20032b2b586fee22805eccebe3059170                                                d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micromonosporales; f__Micromonosporaceae
#> 710324ccc8cafa6fc9b38baf8299d8aa                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> 523e6b60340f4a215877b1c336822e50                                                               d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Microscillaceae
#> 4501f98c83484fd902be11383c8d032b                                                                                                                                     d__Bacteria
#> d4cbf23757aad57c91b058cbb05566e4                                                           d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Opitutales; f__Opitutaceae
#> de7964d5c22112693ed5796bb9f04b0b                                                      d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Frankiales; f__Geodermatophilaceae
#> f4a283be42989b920b83ea73ec72b25b                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> 2f3111649d0850ab86799f8bbadc28f8                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Rhodocyclaceae
#> b8acdd45cf33ab75fe4c0851a1787dc5                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae
#> 5bb7fec39744a6d4e4d61596dabc6886                                                           d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Opitutales; f__Opitutaceae
#> e8eddfb37c4215c86813163a27f202a5                                                             d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridia; f__Hungateiclostridiaceae
#> 826e188c911c0846fd927300127c80c0                                                                                                                                     d__Bacteria
#> ff1fe599c44074ee6bdf919a119df93c                                                          d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Devosiaceae
#> aff00055d3876f624d67fb80bdb59934                                                                                                       d__Bacteria; p__Firmicutes; c__Clostridia
#> 1beded194f3ae4b57715b45e50adef46                                                                                                                                     d__Bacteria
#> 908d91feba7c1f1d8cc8a7145173d942                                                     d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae
#> e7095743695f6d4337d01ceb9e29e224                                                                                              d__Bacteria; p__Acidobacteriota; c__Acidobacteriae
#> 4cd420962855afdf0bc69bb268b846c7                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> 430ba5bea0fc08f65b69497bfbeb5c4b                                                   d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae
#> e12cb2f0d2ce4401413176007dcd4158                                                   d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae
#> ec615213f6c7febb65c9f2e45e438111                                                                                                                                     d__Bacteria
#> 412dd52cd54d0244ba8fc42ad140acda                                                     d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae
#> 89192d9f779341d313f52aad0dc6a06f                                                     d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Beijerinckiaceae
#> ed543f92c61a06fc4e703655ce910a2f                                                                                      d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales
#> fea1a70b0ef0f72e3bf27d8b295331a8                                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae
#> 6afe805b36f5aa29e7f2096aac42ace1                                           d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiales_Incertae_Sedis
#> d9cc5a46831f24ec9192ae82dbc38b0a                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae
#> 0d3eb38eb502c4ad13bd5a7060295411                                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae
#> 69a0ff5538a5490e291c75ce70325066                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> f91d545c815effa3c5afa83a8ad3397d                                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae
#> e6b032dea8e13c34af958fc48af09b58                                                     d__Bacteria; p__Desulfobacterota; c__Desulfobulbia; o__Desulfobulbales; f__Desulfocapsaceae
#> 84a1b89045ab6285328f366224f3ec87                                                                   d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae
#> 4a6f0076e4ac71474b3c364294ad3b33                                                                           d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae
#> 7bdd39e4e3ccc7000862b487510e67d6                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> 998aedf12cdf054f5e504a2e21d20498                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> 567c2640108006c8eca081ad4180e346                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae
#> e50a6625b9176b6856764b213013150f                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> 8d5a5ab7c784b0325a4a86a16fbba80e                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae
#> f58ae28a309fa3d5898d2eb46d9eb9e0                                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae
#> 7ac65a8fd21111b68e7c2d6205112a83                                                             d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridia; f__Hungateiclostridiaceae
#> d7add64c8700f5a8cadf40c8c7f5ff67                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> 6bad9d981d3cf0f90262547180dc6371                                                                     d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales
#> 149c0e63c83113e1246c0e041de67d9c                                                         d__Bacteria; p__Verrucomicrobiota; c__Chlamydiae; o__Chlamydiales; f__Parachlamydiaceae
#> 6f196fae54ea90b5c11627c503c98940                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> b67c995d842e5b4b92fa0ce052603ef1                                           d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Solirubrobacterales; f__Solirubrobacteraceae
#> f44ffae67c7297e2c696b75fbd4e0f4c                                                           d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Opitutales; f__Opitutaceae
#> 33997912cd952feaec18fbd11b31acac                                                          d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Devosiaceae
#> 0a035423c4bf71a17566549f9b69448d                                                                                                                                     d__Bacteria
#> f65a6f1f1e5a4bc5c733ba364c0096d3                                                         d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae
#> f15389a9f245465c33db3de72aeb17ba                                                                  d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae
#> cd2b676d3c0785a60359c0d416701ee4                                                                                          d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria
#> 70ba2a378ebcd9a8f57b024a561cc6b1                                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae
#> 0273707384afdd41e9d022b539094211                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> f9ac98e876ecdb504202388d7ccdda5f                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> bdd9b87b782b8f935ecac62a9824e6a1                                                           d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Pirellulales; f__Pirellulaceae
#> 1f12b97495668f1fc378a349546d1806                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> 83f1b4d776f6e08a6ba6f86634a1e5a9                                           d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Solirubrobacterales; f__Solirubrobacteraceae
#> 49f614420449f6c7ecc8c9f4e137ad05                                                     d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Beijerinckiaceae
#> aba5aa00a1b262f27517eca64ee3333f                                                                   d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae
#> 8971980484eeabeef00145b93b699e30                                                                                      d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales
#> c1347ae54b0fc1dd385a1ccf3b8061ce                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> e9fe19783ebcc7e43a3a07ca8225db62                                                   d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae
#> 7c523a5f0a7c736efd82afef84137901                                                                          d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales
#> e13c07ec336c89db5183cc6e9e57ee6e                                                                                          d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria
#> a017ce431ccadea609b501d46e64e55b                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> fadf38d03d4b944f5782a759591ec3ae                                                                                                                                     d__Bacteria
#> 41bfda1c7f32205246e73a1e3e3ec4f4                                                                              d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Frankiales
#> 9e4f3e697b70758a7a996f3dd4533650                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> 39c4f6caa7510ed108f3ebcf72f79edf                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae
#> 4aa083fe7c1b0639dbc48ec649e68a31                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Rhodocyclaceae
#> fd2ac2e14b72b0e22a052e3ab4580471                                                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales
#> 4ee10294799bcfeab0c8ee567d5d9608                                                     d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae
#> 46869f12c31de7a05781e930c25a1ccd                                                                   d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae
#> 5e2548d6a970645d0ca943a72d207489                                                   d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Ilumatobacteraceae
#> 0538e1beed69f098c3f602130ee11a1c                                                                             d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales
#> 991c563278c4e8c91ffebce3560fe51b                                                   d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae
#> ce05f58164e61517198f7ab0e3304cb5                                                   d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae
#> fdbd0fed21bae91c05b2da268a025a89                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae
#> 7412d11595d77f26792071cd28a26760                                                           d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Opitutales; f__Opitutaceae
#> 4b087e13957eb38128f63e1c99df23c3                                                             d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Cyclobacteriaceae
#> 310a6f0e6841db2ce75005cac942a972                                                                                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria
#> 3df0fada803ac582802b47d4f4efe0c1                                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Acetobacterales; f__Acetobacteraceae
#> 3cd4e34c02ffdd1446df65bb54278c57                                                           d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Pirellulales; f__Pirellulaceae
#> b14a1dea93554257791f16c2e94906b0                                                                                                                                     d__Bacteria
#> 902a4435fc3af101e8a9b5f51550fa4d                                                 d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Propionibacteriales; f__Nocardioidaceae
#> 1995185b31349529f58964c2e5a3becd                                                     d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Beijerinckiaceae
#> 69b247c565afdeacee0240482b1ce536                                                             d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Anaerolineales; f__Anaerolineaceae
#> 95ff1b9528acca018e9c4697a8625504                                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae
#> 132d26d3ce6ac2ab25ccb4dd171879d7                                                  d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Streptomycetales; f__Streptomycetaceae
#> bf1a6e4f94ee8ba60cf79089139641b6                                                      d__Bacteria; p__Desulfobacterota; c__Desulfuromonadia; o__Geobacterales; f__Geobacteraceae
#> 704189aa09f64b992b3a02ee001dc706                                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae
#> dea5ea27839ee67a31d18ba0bf0f7f65                                                                      d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales
#> 3840a4fff0250232f8570cae7c3514ba                                                                          d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales
#> 31ecabc032997aaf41dd0d8def3ec93f                                                                d__Bacteria; p__Firmicutes; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae
#> 024a0e0348ed9abe9ae1c48b1b7ba09e                                                             d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Anaerolineales; f__Anaerolineaceae
#> fc1d9940419420113ce3fbfacc8d703a                                                        d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micrococcales; f__Micrococcaceae
#> 0770fbf28cd307a08bf557af7ffecfd3                                                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales
#> 883cbe8e4d1d98f468355e8ea974d222                                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae
#> cb6d6561383fb68f43f92a88a607749d                                                                                              d__Bacteria; p__Planctomycetota; c__Planctomycetes
#> e07e0f1a9fbc696a9fe4f9b7739ae45b                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae
#> 73d8d8f999acfa13cba2277f7526897e                                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae
#> 5218274c354d41965abbf4b0751d6f3a                                                                                                                                     d__Bacteria
#> 5c31570f706752fa4f7851aca5b5c291                            d__Bacteria; p__Abditibacteriota; c__Abditibacteria; o__Abditibacteriales; f__Abditibacteriaceae; g__Abditibacterium
#> 6a5e71c8b86a52ab1dd5e61164103940                                        d__Bacteria; p__Acidobacteriota; c__Acidobacteriae; o__Bryobacterales; f__Bryobacteraceae; g__Bryobacter
#> 25641b0b803c7d8493ebc02cb104bf3d                                                          d__Bacteria; p__Acidobacteriota; c__Acidobacteriae; o__PAUC26f; f__PAUC26f; g__PAUC26f
#> 8844d70dd1fdb96f2a7faa65bccf78c6                             d__Bacteria; p__Acidobacteriota; c__Acidobacteriae; o__Solibacterales; f__Solibacteraceae; g__Candidatus_Solibacter
#> 371b855afa8bed030e86ec78241b01a3                                                                d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__11-24; f__11-24; g__11-24
#> d0c180a47378d414a1d088ac14db236c                                   d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae; g__Aridibacter
#> a072a588b3940fab3b227e4901c40f7f                                 d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae; g__Blastocatella
#> 977a3523051e6b23da6221dbde795b34                                    d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Blastocatellales; f__Blastocatellaceae; g__uncultured
#> 1c3adc87d2e93b08e576953508eb680d                                                             d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__DS-100; f__DS-100; g__DS-100
#> 74eb0ddf379b00c5032b18fb9fbf1f32                                            d__Bacteria; p__Acidobacteriota; c__Blastocatellia; o__Pyrinomonadales; f__Pyrinomonadaceae; g__RB41
#> c335b34a15db7d532f511e6c952abe62                                                     d__Bacteria; p__Acidobacteriota; c__Holophagae; o__Subgroup_7; f__Subgroup_7; g__Subgroup_7
#> 0c78285397fb853a0e39069213510c54                                                 d__Bacteria; p__Acidobacteriota; c__Subgroup_11; o__Subgroup_11; f__Subgroup_11; g__Subgroup_11
#> 779dec6de9c6a8b7df28fc8ec3c75bee                                                 d__Bacteria; p__Acidobacteriota; c__Subgroup_18; o__Subgroup_18; f__Subgroup_18; g__Subgroup_18
#> 5526624b473bb1f0541dc6fb95c98678                                                 d__Bacteria; p__Acidobacteriota; c__Subgroup_25; o__Subgroup_25; f__Subgroup_25; g__Subgroup_25
#> 649e7fa7b9ce5c5472d34675ed643639                                                     d__Bacteria; p__Acidobacteriota; c__Subgroup_5; o__Subgroup_5; f__Subgroup_5; g__Subgroup_5
#> a9c5addeea542953d1bae4d30d599b02                    d__Bacteria; p__Acidobacteriota; c__Thermoanaerobaculia; o__Thermoanaerobaculales; f__Thermoanaerobaculaceae; g__Subgroup_10
#> 68291fb3b2558d438da4810961c8d53a            d__Bacteria; p__Acidobacteriota; c__Thermoanaerobaculia; o__Thermoanaerobaculales; f__Thermoanaerobaculaceae; g__Thermoanaerobaculum
#> 9c5d92bfefc7a190b93d5b9a6a77697b                                            d__Bacteria; p__Acidobacteriota; c__Vicinamibacteria; o__Subgroup_17; f__Subgroup_17; g__Subgroup_17
#> 3ddab42d701754625cab0c236c5b6e6d                                               d__Bacteria; p__Acidobacteriota; c__Vicinamibacteria; o__Subgroup_9; f__Subgroup_9; g__Subgroup_9
#> 287e840f9ecaef80d7a4762c8d6bf401                              d__Bacteria; p__Acidobacteriota; c__Vicinamibacteria; o__Vicinamibacterales; f__Vicinamibacteraceae; g__Luteitalea
#> 885f411fa54c5576148678f665bb053f                          d__Bacteria; p__Acidobacteriota; c__Vicinamibacteria; o__Vicinamibacterales; f__Vicinamibacteraceae; g__Vicinamibacter
#> b5709a57df9cbc8a741bb78a36f29c06                     d__Bacteria; p__Acidobacteriota; c__Vicinamibacteria; o__Vicinamibacterales; f__Vicinamibacteraceae; g__Vicinamibacteraceae
#> 26698cc60893017ef0b1ab173688d3de                              d__Bacteria; p__Acidobacteriota; c__Vicinamibacteria; o__Vicinamibacterales; f__Vicinamibacteraceae; g__uncultured
#> db64f0028b6518839ab7129c4b7ef72d                                       d__Bacteria; p__Acidobacteriota; c__Vicinamibacteria; o__Vicinamibacterales; f__uncultured; g__uncultured
#> f46fb43262da299ec7f0084298cafa50                                   d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Acidimicrobiales; f__Acidimicrobiaceae; g__uncultured
#> c45e7087de669be270881fc8ace0fad8                                           d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Actinomarinales; f__uncultured; g__uncultured
#> 854a20d5019388e81b1e78419eb017fd                                                   d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__IMCC26256; f__IMCC26256; g__IMCC26256
#> 9a2325ea015a8c5944047e3349cc1565                                                  d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Iamiaceae; g__Iamia
#> c5305ae3c6adc90814d70ab48f4f071a                         d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Ilumatobacteraceae; g__CL500-29_marine_group
#> e230123ec3ce8654942ee1a000a80010                                 d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Ilumatobacteraceae; g__Ilumatobacter
#> 14e3c60cb885509667670151bc4a1df0                                    d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Ilumatobacteraceae; g__uncultured
#> 7fd7b7bece92bbbfcceb2dcf1e835a8e                                        d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Microtrichaceae; g__IMCC26207
#> c098d1fb769dd975c2332dfb70da489e                             d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Microtrichaceae; g__Sva0996_marine_group
#> 95dcdf1eb06b252e67dcdbdaffc505e8                                       d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__Microtrichaceae; g__uncultured
#> c5399d0f9eda8f0815fd28a257952496                                            d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__Microtrichales; f__uncultured; g__uncultured
#> c7388cc314e992819b8f97ea35280982                                                d__Bacteria; p__Actinobacteriota; c__Acidimicrobiia; o__uncultured; f__uncultured; g__uncultured
#> ae8facc45d372f7108bd9499f74d1bfc                                                   d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__0319-7L14; f__0319-7L14; g__0319-7L14
#> 16caf6f5bc653fb6ba668eb821076a21                                    d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Actinomycetales; f__Actinomycetaceae; g__Actinomyces
#> e7e283b72ee15bcf0877c5eeb6137ebc                                d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Corynebacteriales; f__Mycobacteriaceae; g__Mycobacterium
#> 1b3a779b384c73b33207dccd14b14fe6                                      d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Corynebacteriales; f__Nocardiaceae; g__Rhodococcus
#> c8d1afdf7b395ffc22193c5dc2fd5dc3                                       d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Glycomycetales; f__Glycomycetaceae; g__Glycomyces
#> 0e3ffc2a68378086e2480831acd1a025                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micrococcales; f__Intrasporangiaceae; g__Ornithinimicrobium
#> 0ddf7339b608a340252b33f63fcc0019                                       d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micrococcales; f__Microbacteriaceae; g__Agromyces
#> 5fd4b2f7dee3fb0bde6f856ae266dbcb                            d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Micromonosporales; f__Micromonosporaceae; g__Virgisporangium
#> e2f4786bc66410c07bd8b96b5d0e2d92                               d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Propionibacteriales; f__Nocardioidaceae; g__Aeromicrobium
#> c30aae8f5216c9536492d0beae74950f                                d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Propionibacteriales; f__Nocardioidaceae; g__Nocardioides
#> 9087ecee0806a80861788a17dd7e6809                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Propionibacteriales; f__Propionibacteriaceae; g__uncultured
#> c29bc08418318360592b2b2717376a0a                            d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Pseudonocardiales; f__Pseudonocardiaceae; g__Actinophytocola
#> f791c7db1a78d35e84803b5e03c01965                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Pseudonocardiales; f__Pseudonocardiaceae; g__Pseudonocardia
#> 44b44fc4055a58c01f6d1e9de5e0d1d9                                 d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Streptomycetales; f__Streptomycetaceae; g__Streptomyces
#> 342df9401103d274851d3bd80704508d                             d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Streptosporangiales; f__Streptosporangiaceae; g__Nonomuraea
#> 0590e6672bb2c55c0ff54999a2b6dcfe                            d__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Streptosporangiales; f__Thermomonosporaceae; g__Actinomadura
#> d4e9147608bf883fdbbf9d23b17028d1                                                               d__Bacteria; p__Actinobacteriota; c__Coriobacteriia; o__OPB41; f__OPB41; g__OPB41
#> db24794eb8667aed224401386223f671                                                        d__Bacteria; p__Actinobacteriota; c__MB-A2-108; o__MB-A2-108; f__MB-A2-108; g__MB-A2-108
#> fe41e6fd9770030a21418b288930710f                                    d__Bacteria; p__Actinobacteriota; c__Rubrobacteria; o__Rubrobacterales; f__Rubrobacteriaceae; g__Rubrobacter
#> 699f18f3b155bff97de941f2cfbbbadc                                                 d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales; f__Gaiellaceae; g__Gaiella
#> 13870cf319997750e3b96bbb0b2f1aea                                               d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Gaiellales; f__uncultured; g__uncultured
#> 7ca4fb1b1740406241acb862445cf311                                                d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Solirubrobacterales; f__67-14; g__67-14
#> 9c0569fb95451d009fd3754303eaed12                          d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Solirubrobacterales; f__Solirubrobacteraceae; g__Conexibacter
#> 1ac26a3c12288d8fbadd6aa61366b1b7                      d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Solirubrobacterales; f__Solirubrobacteraceae; g__Parviterribacter
#> 89f3ba99552605501b8acb816d288afe                       d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Solirubrobacterales; f__Solirubrobacteraceae; g__Solirubrobacter
#> b4051b5227746f04c1feb9d8979299ff                            d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__Solirubrobacterales; f__Solirubrobacteraceae; g__uncultured
#> 087d358aaf77f81aaf50a623828865bb                                               d__Bacteria; p__Actinobacteriota; c__Thermoleophilia; o__uncultured; f__uncultured; g__uncultured
#> e451c4b3ff1c3dabeee296529b3a7a07                              d__Bacteria; p__Armatimonadota; c__Fimbriimonadia; o__Fimbriimonadales; f__Fimbriimonadaceae; g__Fimbriimonadaceae
#> 012871784f77c9dca7d446ddef2048d1                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Paludibacteraceae; g__Paludibacter
#> 95c6ce2788783bf2f83343deba23b8bd                                             d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Paludibacteraceae; g__uncultured
#> 213227f244d97553a9a14c4e3fd04c98                                              d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Prolixibacteraceae; g__WCHB1-32
#> c0e5e55c5cfbd576e7c085f1e4b69638                                      d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Aurantisolimonas
#> 51e98451d7e3f46127481ec5490f1c91                                          d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Chitinophaga
#> c75285a6c61d8fce20d823769c25c12a                                               d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Cnuella
#> 177966f7619eb24fdbce2266f5bed39a                                       d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Ferruginibacter
#> b89347f0534679f890628364d2c7b4fe                                 d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Flaviaesturariibacter
#> 8dff1cb12c0beb5bc7b2546480a6a6c7                                       d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Flavihumibacter
#> 84f2777cfcb0fbed6281917fcc765022                                       d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Flavisolibacter
#> a5e4d106a2ab947984fdef187996a0b9                                            d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Lacibacter
#> b6cbf6077b64acbb718971a1ffc97a27                                             d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Niastella
#> 27e29495fda424728fd6a80b83a16772                                      d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Parasegetibacter
#> 63c44ab3725dd0b803d43a29db9fd7fd                                            d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__Terrimonas
#> 9c931ab5bf18ef9b10996d59d292e6f5                                            d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Chitinophagaceae; g__uncultured
#> 12fd8efb82c3d4909543f6db9853134e                                              d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Chitinophagales; f__Saprospiraceae; g__uncultured
#> a8a60c20f7ab7d0f2dadd9f8858c19be                                              d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Cyclobacteriaceae; g__uncultured
#> d57e09eb1645f1abba6074d3bd7d4c52                                                   d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Cytophagaceae; g__Cytophaga
#> 7e0f2b982ca2d35b6ac2cf251d78d3c4                                              d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Cytophagaceae; g__Sporocytophaga
#> 3e2e8e938a68a0318bcf7145ad5ed4a9                                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Hymenobacteraceae; g__Adhaeribacter
#> f2117b49bf57d9451bbef6de5b762483                                             d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Hymenobacteraceae; g__Pontibacter
#> be24dfe003725a0b0f28f4136106faf8                                              d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Microscillaceae; g__Chryseolinea
#> e63301d7d0af55d5c43db051dc51522f                                              d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Microscillaceae; g__Ohtaekwangia
#> a6e2e5501507562f2174f65e606b2d43                                                d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Microscillaceae; g__uncultured
#> a06aa93cb50b9a39941a11c7242dde96                                                 d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Cytophagales; f__Spirosomaceae; g__Dyadobacter
#> 6a0a8abb36aa3bfdee14edbe97e9d781                                          d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Flavobacteriales; f__Crocinitomicaceae; g__Fluviicola
#> 4d3def5ffad0c392381591440bfc3668                                      d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Flavobacteriales; f__Flavobacteriaceae; g__Capnocytophaga
#> ce437b81e7cf5ce986f7240c26464b0a                                      d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Flavobacteriales; f__Flavobacteriaceae; g__Flavobacterium
#> e57788f9d2fae76ed9ca566ad4d475a8                                     d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Flavobacteriales; f__NS9_marine_group; g__NS9_marine_group
#> 8c525d0f63137bd068b75c57f5d04c09                                            d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Flavobacteriales; f__Weeksellaceae; g__Empedobacter
#> 6d670aa40dafa152c2d59cf2be9bef74                                                       d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__KD3-93; g__KD3-93
#> e08f8e9e9ffb10a702f2ef617f70eecf                                    d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Lentimicrobiaceae; g__Lentimicrobium
#> 7def1c807ca0767ccdfe5a6001dc51f1                           d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__NS11-12_marine_group; g__NS11-12_marine_group
#> 3b9a17275147adf1dc8c42eb00a06a36                                   d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Arcticibacter
#> 24be0b7e4c85ac53b09fddbca916ab96                                      d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Pedobacter
#> f521a47d3bb4acfec6b261f880b4eea8                                       d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__Sphingobacteriaceae; g__Solitalea
#> 91c24c9b6ffbedcc257d164176060fe0                                               d__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Sphingobacteriales; f__env.OPS_17; g__env.OPS_17
#> 9b45f210647c5fe71b005eb9e25299a5                                       d__Bacteria; p__Bacteroidota; c__Kapabacteria; o__Kapabacteriales; f__Kapabacteriales; g__Kapabacteriales
#> 974077d25fa6b32a7af22670b142a8f1                                                                 d__Bacteria; p__Bacteroidota; c__Kryptonia; o__Kryptoniales; f__BSV26; g__BSV26
#> 203aa2535e54a11e4c48af4ab946686b                                        d__Bacteria; p__Bacteroidota; c__Rhodothermia; o__Rhodothermales; f__Rhodothermaceae; g__Rhodothermaceae
#> 9064afe569448558eb08e5bf5e458bce                                                                        d__Bacteria; p__Bacteroidota; c__SJA-28; o__SJA-28; f__SJA-28; g__SJA-28
#> fad432eeeff7f1d923e6d60aafe18dcd                              d__Bacteria; p__Bdellovibrionota; c__Bdellovibrionia; o__Bacteriovoracales; f__Bacteriovoracaceae; g__Peredibacter
#> a29db0e5f4a04c5c6937745ccccb167b                              d__Bacteria; p__Bdellovibrionota; c__Bdellovibrionia; o__Bdellovibrionales; f__Bdellovibrionaceae; g__Bdellovibrio
#> f179c10f8e86873c90ceb14b974a114e                                                      d__Bacteria; p__Bdellovibrionota; c__Oligoflexia; o__0319-6G20; f__0319-6G20; g__0319-6G20
#> a23e8588a2240d7af6f65519158c902c                                            d__Bacteria; p__Bdellovibrionota; c__Oligoflexia; o__Oligoflexales; f__Oligoflexales; g__Oligoflexus
#> 52034e0da6d5e7962369437956984e18                                                d__Bacteria; p__Bdellovibrionota; c__Oligoflexia; o__Oligoflexales; f__uncultured; g__uncultured
#> 522d20c27cf4237b60f36ffdca056c19                                             d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Anaerolineales; f__Anaerolineaceae; g__Anaerolinea
#> 957615aefb421e6ad4728e610b2bdfa6                                                  d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Anaerolineales; f__Anaerolineaceae; g__UTCFX1
#> 5ccc74e70b8a82a6eebea3385738ece3                                              d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Anaerolineales; f__Anaerolineaceae; g__uncultured
#> e5924fd11a48419c6ec0e2a7cb4e43ee                                   d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Ardenticatenales; f__Ardenticatenaceae; g__Ardenticatenaceae
#> eef7be4b2531cb7e909b58e91af5640d                                          d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Ardenticatenales; f__Ardenticatenaceae; g__uncultured
#> 4416cbc2026260fa4e2ed0a8feff1443                                                 d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Ardenticatenales; f__uncultured; g__uncultured
#> 03f6091f5aacb69f06fc925acd1c9506                                               d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Caldilineales; f__Caldilineaceae; g__Litorilinea
#> f655e5ddf5074d5e71bfa583926dec30                                                d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__Caldilineales; f__Caldilineaceae; g__uncultured
#> 9f733d39e588c7310499dd444cdd7812                                                    d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__RBG-13-54-9; f__RBG-13-54-9; g__RBG-13-54-9
#> 2b413de5e05a1385bab229f50879628f                                                                        d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__SBR1031; f__A4b; g__A4b
#> 178b596f2fe5cb6cefc5513917988d84                                                                      d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__SBR1031; f__A4b; g__OLB13
#> aa94b014c726cf40983045bd606241c1                                                                d__Bacteria; p__Chloroflexi; c__Anaerolineae; o__SBR1031; f__SBR1031; g__SBR1031
#> 55310b3aac4868c4e33b0a9492d7609a                                  d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Chloroflexales; f__Chloroflexaceae; g__Candidatus_Chloroploca
#> c5d8d828df0e1baf5c2db05a6b8c4476                                              d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Chloroflexales; f__Chloroflexaceae; g__Chloronema
#> 6ceb26f60606233a7728d3457af84380                                                d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Chloroflexales; f__Chloroflexaceae; g__FFCH7168
#> 1ff83868d1c9ecf52d9b776e0b3fd896                                        d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Chloroflexales; f__Herpetosiphonaceae; g__Herpetosiphon
#> 3edfccfb9f81b340854e84ee5102a4c1                                               d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Chloroflexales; f__Roseiflexaceae; g__uncultured
#> 963dc552c551872800f2fcbc5663e5fd                                                          d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Kallotenuales; f__AKIW781; g__AKIW781
#> 3fda8c20e5da9fbc3fb7c5c864559019                                                    d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Thermomicrobiales; f__AKYG1722; g__AKYG1722
#> 6fd88fa7e883bdead3c5ba89da853e8d                                            d__Bacteria; p__Chloroflexi; c__Chloroflexia; o__Thermomicrobiales; f__JG30-KF-CM45; g__JG30-KF-CM45
#> 35826fe826b7a9f202e8bef7b39d9a47                                                                      d__Bacteria; p__Chloroflexi; c__Dehalococcoidia; o__S085; f__S085; g__S085
#> cc514f614ba146f489bf62b5a56074c5                                                     d__Bacteria; p__Chloroflexi; c__Gitt-GS-136; o__Gitt-GS-136; f__Gitt-GS-136; g__Gitt-GS-136
#> 26260de3858081ce55e6f220b3ba25df                                                 d__Bacteria; p__Chloroflexi; c__JG30-KF-CM66; o__JG30-KF-CM66; f__JG30-KF-CM66; g__JG30-KF-CM66
#> d6727a85782ffd2a069dccc23adaf521                                                                         d__Bacteria; p__Chloroflexi; c__KD4-96; o__KD4-96; f__KD4-96; g__KD4-96
#> eb1a74c2ca8831ef945775de0827ddda                                                                   d__Bacteria; p__Chloroflexi; c__Ktedonobacteria; o__C0119; f__C0119; g__C0119
#> 4ed3218615bbc1ed4da82d9806e4786b                                                                             d__Bacteria; p__Chloroflexi; c__OLB14; o__OLB14; f__OLB14; g__OLB14
#> 8a3367e36ed4e416d5a77a7946cf6fc5                                                                                 d__Bacteria; p__Chloroflexi; c__TK10; o__TK10; f__TK10; g__TK10
#> 56f45eba0496ed29fb7171041a2a950d                           d__Bacteria; p__Cyanobacteria; c__Sericytochromatia; o__Sericytochromatia; f__Sericytochromatia; g__Sericytochromatia
#> 9d2782f19e60ec0cf72e1dc46b1d1766                                                d__Bacteria; p__Deinococcota; c__Deinococci; o__Deinococcales; f__Deinococcaceae; g__Deinococcus
#> 1559f8b64dc2bdc8cdc26ff8a5ffa4cb                                   d__Bacteria; p__Desulfobacterota; c__Desulfobulbia; o__Desulfobulbales; f__Desulfobulbaceae; g__Desulfobulbus
#> 6d01f384a1cce140edea17dccd0f116d           d__Bacteria; p__Desulfobacterota; c__Desulfobulbia; o__Desulfobulbales; f__Desulfocapsaceae; g__[Desulfobacterium]_catecholicum_group
#> 5950949b71082741302b712e97309b17                          d__Bacteria; p__Desulfobacterota; c__Desulfovibrionia; o__Desulfovibrionales; f__Desulfovibrionaceae; g__Desulfovibrio
#> f6f238f95cc97995c9bf619b0c620b33                      d__Bacteria; p__Desulfobacterota; c__Desulfuromonadia; o__Desulfuromonadia; f__Desulfuromonadaceae; g__Desulfuromonadaceae
#> c7a112ffadf8009efd318625da44b18e                         d__Bacteria; p__Desulfobacterota; c__Desulfuromonadia; o__Desulfuromonadia; f__Geothermobacteraceae; g__Geothermobacter
#> e79f4f30b0dacfcc50f3c1b12804ca3f                                  d__Bacteria; p__Desulfobacterota; c__Desulfuromonadia; o__Geobacterales; f__Geobacteraceae; g__Citrifermentans
#> 297c825dc268f9a7b86b402edb6961b3                                                    d__Bacteria; p__Desulfobacterota; c__uncultured; o__uncultured; f__uncultured; g__uncultured
#> a1ed56ad4e86e0b4adc11baf2c4ab94b                                                           d__Bacteria; p__Elusimicrobiota; c__Elusimicrobia; o__FCPU453; f__FCPU453; g__FCPU453
#> 8afccb8f7f50c5603b471fb0c808c9e6                                                 d__Bacteria; p__Elusimicrobiota; c__Lineage_IIb; o__Lineage_IIb; f__Lineage_IIb; g__Lineage_IIb
#> ce41c973bd4965a400989ec128890d55                    d__Bacteria; p__Entotheonellaeota; c__Entotheonellia; o__Entotheonellales; f__Entotheonellaceae; g__Candidatus_Entotheonella
#> 5d99ba689ad34e4393a905e8a725b72b                           d__Bacteria; p__Entotheonellaeota; c__Entotheonellia; o__Entotheonellales; f__Entotheonellaceae; g__Entotheonellaceae
#> bfaa2d35e6dd13b9b32de811c6f92592                                 d__Bacteria; p__Fibrobacterota; c__Fibrobacteria; o__Fibrobacterales; f__Fibrobacteraceae; g__possible_genus_04
#> 0c544d5f3a49a010971a96d8c39e4905                                          d__Bacteria; p__Firmicutes; c__Bacilli; o__Alicyclobacillales; f__Alicyclobacillaceae; g__Tumebacillus
#> bab3f08e1d75f42987012cae451639ff                                                              d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Bacillus
#> fa97a81b8c33569e2aa1fc9a3137b435                                                         d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Fictibacillus
#> a4d08886c3c6c62efd5678339a07828b                                                           d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Geobacillus
#> 11ed0972601847ec168473a6bc4eeba7                                                          d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Bacillaceae; g__Ureibacillus
#> a6e4d88aa4b40402b1f94760d86b729c                                                  d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Planococcaceae; g__Paenisporosarcina
#> 86deaff47d33946ede1b7c25c6bbe9f9                                       d__Bacteria; p__Firmicutes; c__Bacilli; o__Bacillales; f__Salisediminibacteriaceae; g__Salipaludibacillus
#> a1d3f1f8c151eaeb5cf12179861e0c87                                               d__Bacteria; p__Firmicutes; c__Bacilli; o__Brevibacillales; f__Brevibacillaceae; g__Brevibacillus
#> 3980891c22ccbb387b65f9e437b2f94f                                          d__Bacteria; p__Firmicutes; c__Bacilli; o__Erysipelotrichales; f__Erysipelotrichaceae; g__Turicibacter
#> fc8ec321e0c2fc2a6de9abe98256b3f5                                           d__Bacteria; p__Firmicutes; c__Bacilli; o__Exiguobacterales; f__Exiguobacteraceae; g__Exiguobacterium
#> 84c11e6ffb29ebc2192f4b049b7dbddb                                                 d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Leuconostocaceae; g__Leuconostoc
#> 82104585f0b2617c13eb69f0bbdf7d5a                                               d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus
#> 61cd4918428d3b36910b9cb2fc45946a                                                d__Bacteria; p__Firmicutes; c__Bacilli; o__Paenibacillales; f__Paenibacillaceae; g__Ammoniphilus
#> 0ab9a976935482bce5f8af4cc6685147                                               d__Bacteria; p__Firmicutes; c__Bacilli; o__Paenibacillales; f__Paenibacillaceae; g__Paenibacillus
#> 87be3450216b45c53851fd173582a2d0                                              d__Bacteria; p__Firmicutes; c__Bacilli; o__Paenibacillales; f__Paenibacillaceae; g__Thermobacillus
#> 93a51bb503fadf90121f518db980e158                                       d__Bacteria; p__Firmicutes; c__Bacilli; o__Thermoactinomycetales; f__Thermoactinomycetaceae; g__Laceyella
#> b70cd843161996e075a853900e699d0f                                      d__Bacteria; p__Firmicutes; c__Bacilli; o__Thermoactinomycetales; f__Thermoactinomycetaceae; g__Planifilum
#> 154b8be9eea158f93d27d314d8e0e2c8                               d__Bacteria; p__Firmicutes; c__Bacilli; o__Thermoactinomycetales; f__Thermoactinomycetaceae; g__Thermoactinomyces
#> 70691c2737cf7dc7c194de819f0e6cda                      d__Bacteria; p__Firmicutes; c__Clostridia; o__Christensenellales; f__Christensenellaceae; g__Christensenellaceae_R-7_group
#> 2ae1d80ead92f5991f58da14d000acd6                                         d__Bacteria; p__Firmicutes; c__Clostridia; o__Christensenellales; f__Christensenellaceae; g__uncultured
#> 44fa299492152518fec34b9fa4554648                                       d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridia; f__Hungateiclostridiaceae; g__Pseudobacteroides
#> e15c8f9df87f363b0dbbf49de0bed249                                                 d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Caloramatoraceae; g__Fonticella
#> 73689cad6eb4e86d9d0bb48c0da24463                                  d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium_sensu_stricto_1
#> 234e15e82d64f2e1d61802fa8ea1998a                                 d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium_sensu_stricto_10
#> 88c76ddb1b24acffe65d8128e2f48244                                  d__Bacteria; p__Firmicutes; c__Clostridia; o__Clostridiales; f__Clostridiaceae; g__Clostridium_sensu_stricto_8
#> 5544139a16716cc6ef48b2c16c9613e9                                          d__Bacteria; p__Firmicutes; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae; g__Lachnoclostridium
#> 5b7fe68483e701059142cb531397b547                                                   d__Bacteria; p__Firmicutes; c__Clostridia; o__Peptococcales; f__Peptococcaceae; g__uncultured
#> 0a0875a0acf4a28d1bc1c4d1d11aaa64                          d__Bacteria; p__Firmicutes; c__Clostridia; o__Peptostreptococcales-Tissierellales; f__Anaerovoracaceae; g__Anaerovorax
#> 6b47ec9d1c70f07d16834e9a4e54aa70                 d__Bacteria; p__Firmicutes; c__Clostridia; o__Peptostreptococcales-Tissierellales; f__Peptostreptococcaceae; g__Sporacetigenium
#> f8752b4824d239cfb19f3216cc34941b                  d__Bacteria; p__Firmicutes; c__Clostridia; o__Peptostreptococcales-Tissierellales; f__Sedimentibacteraceae; g__Sedimentibacter
#> ec8a0768cd2fb426d7e6bfc2c4f50e5b                                                d__Bacteria; p__Firmicutes; c__Clostridia; o__Thermincolales; f__Thermincolaceae; g__Thermincola
#> 8a6b88f9515bf2dfac09e6a69045a4fe                                                d__Bacteria; p__Firmicutes; c__Limnochordia; o__Limnochordia; f__Limnochordia; g__Hydrogenispora
#> f5151c67af19c8471743ef3e2914f6b1                               d__Bacteria; p__Firmicutes; c__Negativicutes; o__Veillonellales-Selenomonadales; f__Sporomusaceae; g__Anaerospora
#> c3842619fcbde52781c512c2561fca0b                             d__Bacteria; p__Firmicutes; c__Negativicutes; o__Veillonellales-Selenomonadales; f__Veillonellaceae; g__Veillonella
#> b3b00faa9dee9fa6b3b33744f51da00e                                   d__Bacteria; p__Firmicutes; c__Negativicutes; o__Veillonellales-Selenomonadales; f__uncultured; g__uncultured
#> c4f04e28c6267ac6c5b05d403e8f91c3                                                       d__Bacteria; p__Firmicutes; c__Negativicutes; o__uncultured; f__uncultured; g__uncultured
#> 7fc06f107a6ddf6486ef2d43fcf23626                                  d__Bacteria; p__Firmicutes; c__Symbiobacteriia; o__Symbiobacteriales; f__Symbiobacteraceae; g__Symbiobacterium
#> 8d1dffcd377b338d2a0388637507b592                                       d__Bacteria; p__Firmicutes; c__Symbiobacteriia; o__Symbiobacteriales; f__Symbiobacteraceae; g__uncultured
#> 15410d762889d51818986bfcdb5d8965                                                             d__Bacteria; p__Gemmatimonadota; c__AKAU4049; o__AKAU4049; f__AKAU4049; g__AKAU4049
#> 48f68e44e258de63575b8bf6cf565d52                                d__Bacteria; p__Gemmatimonadota; c__Gemmatimonadetes; o__Gemmatimonadales; f__Gemmatimonadaceae; g__Gemmatimonas
#> 4712befe3685661e0fc0968859c369d5                             d__Bacteria; p__Gemmatimonadota; c__Gemmatimonadetes; o__Gemmatimonadales; f__Gemmatimonadaceae; g__Roseisolibacter
#> be9a3727a0422aea147100370d046fa8                                  d__Bacteria; p__Gemmatimonadota; c__Gemmatimonadetes; o__Gemmatimonadales; f__Gemmatimonadaceae; g__uncultured
#> 3c5942781761830f7b5406f551574424                              d__Bacteria; p__Gemmatimonadota; c__Longimicrobia; o__Longimicrobiales; f__Longimicrobiaceae; g__Longimicrobiaceae
#> 291d9f4d3d3c5cc3b0a1804770ed1e77                                  d__Bacteria; p__Gemmatimonadota; c__Longimicrobia; o__Longimicrobiales; f__Longimicrobiaceae; g__YC-ZSS-LKJ147
#> 784764c519a64dc9a147025a23222e93 d__Bacteria; p__Gemmatimonadota; c__S0134_terrestrial_group; o__S0134_terrestrial_group; f__S0134_terrestrial_group; g__S0134_terrestrial_group
#> b9e955ea5254dd68b54775aad2a67829                                       d__Bacteria; p__Halanaerobiaeota; c__Halanaerobiia; o__Halanaerobiales; f__Halanaerobiaceae; g__Halocella
#> 8197f1594e8a17f3fc989fea64f54b8f                                    d__Bacteria; p__Halanaerobiaeota; c__Halanaerobiia; o__Halanaerobiales; f__Halobacteroidaceae; g__uncultured
#> bf4a4f2383a16d2ecf27bf8e162b9ba4                            d__Bacteria; p__Latescibacterota; c__Latescibacterota; o__Latescibacterota; f__Latescibacterota; g__Latescibacterota
#> 4e7d22d793d30bdedbd32a5b73b564b8                              d__Bacteria; p__Methylomirabilota; c__Methylomirabilia; o__Rokubacteriales; f__Rokubacteriales; g__Rokubacteriales
#> e88959b739c8cceb94b94b7887dd1df3                                                    d__Bacteria; p__Methylomirabilota; c__Methylomirabilia; o__Rokubacteriales; f__WX65; g__WX65
#> 5403445694339e2ae5ecac06d9f3ff47                                      d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Anaeromyxobacteraceae; g__Anaeromyxobacter
#> ba3943fb075bbf63fb699638292826a1                                                    d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__Archangium
#> 67b2d836db093cead5a8286142669613                                                        d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__KD3-10
#> 6c3a0ea164d3eb098dd2baac5519b38f                                                    d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__Myxococcus
#> c643a776c36bcc00dfb5d790e1ea3dfe                                                       d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__P3OB-42
#> d53cc81a9d57a64d194c876b76bd66e0                                                    d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Myxococcaceae; g__uncultured
#> 7dc00d40971bb454b39fdd9b2dd9aa25                                            d__Bacteria; p__Myxococcota; c__Myxococcia; o__Myxococcales; f__Vulgatibacteraceae; g__Vulgatibacter
#> b4d0cec8543e940be1b5096370fdfb80                                                                   d__Bacteria; p__Myxococcota; c__Polyangia; o__Blfdi19; f__Blfdi19; g__Blfdi19
#> 6c54ceb94bc5e8fbd65224cffc60990a                                                     d__Bacteria; p__Myxococcota; c__Polyangia; o__Haliangiales; f__Haliangiaceae; g__Haliangium
#> b12ee9ca5a2320957213a114788cf2b7                                                                d__Bacteria; p__Myxococcota; c__Polyangia; o__MSB-4B10; f__MSB-4B10; g__MSB-4B10
#> 5e040c2397f959e7f0dfc6ec2855b812                                                  d__Bacteria; p__Myxococcota; c__Polyangia; o__Nannocystales; f__Nannocystaceae; g__Nannocystis
#> 3105a594bba2604ae5ed19f5e0495d2e                                                              d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__BIrii41; g__BIrii41
#> d75c040bdaa1cca96cd0655c9af7ab28                                             d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__Phaselicystidaceae; g__Phaselicystis
#> 9b909edaebabc90b9c1a703626e6bd0c                                               d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__Polyangiaceae; g__Pajaroellobacter
#> a6593da81035cf121dbe46346a4ea960                                                     d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__Polyangiaceae; g__Polyangium
#> 60e4c88400a750a8122b97476313a67e                                                     d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__Polyangiaceae; g__uncultured
#> 7854a78e3b881cee1c63ea0a50ba29ce                                                 d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__Sandaracinaceae; g__Sandaracinus
#> 1796734b56d45bdc997ccc09ab8731d1                                                   d__Bacteria; p__Myxococcota; c__Polyangia; o__Polyangiales; f__Sandaracinaceae; g__uncultured
#> 8243cd669ae7327f4326be78429ab05a                                                     d__Bacteria; p__Myxococcota; c__bacteriap25; o__bacteriap25; f__bacteriap25; g__bacteriap25
#> ebb93dc32c63402289114859d476cc5b                                                                                   d__Bacteria; p__NB1-j; c__NB1-j; o__NB1-j; f__NB1-j; g__NB1-j
#> 8a01af5a531e851e2c50fc01131cda33                                                d__Bacteria; p__Nitrospirota; c__Nitrospiria; o__Nitrospirales; f__Nitrospiraceae; g__Nitrospira
#> 03d2cc3f7b2bc352752d39e3738d0af9                                     d__Bacteria; p__Patescibacteria; c__Berkelbacteria; o__Berkelbacteria; f__Berkelbacteria; g__Berkelbacteria
#> ec378af6e32960b11b5f7391da4ae4fb                                 d__Bacteria; p__Patescibacteria; c__Gracilibacteria; o__Gracilibacteria; f__Gracilibacteria; g__Gracilibacteria
#> b407f2de868c064f65eb5ad8f82a8d0d        d__Bacteria; p__Patescibacteria; c__Parcubacteria; o__Candidatus_Moranbacteria; f__Candidatus_Moranbacteria; g__Candidatus_Moranbacteria
#> 73a95b6ba34fcd260989d848fe15d377                                         d__Bacteria; p__Patescibacteria; c__Parcubacteria; o__Parcubacteria; f__Parcubacteria; g__Parcubacteria
#> 3abb68202d15781590a2d1f9c5b37e15                                                     d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales; f__LWQ8; g__LWQ8
#> 77923ef8e9ed01c9b6a8ac85eb028163                                       d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales; f__Saccharimonadaceae; g__TM7a
#> 9981b8fd20f400e052527f84f5db1b66                           d__Bacteria; p__Patescibacteria; c__Saccharimonadia; o__Saccharimonadales; f__Saccharimonadales; g__Saccharimonadales
#> 597bdec49adf7675912f9fba4e262e70                                                                     d__Bacteria; p__Planctomycetota; c__BD7-11; o__BD7-11; f__BD7-11; g__BD7-11
#> 0fd3126118b13fb187bc55476a67c95a                                                                         d__Bacteria; p__Planctomycetota; c__OM190; o__OM190; f__OM190; g__OM190
#> 2a0c511ecc2513c42fe624e8e99d07da                                                              d__Bacteria; p__Planctomycetota; c__Phycisphaerae; o__CCM11a; f__CCM11a; g__CCM11a
#> d2912fea2dcf2e5b97e56c6a7d24295d                                          d__Bacteria; p__Planctomycetota; c__Phycisphaerae; o__Phycisphaerales; f__Phycisphaeraceae; g__AKYG587
#> 2976af6aaf4261b6c274fcc3229a4d65                                           d__Bacteria; p__Planctomycetota; c__Phycisphaerae; o__Phycisphaerales; f__Phycisphaeraceae; g__SM1A02
#> f097c45ce419da576aa111c98a4c7272                                  d__Bacteria; p__Planctomycetota; c__Phycisphaerae; o__Tepidisphaerales; f__Tepidisphaeraceae; g__Tepidisphaera
#> 3abc8260d34008f263e232ea25426c3f                              d__Bacteria; p__Planctomycetota; c__Phycisphaerae; o__Tepidisphaerales; f__Tepidisphaeraceae; g__Tepidisphaeraceae
#> d7e5b6432dc178aaeed9c81e0ec57afb                                d__Bacteria; p__Planctomycetota; c__Phycisphaerae; o__Tepidisphaerales; f__Tepidisphaerales; g__Tepidisphaerales
#> dfe71ce78ddfdbc7b453d17d628365f2                              d__Bacteria; p__Planctomycetota; c__Phycisphaerae; o__Tepidisphaerales; f__WD2101_soil_group; g__WD2101_soil_group
#> 289b18a86084273ce3b4038bf37e2840                                             d__Bacteria; p__Planctomycetota; c__Pla3_lineage; o__Pla3_lineage; f__Pla3_lineage; g__Pla3_lineage
#> 844a67032c52dc49019bf5dff248bf9c                                             d__Bacteria; p__Planctomycetota; c__Pla4_lineage; o__Pla4_lineage; f__Pla4_lineage; g__Pla4_lineage
#> 38cb3e3ce5ac5aed7dec50247ecfb267                                                   d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Gemmatales; f__Gemmataceae; g__Gemmata
#> 74e88e71e9929b1ddfcafccfab3b7bba                                               d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Gemmatales; f__Gemmataceae; g__Telmatocola
#> 66bc0737fc7fa1463ac7189ab7073998                                                d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Gemmatales; f__Gemmataceae; g__uncultured
#> b5a7aa588c9cae6ae8ca2f7d87ac53e8                                          d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Isosphaerales; f__Isosphaeraceae; g__uncultured
#> db0540920457b21b170947aac9137771                                       d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Pirellulales; f__Pirellulaceae; g__Blastopirellula
#> ebbaf260510ab16ac6b8c2b9866d05b9                                          d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Pirellulales; f__Pirellulaceae; g__Pir4_lineage
#> c1f461b67419d4a80a61db0867da9ed6                                             d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Pirellulales; f__Pirellulaceae; g__Pirellula
#> 0c6eba6452635eecfb865a1049299176                                            d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Pirellulales; f__Pirellulaceae; g__uncultured
#> 7ecabb10f29a0c44eb53e9f9852fe22a                              d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Planctomycetales; f__Rubinisphaeraceae; g__Planctomicrobium
#> a33cb3037e7825920073b4dddaaa41e6                                       d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Planctomycetales; f__Rubinisphaeraceae; g__SH-PL14
#> d625812c15c074c2e2b18099f97412c9                                    d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Planctomycetales; f__Schlesneriaceae; g__Planctopirus
#> 17ef77f12e2a55717e910a5f6a625448                                           d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__Planctomycetales; f__uncultured; g__uncultured
#> 7baf638fdd3f71166aafd22ecf410348                                                 d__Bacteria; p__Planctomycetota; c__Planctomycetes; o__uncultured; f__uncultured; g__uncultured
#> 9436cda211f60f722865cd23ef75d4ed                                  d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Acetobacterales; f__Acetobacteraceae; g__Roseomonas
#> c8ab0b20e740eb6b35e68d58f69736e4                                  d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Azospirillales; f__Azospirillaceae; g__Azospirillum
#> 5e43f2d8be71399fc8d6736754cb42ac                                   d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Azospirillales; f__Azospirillaceae; g__Skermanella
#> 19ad90ee1081a003406b02179ae63709                                         d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Azospirillales; f__uncultured; g__uncultured
#> 0a25cdbd0660e68970c563c43a18f54a                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae; g__Brevundimonas
#> 0aa6b210e902d8ab9d6f267407258ebb                            d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae; g__Phenylobacterium
#> a864c36d0d898cc5abd94d6d6f0bd5be                                  d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Caulobacteraceae; g__uncultured
#> dea0dcb305ff0f27e3c36c7133cf933d                                     d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Caulobacterales; f__Hyphomonadaceae; g__Hirschia
#> 463977baa69c7724fdc9407ad162714c                            d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Defluviicoccales; f__Defluviicoccaceae; g__Defluviicoccus
#> fcb72326c29342bcfb052953466eafc6                                                  d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Dongiales; f__Dongiaceae; g__Dongia
#> 5a554353f9ed728f0ab0631f29c8be2e                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Ferrovibrionales; f__Ferrovibrionaceae; g__Ferrovibrio
#> a9dd27b48e443415a743bf9902f515c4                                        d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Micavibrionales; f__uncultured; g__uncultured
#> 77a40bcf979c0099c52f19a3f1fe5217                                      d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Micropepsales; f__Micropepsaceae; g__uncultured
#> b4007ac525a9b409d70426e24421cff0                                      d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Reyranellales; f__Reyranellaceae; g__Reyranella
#> 406ba18a175c4c79efb25f48a6beafe4                                      d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Reyranellales; f__Reyranellaceae; g__uncultured
#> fe1c166755bc722cefd2e2a3b2b70672                                        d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Amb-16S-1323; g__Amb-16S-1323
#> 8afd274849ac07c0ca66f8f5f143df3c                                           d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Beijerinckiaceae; g__Bosea
#> 197301ed9d8728b5c061b99a17973161                                      d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Beijerinckiaceae; g__Microvirga
#> 14779f9dc3505598ad3d799e4ee676c4                                                      d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__D05-2; g__D05-2
#> 4085c246640eee043d523ce2d799023a                                              d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Devosiaceae; g__Devosia
#> bba09a403e67c1381de01c86511f1057                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Hyphomicrobiaceae; g__Hyphomicrobium
#> 392c59fee4e6b91d7654f2d3e3c65be8                                  d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Hyphomicrobiaceae; g__Pedomicrobium
#> d2b05ae73ad5ed296d7add9e2bb2cce4                                            d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__KF-JG30-B3; g__KF-JG30-B3
#> 8306a9473dabfd5918ea869dbfcd1291                                    d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Methyloligellaceae; g__uncultured
#> 06064c372bb624503c2905d13f7691f5  d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae; g__Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium
#> 3265692341e42dbd5cf4ed1938b906e3                                       d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiaceae; g__Mesorhizobium
#> da96a0797b9881a123795838b9f01140                              d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiales_Incertae_Sedis; g__Nordella
#> 0bdddf10e577657becec88e42001fd16                         d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiales_Incertae_Sedis; g__Phreatobacter
#> 168999a950ebe17e6d2ea1f90e557d91                            d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Rhizobiales_Incertae_Sedis; g__uncultured
#> a266ffe9d38f626d106fa615da4618ef                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae; g__Bradyrhizobium
#> 7706405d5bce7715f01d8b54188def42                                   d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae; g__Pseudolabrys
#> 4b0cfc8899b466bd9992f7b22a5a06be                              d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae; g__Pseudorhodoplanes
#> b81b2b4a7ca9eb2c23b170648dfaf22d                                    d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae; g__Rhodoplanes
#> 75cf40872201f150cb5b191676f64b78                                     d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__Xanthobacteraceae; g__uncultured
#> c8e5794c58de9326084a62d3282d257e                                            d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales; f__uncultured; g__uncultured
#> b0706a571f169c0bac7b4864052fcaa7                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Actibacterium
#> 0388f969e3527b2938acc676d6757bdb                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Amaricoccus
#> 2dfe751390aedfbd9f62dbb77d1f10b4                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Limibaculum
#> 2548dc9697e1f9afe04a9012f40e2f0a                           d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Pseudorhodobacter
#> 2f019615d5141a4ab00b1aab98882a33                                 d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Rhodobacter
#> 629e4f0cf1826904ff74907b2507f58b          d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Roseobacter_clade_CHAB-I-5_lineage
#> 774861f348cec8451102730c21a790b0                            d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodobacterales; f__Rhodobacteraceae; g__Rubellimicrobium
#> 11c0d6de28ab026a9ae17abe4efdd6bf                        d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Magnetospirillaceae; g__Magnetospirillum
#> c11381dda0841c0d08cbfcef3a760aa3                        d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhodospirillales; f__Magnetospirillaceae; g__Telmatospirillum
#> 3c9619d4f9db4ac6ad7151a0c9d86ad3                           d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Novosphingobium
#> 7e2c72204a8c45b1b0ac39e0e6047098                               d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Sphingobium
#> 8254e0f9f07709b62dd8483f3b25fad0                              d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Sphingomonadales; f__Sphingomonadaceae; g__Sphingomonas
#> e78135f2886d427e678d12c446bf070e                                      d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Thalassobaculales; f__uncultured; g__uncultured
#> abec938a8d53e7471f3f476d50694be9                        d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Tistrellales; f__Geminicoccaceae; g__Candidatus_Alysiosphaera
#> 4963567d2979bb761291aa98e0b8595b                                      d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Tistrellales; f__Geminicoccaceae; g__uncultured
#> d4de9d1453102564f35262faa3952266                                             d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__uncultured; f__uncultured; g__uncultured
#> efa4a78c4b993e9d115841c8a02c795e                                       d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Aeromonadales; f__Aeromonadaceae; g__Aeromonas
#> fcacd7424bf3543322ea3e9842d14e86                               d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Alteromonadales; f__Alteromonadaceae; g__Alishewanella
#> f941e15ea7f7560191d57abc4d2bedd8                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__Achromobacter
#> ac294f5242de3c9b26e75502c48f6383                                     d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Alcaligenaceae; g__Advenella
#> 7b6a8e3f8bfa1bcd1b5fc55db048d005                                                d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__B1-7BS; g__B1-7BS
#> c9e776ab599fe75db286675b698c9067                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Cupriavidus
#> ebe5cb945f6d6d849c07424591f80753                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Burkholderiaceae; g__Lautropia
#> f116387f77802409103a559e4205df06                                    d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Acidovorax
#> 93f9939b261cddc9a57071b42771c1d4                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Aquabacterium
#> 971a5f48b3418449ca7cf3e190a5e831                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Azohydromonas
#> 730de992ca06a95e7193a93146a97cf2                                     d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Comamonas
#> 8fa8b4ac69784c26901e036d13ab2332                                       d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Delftia
#> c363c9efaa85023b488b732b97e47270                                d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Hydrogenophaga
#> 2e2181cc893ca1e1e677e888bf8af441                                     d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Mitsuaria
#> 86fa662da6869634e7f90f85e9ee4933                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Piscinibacter
#> 9a02177f2c235a6eb62984d5aab151ff                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Pseudorhodoferax
#> d6854be498e22b65da3c00572ee178f0                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Ramlibacter
#> 015f364f150b230450da0a123e7bb8ab                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Rhizobacter
#> af2c1b22758cd6653e8e841497d096e7                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Sphaerotilus
#> e7f2fab8c0e7c365768a43ef71406605                                    d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Comamonadaceae; g__Variovorax
#> fb57cf192bce1876471f2f6181bf4b55                             d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Methylophilaceae; g__Methylobacillus
#> 6c740ce69c9f45bea684dc9e36054006                               d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Methylophilaceae; g__Methylotenera
#> 3b27c7bf279535d8e39664d8dc1d0daf                                     d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Methylophilaceae; g__UBA6140
#> ceea3c29ad012c2918fb6aa21a44e8cb                                      d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Neisseriaceae; g__Neisseria
#> fa99ecdb27c1cc44fd866ee1da250c4c                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Nitrosomonadaceae; g__Ellin6067
#> 387cc83573db0e444e81c89a8a4f8f9c                                       d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Nitrosomonadaceae; g__MND1
#> c7a8d1af5f7f472671f700775de3e0a2                                     d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Nitrosomonadaceae; g__mle1-7
#> 98a6a84550e07ea252cfd8a0f51042d4                                    d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae; g__Massilia
#> 08b809a68126c38db786cafe9710b013                          d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Oxalobacteraceae; g__Noviherbaspirillum
#> 8db4cc0a5d2f03fce8efd2b8374ff359                                      d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Rhodocyclaceae; g__Azoarcus
#> 3e7c7d38423c93b43a1d9865cf0bc4a6                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__SC-I-84; g__SC-I-84
#> 6f02050760c030ed8cd1ed836697fa6b                                    d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__Sutterellaceae; g__uncultured
#> 8fc1865cfcae32f72c1b1570299e2ef3                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Burkholderiales; f__TRA3-20; g__TRA3-20
#> 38cc138fac95b006929361093d840a00                                                            d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__CCD24; f__CCD24; g__CCD24
#> f5f8eb63e22362509e2c7c78271305f2                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Cellvibrionales; f__Cellvibrionaceae; g__Cellvibrio
#> 0bb238aaee3ccb7dc8dfa86dbc1ba9f7                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Cellvibrionales; f__Halieaceae; g__OM60(NOR5)_clade
#> ca4546caf07edf5200cf79cea04aebcd                             d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Diplorickettsiales; f__Diplorickettsiaceae; g__Aquicella
#> 893d12123eb7f623489a4f025a15ae7f                                         d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Erwiniaceae; g__Pantoea
#> 192c5a817f48368339e519e0a0a907da                                       d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Yersiniaceae; g__Serratia
#> 3467185ca92b7bec2523214605adabf0                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Gammaproteobacteria; f__Gammaproteobacteria; g__Gammaproteobacteria
#> a2aae80d8fb922fe1074d243b59ce237                d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Gammaproteobacteria_Incertae_Sedis; f__Unknown_Family; g__Acidibacter
#> f11233f16d4177d38ab1c0bfc649e056       d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Gammaproteobacteria_Incertae_Sedis; f__Unknown_Family; g__Candidatus_Berkiella
#> b219ca7793ab9f5e441f2a68bc76152c                                d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Oceanospirillales; f__Alcanivoracaceae; g__Ketobacter
#> e6e1848b51e705e55b3805703f088311                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Oceanospirillales; f__Halomonadaceae; g__Halomonas
#> bd074b9bf365d11468bd76268d63b15b                                                         d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__PLTA13; f__PLTA13; g__PLTA13
#> f381b492d6651c52fd09ce9351e7eb2d                                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pasteurellales; f__Pasteurellaceae; g__Haemophilus
#> d9b31e1c30facc24fa778b1d65f5c457                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Acinetobacter
#> 14a67c8c8cc31d203eac8109443fab5f                                      d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__Cavicella
#> 026ad0ea53fb202ebb770acf65705d6e                                     d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Moraxellaceae; g__uncultured
#> 19ba9fea04969f975d6dc8d3da20d5e6                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Azotobacter
#> 0fc4426ed7dccc48ce27cfff8083433d                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas
#> 45ad90e179842705fa984962f131fd61                                                            d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__R7C24; f__R7C24; g__R7C24
#> d284afab226d22fd0596221d17db23d3                        d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Steroidobacterales; f__Steroidobacteraceae; g__Steroidobacter
#> 38194655ba17dbfd7fcc1ed42a57c169                            d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Steroidobacterales; f__Steroidobacteraceae; g__uncultured
#> e30be5268bbebe597424bdf171cc1607                   d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Steroidobacterales; f__Woeseiaceae; g__JTB255_marine_benthic_group
#> cb7ccdb4956134b452a3cb78d1751c57                                       d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Steroidobacterales; f__Woeseiaceae; g__Woeseia
#> b585436677ac88dd4ac6239655f2a268                                              d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Vibrionales; f__Vibrionaceae; g__Vibrio
#> f5d770eff5c9b9b0c859db57d175eb34                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Rhodanobacteraceae; g__Ahniella
#> 1b45ffbc88c788dcee2239e68901cf05                               d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Rhodanobacteraceae; g__Luteibacter
#> 3fd951cf7cc3a0740a1f20f8a9c6d2c2                                d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Rhodanobacteraceae; g__uncultured
#> 2027d3e90ccc4003f1954322e99db6a9                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__Arenimonas
#> 6157cb92771830955413f808e83fc420                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__Luteimonas
#> 6349506cac32056a214b04239bdc3a33                                  d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__Lysobacter
#> 2254a8d37cf537b1e35097d70b6ac91b                           d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__Pseudoxanthomonas
#> 38f3bc8ef8836574de0fcbeafd348b65                            d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__Stenotrophomonas
#> fec23c70ae0c70ecceea1487cd978ee5                                 d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Xanthomonadales; f__Xanthomonadaceae; g__Thermomonas
#> 503d123a6c78c6fb280a6472bd61b51e                                                     d__Bacteria; p__Sumerlaeota; c__Sumerlaeia; o__Sumerlaeales; f__Sumerlaeaceae; g__Sumerlaea
#> 9f33de8570c639e727465f1ee4d750f6                                          d__Bacteria; p__Verrucomicrobiota; c__Chlamydiae; o__Chlamydiales; f__Parachlamydiaceae; g__uncultured
#> d68e1899fa65515135d9bd8d356b7e3d                d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Chthoniobacterales; f__Chthoniobacteraceae; g__Candidatus_Udaeobacter
#> e15b88ccac7a4bef64215e915f6e0d89                        d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Chthoniobacterales; f__Chthoniobacteraceae; g__Chthoniobacter
#> 767042b4b19d8bfaaa426dad7dbce378                          d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Chthoniobacterales; f__Terrimicrobiaceae; g__Terrimicrobium
#> 4a3c4abee502d6185abd750522566bf3       d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Chthoniobacterales; f__Xiphinematobacteraceae; g__Candidatus_Xiphinematobacter
#> c254aa5c00b7f27d756f7bbe96b3d4ad                                             d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Opitutales; f__Opitutaceae; g__IMCC26134
#> fe8d8af791eb6e5d4335bb7604c74e9d                                              d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Opitutales; f__Opitutaceae; g__Opitutus
#> d81800d7e31b69877ce1eab9858b8f3a                                d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae; g__ADurb.Bin063-1
#> c55c7cb80c5700ff9aded76192543b5a                                        d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae; g__DEV008
#> 91dee7b3b90bbf2e53b2472b9bbde0e9                                      d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae; g__Ellin517
#> 6b28df61cb5dee74b76b1a13aaa0027d                                    d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae; g__Oikopleura
#> 9a830fa26205ad91a02b5bcaebb0221a                               d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae; g__Pedosphaeraceae
#> a54b8142d2a09d058b81c377e8e286b8                                    d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Pedosphaerales; f__Pedosphaeraceae; g__uncultured
#> 69ef3a03a4d733dbffd99fdaa1d2e62a                                             d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Verrucomicrobiales; f__DEV007; g__DEV007
#> 083cbb131dd305a88b7d338bc371c35f                              d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Verrucomicrobiales; f__Rubritaleaceae; g__Luteolibacter
#> ef1a92caffb4d7cae286fdc77b0a93a6                           d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Verrucomicrobiales; f__Verrucomicrobiaceae; g__Brevifollis
#> 2cddd5dca166486f50a0828b28001e1e                       d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Verrucomicrobiales; f__Verrucomicrobiaceae; g__Prosthecobacter
#> bb899d1dcfd89cb200a76cb8951904f8                        d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Verrucomicrobiales; f__Verrucomicrobiaceae; g__Roseimicrobium
#> 7817f6851e9df03a3d992cf915310dba                            d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae; o__Verrucomicrobiales; f__Verrucomicrobiaceae; g__uncultured
#> 54eaf6571de743a86ccb531bee8aad7e                                                                                             d__Bacteria; p__WS2; c__WS2; o__WS2; f__WS2; g__WS2
#>                                  3.16RI 2.14RO 3.10RO 3.10RI 4.16RI 1.8RO 4.1RI
#> 374a5eb6496d14e8e9d6f3d8c7e34d9d      0      0      0      0     25     0    35
#> 4c47d5cf81df1ac6e2a0560029a81578      0      0      0     38      0     0     0
#> 8fb6741a6685fad0374f85ad3e16156c      0      0      0      0      0     0     0
#> 310f3de009c95de5c938b8e811dfc96f      0      0      0     31     33    27     0
#> 7d3b39e8c0f0fdcb81b0c898a98c38f9      0      0      0      0      0     0     0
#> fc211549300b0954dbb4a4bf57e8a605      0      0      0      0      0     0     0
#> 53fad3d8c52f4b5a023ff92aee9c0f1a      0      0      0      0     22     0     0
#> 5e37342d650cde6cfc5bbe2726b12e3e      0      0      0      0      0     0     0
#> 56025f74cc8e54d4a054442d786ef096      0      0      0      0      0     0     0
#> d9e51d98f1c86a1c426552d252c944e4      0      0      0      0      0     0     0
#> bc905dcc609da83d3620f1986e63a9f7      0      0      0      0      0     0     0
#> 90009df763e86e6eca7baa54f01551e6      0      0      0      0      0     0     0
#> 6bbd12b392461134369d835715cd2765      0      0      0      0     19     0     0
#> 83c88a72c2dab171c38aadd98aa2d389     11      0      0      0      0     0     0
#> 2b93708de9f1c24ff6140814a599a3a8      0      0      0      0      0     0     0
#> 671a53a2acbf0c8c7e240e446daf0000      0      0      0      0     28     0     0
#> bfe35644221421300c0c17e00e465c05     18      0     26      0     82     8    74
#> d12bdf20f24b92aa8ec94c6c3f60497f      0      0      0      0      0     0     0
#> 3bdb72ce060bd89f107c5c2bc95399ed      0      0      0      0      0     0     0
#> e6a82be4b75aac8b928c1dfbe9ee9287      0      0      0      0      0     0    30
#> e3a63ba0a0ec40fbf270a57c180e0d3c      0      0      0      0      0     0     0
#> 2877cf981b128a18f788ba5437c1afa6      0      0      0      0      0     0     0
#> c015bc50696fe3e3d94ace01d9dcd691      0      0      0      0      0     0     0
#> eb2eeddf829080f8d9f25c6be5d507cd      0      0      0      0      0     0    16
#> b1fcda6d8df4c25c8e6d7124c076b1fe      0      0      0      0      0     0     0
#> de8bb5c39fb121fd048dfe0d639e6e8a      0      0      0      0      0     0     0
#> a234ac03223c275785dedf6df3b2aeff      0      0     11      0      0     0    12
#> 9ea1641eeed1f51fc228c653076d1d08      0      0      0      0      0     0     0
#> 694b926113eb8291fea166446890c8fe      0      0      0      0      0     0     0
#> b2af163540d17007833dc819704afb28      0      0      0      0      0     0     0
#> 6fc0ba1f8ff8f9259c9d492e59771755      0      0      0      0     20     0     0
#> 7b8c1dba059eac7ec7f3cefa203f6b2d      0      0      0      0      0     0     0
#> 0b06028193ef1ce19f5339a8ba7d4448      0      0      0      0      0     0    35
#> 56b4ec05ef3bbded515ec3dd2d573197      0      0      0      0      0     0     0
#> fc9c84f8767611251f0eb93af4ca10db      0      0      0     53      0     0     0
#> d3616bc0d8925a452274c484c73406c3      0      0      0      0    101     0     0
#> 7986b8b93f3afce920bd19a5d2229c09      0      0      0      0      0     0    15
#> bc78ee7297b6912452be0865d113f7ad      0      0     21      0      0     0     0
#> c472668b46b0d0f575108aa4170dd81f      0      0      0      0      0     0     0
#> 2658a78c5ff06760389ce2e3a548b4a2      0      0      0     53     77     0     0
#> 600c9e81fd1d1ddd5a37a5aa0f9f796d      0      0      0      0      0     0     0
#> f8d0bb97f093d8e8630ea62ac8a1bf13      0      0      0      0      0     0    45
#> 6c685b83d63d6d56a9dcaba34215ae6b      0      0      0      0      0     0    19
#> e0e7e6208280a5846a8b66f00667f8ca      0      0     49      0    157     0    99
#> 6d1aa1c1f93c423a856a5dd3b62e5c70      0      0      0      0      0     0     0
#> e553e766ea0e08f1f9fa2debe7946af6      0      0      0      0      0     0     0
#> 0f1ef43f6b23a567aba80599766d70bc      0      0      0      0      0     0     0
#> 3c4ccaa45e8d148c9a17b5d0f3720cd2     42      0    131    173     99     0   141
#> d9df6265d046a558d87d397b8a177d54      0      0      0      0      0     0     0
#> 4bbad82a39655a26a2de61b88c5af7fa     41      0      0      0      0     0     0
#> f5a6865451a3b138f4915c5ab445a0a7      0      0      0      0      0     0     0
#> 4a86fc550fa1ac37aa480d37edea875f      0      0      0      0      0     0     0
#> 02fc21e72f66adea757b355edea4a369      0      0      0     22      0     0     0
#> 25cf10f4f15dab00925cef52b9f08442      0      0      0      0      0     0     0
#> fa8c8b98d197c252368c160b43cefb98      0      0     33      0     45     0     0
#> c0b07d29aed1a875641aede53bc0a20c      0      0      0      0     33     0     0
#> 341433c1c8cc65ce47dcf32179d08fd8      0      0      0      0      6     0     0
#> 8bf175f329d16e4aa732cf2b32279df3      0      0      0      0      0     0     0
#> 062f7552b747dccab9586dff7275b5d9      0      0      0      0      0     0     0
#> 059e124ad6ba8f98c3d21db28a075ea8      0     53    104    147     99    51    58
#> 3d7a81faea3f783158352c9f739c8d17      0     67      0     93      0     0     0
#> 90be2fb2c82a98002d33ed7531e877e3      0      0      0      0      0     0     0
#> b0e263cdd1c189fc1757d8189db14be2      0     41     36     57    431    67    79
#> b962ce332fea338b05293d1a5a29cbca      0      0      0      0      0     0     0
#> 7a310a8f8e4bcbfb717263bd5376994d      0      0      0      0      6     0     0
#> 57136cbdae8fe066c2bb4a661003719b      0      0      0      0      0     0    18
#> 69fe46cab1cc9407470a2397bfc84103      0      0      0      0     12     0     0
#> b88efba97658a6936bb053985d026ace      0      0      0      0     18     0     0
#> 4d691fbbe2dbfb2815649fb72c757977      0      0      0      0      0     0     0
#> 49192cc8a222d97b965b839ff5d0f333      0      0      0     17      0     0     0
#> 859873fcd088c294a92a13340744eb6f      0      0      0      0      0     0     0
#> 0a8c44ca4d744e313e6dc5b239fca562      0      0      0      0      0     0    15
#> 09f94345824c82e8fe0d824cc4478fb5      0      0      0      0      0     0     0
#> 9dfbdb741fb0dc25d2206b8680bbf579      0      0      0      0      0     0     0
#> dc2bea463c720874cccd189ccdcf0e3e      0      2      0      0      0     0     0
#> 6273eb52991d528cf5281215f79b40a3      0      0      0      0      0     0     0
#> dc4bcbe74986e44ec3a040b866dfda9e      0      0      0      0     22     0     0
#> 96ec4781ee61f4083792f252f1b6cd3a      0      0      0      0      0     0     0
#> 9657d52acca51e701e7a3b7cedcfcc12      0      0      0      0      0     0     0
#> f65f183293d4710041057ab98f2c94ed      0      0      0      0      0     0     0
#> 8f6b3783b31a2dd8640de3df92dbb9c4      0      0      0      0      0     0     0
#> 0d90eb4842c9ad9126ed33c538c768ac      0      0      0      0      0     0     0
#> 8d3ee46aaaf7728594014be61c649347      0      0      0     30      0     0     0
#> e011a95fda1f1f2bafc4e5b83b52d989      0      0      0      0      0     0     0
#> ef8ea65ba4c2e2da3ba9556a87d2385d      0      0     13      0      0     0     0
#> dbbb92eb3a8a3feaf96092d46f2b6a45      0      0      0      0     21     0    77
#> 06cb03cddd8a12dbfaf27e8984b395c4     88      0     99    158    477     0   188
#> cdf134bb501c54d997945dabd44e35ff      0      0      0      0     19     0     6
#> f3c4561cc01a3b45ee3f134d709c6b94      0      0      0      0      0     0     0
#> b121442fd1eefe78ce4f4602aac1842e      0      0      0      0      0     0     0
#> ba8707aa88a2a62939d9338189b2a733      0      0      0      0      0     0     0
#> ecf7289628d49f9244469337ee1622cb      0      0     36      0     64     0     0
#> 9f06b400b93f859b11874b0bd7f09295      3      0      0      0      0     0     0
#> f775be5cf7c0dd2177006a44913503d5      0      0      0      0     25     0     0
#> e87ceb10f2becb2a6c4cf0692b1923ee      0      0     31     32     69     0    53
#> 036574bb8e05f3b9a11ab7a4090d6950      0      0      0      0      0     0     0
#> 0267ba9db97e00ae01fc2cccbb6263e0      0      0      0      0      0     0     0
#> be4f701c1dbfbda16c5150af29430141      0      0    126     54     39    29    81
#> 655a66e6566798ade91f73b2942c0b63      0      0      0      0      0     0     0
#> cd79f2696509d0ad906f5addba264079      0      0      0      0      0     0     0
#> 07ea3284bda2bf7daf1251ad2744e4c4      0      0      0     88     64     0     0
#> e3870da9de4a2b18a7ae9e9cae208e3e      0      0      0      0      0     0     0
#> e00bd59c1ea0ca9b1cc5f976801d4e1b      0      0      0      0      0     0     0
#> b39b381d0f533c9923c3b898ffab5614      0      0      0      0      0     0     0
#> e2959f0552fa6df5707d9ec8b57816d6      0      0      0     41      0     0     0
#> 53157acd1642f973842d01e1cd287421      0      0      0      0      0     0     7
#> 5991790d1c0aa8b73b2d7f906449a1f9      0      0      0      0      0     0     0
#> adee3552c211d9e33aaa28e82532bdea      0      0      0      0      0     0     0
#> 0089040d041888e3ae4aa4e0010a0e36      0      0      0      0      0     0     7
#> 8e876e13d051992c4c8616cf2db7d73a      0      0      0      0      0     0     0
#> 179abbc5c1f1549ffd0d7ae6963e3783      0      0      0      0      0     0     0
#> 340d095f9a9fb445582b60e1282d985d      0      0      0      0     22     0     0
#> 3d7e5432d8eaa3fa56ab7138605ea43d      0      0      0      0      0     0     0
#> d4b233740db907d69035915e8657fa8d      0      0      0      0     22     0    10
#> b22119df3b1de76745abaac585b646f8      0      0      0      0      0     0     0
#> 64761a24fc66e7c93145217dd153b1ee      0      0      0      0      0     0     0
#> eb6d7f9e0601207fdcfc4a88947d513c      0      0      0      0      0     0     0
#> 0d20271e3e2dbc68e1de868599a9a1ae      0      0      0      0      0     0     0
#> 369ea22772cb774950c830d0af7e8235      0      0      0      0     12     0     0
#> cae63066cfe088a49fcac50f119c9e2d      0      0      0      0      0     0     0
#> d9c57ec20ed90e2cebf6cacc4f9a5511      0      0      0      0      0     0     0
#> f5828cc3078b4c219a44c67e96e5d455      0      0      0      0      0     0     0
#> 73847ac2778d4cba27fe9758ca129d6c      0      0      0      0      0     0     0
#> d78003f13712da48d0bb36fd45116b91      0      0     23     52     25     0    37
#> eb40853d986ad72c82658d939f291ff2      0      0      0      0      0     0     0
#> 5d675a3518222fc99c8cba34feaeac81      0      0      0     88     94     0   123
#> 3ba390ae11c4e76375af082eefe7e41f      0      0      0      0      0     0     0
#> 3aa66daa4e560cf62c1833697bcd3837      0      0      0      0     26     0     0
#> 6f4388510be5ea70356f2257d8c57ddf      0      0      0      0      0     0     0
#> b33390d034148edf07d8a742778ddcbd      0      0      0      0      0     0     0
#> 09225ec7044394bbd3b6811a6799e481      0      0      0      0      0     0     0
#> fe6fbf9b86e0f487911a856009fa4a1e     36      0      0      0     69     0    40
#> 0aee6329e0b5b81449980151bf40205a      0      0      0      0      0     0     0
#> 0ef66a78236f90916afc2305f3ff6c80      0      0      0      0     39     0     0
#> fb1eb41359005c2e3e308fb668f07220      0      0      0      0     36     0     0
#> cb59c57800227006aead644df42b8986      0      0      0      0      0     0     0
#> 138da9dc244c3292536cacb8df79467b      0      0      0      0      0     0     0
#> 72f92361ec0696cb8df985b93d277a25     25      0     39     38     26     8     0
#> 4d5ab121763879fd59e8f00d35473540      0      0      0      0      0     0     0
#> 20032b2b586fee22805eccebe3059170      0      0      0      0     48     0    18
#> 710324ccc8cafa6fc9b38baf8299d8aa      0      0      0      0     44     0     0
#> 523e6b60340f4a215877b1c336822e50      0      0      0      0      0     0     0
#> 4501f98c83484fd902be11383c8d032b      0      0      0      0      0     0     0
#> d4cbf23757aad57c91b058cbb05566e4      0      0      0      0      0     0     0
#> de7964d5c22112693ed5796bb9f04b0b      0      0      0      0      0     0   112
#> f4a283be42989b920b83ea73ec72b25b      0      0      0      0      0     0     0
#> 2f3111649d0850ab86799f8bbadc28f8      0      0      0      0      0     0    17
#> b8acdd45cf33ab75fe4c0851a1787dc5      0      0      0      0      0     0     0
#> 5bb7fec39744a6d4e4d61596dabc6886      0      0      0     15      0     0     0
#> e8eddfb37c4215c86813163a27f202a5      0      0      0      0      0     0    15
#> 826e188c911c0846fd927300127c80c0      0      0     24      0      0    10     0
#> ff1fe599c44074ee6bdf919a119df93c      0      0      0      0      0     0     0
#> aff00055d3876f624d67fb80bdb59934      0      0      0      0      0     0     0
#> 1beded194f3ae4b57715b45e50adef46      0      0      0      0      0     0     4
#> 908d91feba7c1f1d8cc8a7145173d942      0      0      0      0      0     0     0
#> e7095743695f6d4337d01ceb9e29e224      0      0      0      0      9     0     5
#> 4cd420962855afdf0bc69bb268b846c7      0     46     40     77    202    36    74
#> 430ba5bea0fc08f65b69497bfbeb5c4b      0      0      0      0      0     0     0
#> e12cb2f0d2ce4401413176007dcd4158      0      0      0      0     73     0    40
#> ec615213f6c7febb65c9f2e45e438111      0      0      0      0      0     0     0
#> 412dd52cd54d0244ba8fc42ad140acda      0      0      0      0      0     0     0
#> 89192d9f779341d313f52aad0dc6a06f      0      0      0      0      0     0     0
#> ed543f92c61a06fc4e703655ce910a2f      0      0      0      0      0     0     0
#> fea1a70b0ef0f72e3bf27d8b295331a8      0      0      0      0      0     0     0
#> 6afe805b36f5aa29e7f2096aac42ace1      0      0      0      0      0     0     0
#> d9cc5a46831f24ec9192ae82dbc38b0a      0      0      0      0      0     0     0
#> 0d3eb38eb502c4ad13bd5a7060295411      0      0      0      0      0     0     0
#> 69a0ff5538a5490e291c75ce70325066      0      0      0      0      0     0     0
#> f91d545c815effa3c5afa83a8ad3397d      0      0      5      0      0     0     0
#> e6b032dea8e13c34af958fc48af09b58      0      0      0      0      0     0     0
#> 84a1b89045ab6285328f366224f3ec87      0      0      0      0      0     0     0
#> 4a6f0076e4ac71474b3c364294ad3b33      0      0      0      0      0     0     0
#> 7bdd39e4e3ccc7000862b487510e67d6      0      0      0      0      0     0     0
#> 998aedf12cdf054f5e504a2e21d20498      0      0      0     23      0     0     0
#> 567c2640108006c8eca081ad4180e346      0      0      0     27      0     0     0
#> e50a6625b9176b6856764b213013150f      0      0      0      0      0     0    31
#> 8d5a5ab7c784b0325a4a86a16fbba80e      0      0      0      0      0     0     0
#> f58ae28a309fa3d5898d2eb46d9eb9e0      0      0      0      0      0     0     0
#> 7ac65a8fd21111b68e7c2d6205112a83      0      0      0      0      0     0     0
#> d7add64c8700f5a8cadf40c8c7f5ff67      0      0      0     13     14     0    16
#> 6bad9d981d3cf0f90262547180dc6371      0      0      0      0      0     0     0
#> 149c0e63c83113e1246c0e041de67d9c      0      0      0      0      0     0     0
#> 6f196fae54ea90b5c11627c503c98940      0      0      0      0      0     0     0
#> b67c995d842e5b4b92fa0ce052603ef1      0      0      0      0     43     0     0
#> f44ffae67c7297e2c696b75fbd4e0f4c      0      0      0      0      0     0     0
#> 33997912cd952feaec18fbd11b31acac      0      0      0      0     21     0     0
#> 0a035423c4bf71a17566549f9b69448d      0      0      0      0      0     0     0
#> f65a6f1f1e5a4bc5c733ba364c0096d3      0      0     26      0     23     0     0
#> f15389a9f245465c33db3de72aeb17ba      0      0     19      0     29     0    27
#> cd2b676d3c0785a60359c0d416701ee4      0      0      0      0      0     0     0
#> 70ba2a378ebcd9a8f57b024a561cc6b1      0      0      0      0      0     0     0
#> 0273707384afdd41e9d022b539094211      0      0      0      0      0     0     0
#> f9ac98e876ecdb504202388d7ccdda5f      0      0      0      0      0     0     0
#> bdd9b87b782b8f935ecac62a9824e6a1      0      0      0      0      0     0     0
#> 1f12b97495668f1fc378a349546d1806      0      0      0      0      0     0     0
#> 83f1b4d776f6e08a6ba6f86634a1e5a9      0      0      0      0      0     0     0
#> 49f614420449f6c7ecc8c9f4e137ad05      0      0      0      0     41     0     0
#> aba5aa00a1b262f27517eca64ee3333f      0      0      0     12      0     0     0
#> 8971980484eeabeef00145b93b699e30      0      0      0      0     38     0     0
#> c1347ae54b0fc1dd385a1ccf3b8061ce      0      0      0      0      0     0     0
#> e9fe19783ebcc7e43a3a07ca8225db62      0      0      0      0      0     0     0
#> 7c523a5f0a7c736efd82afef84137901      0      0      0      0      0     0     0
#> e13c07ec336c89db5183cc6e9e57ee6e      0      0      0      0      0     0     0
#> a017ce431ccadea609b501d46e64e55b      0      0      0      0      0     0     0
#> fadf38d03d4b944f5782a759591ec3ae      0      0      0     49      0     0     0
#> 41bfda1c7f32205246e73a1e3e3ec4f4      0      0      0      0      0     0     0
#> 9e4f3e697b70758a7a996f3dd4533650      0      0      0      0      0     0    40
#> 39c4f6caa7510ed108f3ebcf72f79edf      0      0      0      0      0     0     0
#> 4aa083fe7c1b0639dbc48ec649e68a31      0      0      0      0      0     0     0
#> fd2ac2e14b72b0e22a052e3ab4580471      0      0      0      0     16     0     0
#> 4ee10294799bcfeab0c8ee567d5d9608      0      0      0      0      0     0     0
#> 46869f12c31de7a05781e930c25a1ccd      0      0      0      0      0     0     0
#> 5e2548d6a970645d0ca943a72d207489      0      0      0     39      0     0     0
#> 0538e1beed69f098c3f602130ee11a1c      0      0      0      0      0     0     0
#> 991c563278c4e8c91ffebce3560fe51b      0      0      0      0      0     0     0
#> ce05f58164e61517198f7ab0e3304cb5      0      0      0     12     40     0    22
#> fdbd0fed21bae91c05b2da268a025a89      0      0      0      0     26     0     0
#> 7412d11595d77f26792071cd28a26760      0      0      0      0      0     0     0
#> 4b087e13957eb38128f63e1c99df23c3      0      0      0      0     11     0     0
#> 310a6f0e6841db2ce75005cac942a972      0      0      0      0      0     0     0
#> 3df0fada803ac582802b47d4f4efe0c1      0      0      0      7      0     0     0
#> 3cd4e34c02ffdd1446df65bb54278c57      0      0      0      0      0     0     0
#> b14a1dea93554257791f16c2e94906b0      0      0      0      0      0     0     0
#> 902a4435fc3af101e8a9b5f51550fa4d      0      0      0      0      0     0     0
#> 1995185b31349529f58964c2e5a3becd      0      0      0      0      0     0     0
#> 69b247c565afdeacee0240482b1ce536      0      0      0      0      0     0     0
#> 95ff1b9528acca018e9c4697a8625504      0      0      0     33      0     0     0
#> 132d26d3ce6ac2ab25ccb4dd171879d7      0      0     98      0      0     0     0
#> bf1a6e4f94ee8ba60cf79089139641b6      0      0     13      0      0     0     0
#> 704189aa09f64b992b3a02ee001dc706      0      0      0      0      0     0    21
#> dea5ea27839ee67a31d18ba0bf0f7f65      0      0      0     19      0     0    11
#> 3840a4fff0250232f8570cae7c3514ba      0      0      0      6      0     0     0
#> 31ecabc032997aaf41dd0d8def3ec93f      0      0      0      4      0     0     7
#> 024a0e0348ed9abe9ae1c48b1b7ba09e      0      0      0      0      0     0    13
#> fc1d9940419420113ce3fbfacc8d703a      0      0      0      0     51    41    52
#> 0770fbf28cd307a08bf557af7ffecfd3      0      0      0      0      0     0     0
#> 883cbe8e4d1d98f468355e8ea974d222      0      0      0      0     29     0     0
#> cb6d6561383fb68f43f92a88a607749d      0      3      0      0      0     0     0
#> e07e0f1a9fbc696a9fe4f9b7739ae45b      0      0      0      0      0     0     0
#> 73d8d8f999acfa13cba2277f7526897e      0      0      0      0     45     0     0
#> 5218274c354d41965abbf4b0751d6f3a      0      0      0      0      0     0     0
#> 5c31570f706752fa4f7851aca5b5c291      2      3      0      0     15     0     0
#> 6a5e71c8b86a52ab1dd5e61164103940     13      0     23     40     71     6     0
#> 25641b0b803c7d8493ebc02cb104bf3d      0      0      0      0      3     0     0
#> 8844d70dd1fdb96f2a7faa65bccf78c6      0      0      0      0      0     0    41
#> 371b855afa8bed030e86ec78241b01a3      6      0      0     11     13     0     4
#> d0c180a47378d414a1d088ac14db236c      0      0      0      0     23     0     0
#> a072a588b3940fab3b227e4901c40f7f      0      0      0      0      0     0     0
#> 977a3523051e6b23da6221dbde795b34      0      0      0      0      0     0     0
#> 1c3adc87d2e93b08e576953508eb680d      0      0      0      0      0     0     0
#> 74eb0ddf379b00c5032b18fb9fbf1f32     14      0     28      0     55    11    35
#> c335b34a15db7d532f511e6c952abe62      0      0      0     24     19     7     0
#> 0c78285397fb853a0e39069213510c54      0      0      0      0      0     0     0
#> 779dec6de9c6a8b7df28fc8ec3c75bee      0      0      7     12      0     0     0
#> 5526624b473bb1f0541dc6fb95c98678      0      0      0      0      0     0     2
#> 649e7fa7b9ce5c5472d34675ed643639      0      0      0     18     13     0     0
#> a9c5addeea542953d1bae4d30d599b02      0      0      0     33     67     0    36
#> 68291fb3b2558d438da4810961c8d53a      0      0      0      0      0     2     0
#> 9c5d92bfefc7a190b93d5b9a6a77697b      0      0     28     34     60     0    46
#> 3ddab42d701754625cab0c236c5b6e6d      0      0      6      6     10     0     9
#> 287e840f9ecaef80d7a4762c8d6bf401      0      0      0      0     46     0    32
#> 885f411fa54c5576148678f665bb053f      0      0      0      0      0     0     0
#> b5709a57df9cbc8a741bb78a36f29c06     62     15    104    186    326     0   163
#> 26698cc60893017ef0b1ab173688d3de      0      0     17     33     28    31    16
#> db64f0028b6518839ab7129c4b7ef72d    107     22    140    314    502    53   275
#> f46fb43262da299ec7f0084298cafa50      0      0      0      0      0     0     0
#> c45e7087de669be270881fc8ace0fad8     19      0     76    111     81    31    97
#> 854a20d5019388e81b1e78419eb017fd      0      0      0     60      0     0    86
#> 9a2325ea015a8c5944047e3349cc1565      0      0     33     49     57     0    43
#> c5305ae3c6adc90814d70ab48f4f071a      0      0      0     25     21     0     0
#> e230123ec3ce8654942ee1a000a80010     24      0     79    103    124    22    91
#> 14e3c60cb885509667670151bc4a1df0      0      0      0      0      0     0     0
#> 7fd7b7bece92bbbfcceb2dcf1e835a8e      0      0      0      0      0     0     0
#> c098d1fb769dd975c2332dfb70da489e      0      0      0     40      0     0     0
#> 95dcdf1eb06b252e67dcdbdaffc505e8      0      0      0      0      0     0     0
#> c5399d0f9eda8f0815fd28a257952496     30      0      0     23    160     0     0
#> c7388cc314e992819b8f97ea35280982      0      0      0      0     34     0     0
#> ae8facc45d372f7108bd9499f74d1bfc      0      0      0     19      0     0     0
#> 16caf6f5bc653fb6ba668eb821076a21      0      0      0      0      0     0     0
#> e7e283b72ee15bcf0877c5eeb6137ebc      0      0      0      0      0     0     0
#> 1b3a779b384c73b33207dccd14b14fe6      0      0      0      0      0     0     0
#> c8d1afdf7b395ffc22193c5dc2fd5dc3      0      0      0      0      0     0     0
#> 0e3ffc2a68378086e2480831acd1a025      0      0      0      0      0     0     0
#> 0ddf7339b608a340252b33f63fcc0019      0      0     18     17     34     0    37
#> 5fd4b2f7dee3fb0bde6f856ae266dbcb      0      0      0     19      0     0    12
#> e2f4786bc66410c07bd8b96b5d0e2d92      0      0      0     20    144     0    25
#> c30aae8f5216c9536492d0beae74950f     15     11    107    189    262    28   238
#> 9087ecee0806a80861788a17dd7e6809      0      0      0      0     18     0     0
#> c29bc08418318360592b2b2717376a0a      0      0      0     15      0     0     0
#> f791c7db1a78d35e84803b5e03c01965      0      0      0      0     26     0     0
#> 44b44fc4055a58c01f6d1e9de5e0d1d9      0      0    100    161     79     0     0
#> 342df9401103d274851d3bd80704508d     15      0      0      0      0     0     0
#> 0590e6672bb2c55c0ff54999a2b6dcfe      0      0      0     43      0     0     0
#> d4e9147608bf883fdbbf9d23b17028d1      0      0      0      0     36     0     0
#> db24794eb8667aed224401386223f671     39     13     18     68     53    23    66
#> fe41e6fd9770030a21418b288930710f      0      0      0      0      0     0     0
#> 699f18f3b155bff97de941f2cfbbbadc      0      0      0     39     44     0     6
#> 13870cf319997750e3b96bbb0b2f1aea      0      0      0      9      0     0     0
#> 7ca4fb1b1740406241acb862445cf311     20      0      0     66     60    10    19
#> 9c0569fb95451d009fd3754303eaed12      0      0      0      0      0     0     0
#> 1ac26a3c12288d8fbadd6aa61366b1b7      0      0      0      0      0     0     0
#> 89f3ba99552605501b8acb816d288afe      0      0      0      0      0     0     0
#> b4051b5227746f04c1feb9d8979299ff      0      0      0      0      0     0     0
#> 087d358aaf77f81aaf50a623828865bb      0      0      0      0      0     0     0
#> e451c4b3ff1c3dabeee296529b3a7a07      0      0      7      0     28     0    11
#> 012871784f77c9dca7d446ddef2048d1      0      0      0      0      0     0     0
#> 95c6ce2788783bf2f83343deba23b8bd      0      0      0      0      8     0     0
#> 213227f244d97553a9a14c4e3fd04c98      0      0      5      0      0     0    15
#> c0e5e55c5cfbd576e7c085f1e4b69638      0      0      0      0      0     0     0
#> 51e98451d7e3f46127481ec5490f1c91      0      0      0      0      0     0     0
#> c75285a6c61d8fce20d823769c25c12a      0      0      0      0      0     0     0
#> 177966f7619eb24fdbce2266f5bed39a      0      0      0      0      0     0     0
#> b89347f0534679f890628364d2c7b4fe      0      0      0      0     20     0     0
#> 8dff1cb12c0beb5bc7b2546480a6a6c7      0      0      0      0     38     0     0
#> 84f2777cfcb0fbed6281917fcc765022      0      0      0      0     10     0     0
#> a5e4d106a2ab947984fdef187996a0b9      0      0      0      0      0     0    22
#> b6cbf6077b64acbb718971a1ffc97a27      0      0      0      0      0     0     0
#> 27e29495fda424728fd6a80b83a16772      0      0      0     32     42     0    33
#> 63c44ab3725dd0b803d43a29db9fd7fd      0      0      0     21     43     0     0
#> 9c931ab5bf18ef9b10996d59d292e6f5      0      0     37      0    117     0    37
#> 12fd8efb82c3d4909543f6db9853134e      0      0      0      0     35     0    12
#> a8a60c20f7ab7d0f2dadd9f8858c19be      0      0      0      0      0     0     0
#> d57e09eb1645f1abba6074d3bd7d4c52      0      0      0      2      0     0     0
#> 7e0f2b982ca2d35b6ac2cf251d78d3c4      0      0      8      0      0     0     0
#> 3e2e8e938a68a0318bcf7145ad5ed4a9      0      0      0      0      0     0     0
#> f2117b49bf57d9451bbef6de5b762483      0      0      0      0      0     0     0
#> be24dfe003725a0b0f28f4136106faf8      0      0      0     16     14     0     0
#> e63301d7d0af55d5c43db051dc51522f      0      5      0      0     36     0    34
#> a6e2e5501507562f2174f65e606b2d43      9      0     46     23     36     0    19
#> a06aa93cb50b9a39941a11c7242dde96      0      0      0      4      4     0     0
#> 6a0a8abb36aa3bfdee14edbe97e9d781      0      0      0      0      0     0     9
#> 4d3def5ffad0c392381591440bfc3668      0      0      0      0      0     0     0
#> ce437b81e7cf5ce986f7240c26464b0a      0      0      0      0      0     0     0
#> e57788f9d2fae76ed9ca566ad4d475a8      0      0      0      0      0     0     0
#> 8c525d0f63137bd068b75c57f5d04c09      0      0      0      0      0     0     0
#> 6d670aa40dafa152c2d59cf2be9bef74      0      0      0      0      0     0     0
#> e08f8e9e9ffb10a702f2ef617f70eecf      0      0      0      0      0     0     0
#> 7def1c807ca0767ccdfe5a6001dc51f1      0      0      0      0      0     0     0
#> 3b9a17275147adf1dc8c42eb00a06a36      0      0      0      0      0     0     0
#> 24be0b7e4c85ac53b09fddbca916ab96      0      0      0      7     20     0     0
#> f521a47d3bb4acfec6b261f880b4eea8      0      0      0      0      0     0     0
#> 91c24c9b6ffbedcc257d164176060fe0      0      0      0      0      0     0     0
#> 9b45f210647c5fe71b005eb9e25299a5      0      0      0      0      0     0     3
#> 974077d25fa6b32a7af22670b142a8f1      0      0      0      0      0     0     0
#> 203aa2535e54a11e4c48af4ab946686b      0      0      0      0     18     0     3
#> 9064afe569448558eb08e5bf5e458bce      0      0      0      0      0     0     0
#> fad432eeeff7f1d923e6d60aafe18dcd      0      0      0      0      0     0     0
#> a29db0e5f4a04c5c6937745ccccb167b      5      0     20     13     86     0    34
#> f179c10f8e86873c90ceb14b974a114e      0      0      0      0      0     0     0
#> a23e8588a2240d7af6f65519158c902c      0      0      0      0      0     0    10
#> 52034e0da6d5e7962369437956984e18      0      0      0      0      0     0     0
#> 522d20c27cf4237b60f36ffdca056c19      0      0      0      0      0     0     0
#> 957615aefb421e6ad4728e610b2bdfa6      0      0      0      0      0     0     0
#> 5ccc74e70b8a82a6eebea3385738ece3      2      0      3      0     25     0     0
#> e5924fd11a48419c6ec0e2a7cb4e43ee      0      0      0      0      5     0     0
#> eef7be4b2531cb7e909b58e91af5640d      0      0      0      0      0     0     0
#> 4416cbc2026260fa4e2ed0a8feff1443      0      0      0      0      0     0     0
#> 03f6091f5aacb69f06fc925acd1c9506      0      0      0      9      0     0     0
#> f655e5ddf5074d5e71bfa583926dec30      0      0     11      8     29     5     0
#> 9f733d39e588c7310499dd444cdd7812      0      0      0      0      0     0     0
#> 2b413de5e05a1385bab229f50879628f      0      0      0      0      0     0     0
#> 178b596f2fe5cb6cefc5513917988d84      0      0      0      0     14     0     0
#> aa94b014c726cf40983045bd606241c1      0      0      0      5      2     0     0
#> 55310b3aac4868c4e33b0a9492d7609a      0      0      0      0      0     0     0
#> c5d8d828df0e1baf5c2db05a6b8c4476      0      0      0      0      0     0     0
#> 6ceb26f60606233a7728d3457af84380      0      0      0      0     14     0     0
#> 1ff83868d1c9ecf52d9b776e0b3fd896      0      0      0      0      0     0     0
#> 3edfccfb9f81b340854e84ee5102a4c1      0      0      0      0     19     0    13
#> 963dc552c551872800f2fcbc5663e5fd      0      0      0      0      0     0     0
#> 3fda8c20e5da9fbc3fb7c5c864559019      0      0      6     11     17     0    20
#> 6fd88fa7e883bdead3c5ba89da853e8d      0      0      0      0      0     0     0
#> 35826fe826b7a9f202e8bef7b39d9a47      0      0      0      0     22     0    26
#> cc514f614ba146f489bf62b5a56074c5     25     13     20     38     48    25    36
#> 26260de3858081ce55e6f220b3ba25df      0      0     12     43      0     6    25
#> d6727a85782ffd2a069dccc23adaf521     39      0     55    121    265     0   100
#> eb1a74c2ca8831ef945775de0827ddda      0      0      0      7      6     0     0
#> 4ed3218615bbc1ed4da82d9806e4786b      0      0      8     10     11     0    14
#> 8a3367e36ed4e416d5a77a7946cf6fc5      3      0      0      0     16     0     0
#> 56f45eba0496ed29fb7171041a2a950d      0      0     13      0     32     0    26
#> 9d2782f19e60ec0cf72e1dc46b1d1766      0      0      0      0      0     0     0
#> 1559f8b64dc2bdc8cdc26ff8a5ffa4cb      0      0      8      0     14     0    18
#> 6d01f384a1cce140edea17dccd0f116d      0      0      0      0     12     0     0
#> 5950949b71082741302b712e97309b17      0      0      0      0      0     0     0
#> f6f238f95cc97995c9bf619b0c620b33      0      0      0     10     13     0     0
#> c7a112ffadf8009efd318625da44b18e      0      0      0      0      0     0     0
#> e79f4f30b0dacfcc50f3c1b12804ca3f      0      0      0     15      8     0    19
#> 297c825dc268f9a7b86b402edb6961b3      0      0      0      0     18     7     0
#> a1ed56ad4e86e0b4adc11baf2c4ab94b      0      0      0      0      0     0     0
#> 8afccb8f7f50c5603b471fb0c808c9e6      0      0      2      0      0     0     0
#> ce41c973bd4965a400989ec128890d55      0      0      0      0      7     0    14
#> 5d99ba689ad34e4393a905e8a725b72b      0      0      0      0      0     0     6
#> bfaa2d35e6dd13b9b32de811c6f92592      0      0      0      0      0     0     2
#> 0c544d5f3a49a010971a96d8c39e4905      0      0      0      8      0     0     0
#> bab3f08e1d75f42987012cae451639ff     25     12    117    202    389    24   244
#> fa97a81b8c33569e2aa1fc9a3137b435      0      0     17     43      0     0    29
#> a4d08886c3c6c62efd5678339a07828b      0      0      5     28      0     0    16
#> 11ed0972601847ec168473a6bc4eeba7      0      0     10      0      0     0    12
#> a6e4d88aa4b40402b1f94760d86b729c      0      0      0      0      0     0     9
#> 86deaff47d33946ede1b7c25c6bbe9f9      0      0      0      0      0     0     0
#> a1d3f1f8c151eaeb5cf12179861e0c87      0      0      0      0      0     0     5
#> 3980891c22ccbb387b65f9e437b2f94f      0      0      0      3      0     0     0
#> fc8ec321e0c2fc2a6de9abe98256b3f5      0      0      0      0     25     0     0
#> 84c11e6ffb29ebc2192f4b049b7dbddb      0      0      0      0      0     0     0
#> 82104585f0b2617c13eb69f0bbdf7d5a      0      0      0      0      5     0     0
#> 61cd4918428d3b36910b9cb2fc45946a      0      0      0      0      0     0     0
#> 0ab9a976935482bce5f8af4cc6685147      0      0      0      0      0     0    16
#> 87be3450216b45c53851fd173582a2d0      0      0      0     15      0     0     0
#> 93a51bb503fadf90121f518db980e158      0      0      0     13      0     0     0
#> b70cd843161996e075a853900e699d0f      0      0     26     95      0     0    37
#> 154b8be9eea158f93d27d314d8e0e2c8      0      0      6     41      0     2    13
#> 70691c2737cf7dc7c194de819f0e6cda      0      0     16      0     40     0    21
#> 2ae1d80ead92f5991f58da14d000acd6      0      0      0      0      9     0     0
#> 44fa299492152518fec34b9fa4554648      0      0      0      0      0     0     0
#> e15c8f9df87f363b0dbbf49de0bed249      0      0      0      0     18     0     0
#> 73689cad6eb4e86d9d0bb48c0da24463      0      0     11      0     21     0     0
#> 234e15e82d64f2e1d61802fa8ea1998a      0      0      0      0     49     0    29
#> 88c76ddb1b24acffe65d8128e2f48244      0      0      0      0      0     0     0
#> 5544139a16716cc6ef48b2c16c9613e9      0      0      0      0      7     0     0
#> 5b7fe68483e701059142cb531397b547      0      0      0      0      9     0     0
#> 0a0875a0acf4a28d1bc1c4d1d11aaa64      0      0      0      0     19     0     0
#> 6b47ec9d1c70f07d16834e9a4e54aa70      0      0      0      0      0     0     0
#> f8752b4824d239cfb19f3216cc34941b      0      0      0      0     61     0     0
#> ec8a0768cd2fb426d7e6bfc2c4f50e5b      0      0      0      0      0     0     0
#> 8a6b88f9515bf2dfac09e6a69045a4fe      0      0      0      0      0     0     0
#> f5151c67af19c8471743ef3e2914f6b1      0      0      0      0      0     0     0
#> c3842619fcbde52781c512c2561fca0b      0      0      0      0      0     0     0
#> b3b00faa9dee9fa6b3b33744f51da00e      0      0      0      0      0     0     6
#> c4f04e28c6267ac6c5b05d403e8f91c3      0      0      0      0      0     0     0
#> 7fc06f107a6ddf6486ef2d43fcf23626      0      0      0      0      0     0     5
#> 8d1dffcd377b338d2a0388637507b592      0      0      0      0      0     0     9
#> 15410d762889d51818986bfcdb5d8965      0      0      2      0      0     0     0
#> 48f68e44e258de63575b8bf6cf565d52      0      0      0     16      0     0    19
#> 4712befe3685661e0fc0968859c369d5      0      0      0      0      0     0     0
#> be9a3727a0422aea147100370d046fa8     20      0     13     93    157     0    32
#> 3c5942781761830f7b5406f551574424      0      0      0      0      0     0     0
#> 291d9f4d3d3c5cc3b0a1804770ed1e77      0      0      0      8      0     0     0
#> 784764c519a64dc9a147025a23222e93      0      0      0     15      0     2     8
#> b9e955ea5254dd68b54775aad2a67829      0      0      0      0      0     0     2
#> 8197f1594e8a17f3fc989fea64f54b8f      0      0      0      0      8     0     0
#> bf4a4f2383a16d2ecf27bf8e162b9ba4      0      0      0      0      0     0     0
#> 4e7d22d793d30bdedbd32a5b73b564b8     13      0     23     35     20     9    32
#> e88959b739c8cceb94b94b7887dd1df3      0      0      0      0      0     0     0
#> 5403445694339e2ae5ecac06d9f3ff47      8      0      0      0     19     0     0
#> ba3943fb075bbf63fb699638292826a1      0      0      0      0      0     0     0
#> 67b2d836db093cead5a8286142669613      0      0      0      0      0     0     0
#> 6c3a0ea164d3eb098dd2baac5519b38f      0      0      0      0      0     0     0
#> c643a776c36bcc00dfb5d790e1ea3dfe      0      0      0      0      0     0     0
#> d53cc81a9d57a64d194c876b76bd66e0      0      0      0      0      0     0    19
#> 7dc00d40971bb454b39fdd9b2dd9aa25      0      0      0      0      0     0     0
#> b4d0cec8543e940be1b5096370fdfb80      0      0      0      0      8     0     0
#> 6c54ceb94bc5e8fbd65224cffc60990a      0     12     14     46     22     0     0
#> b12ee9ca5a2320957213a114788cf2b7      0      0      0      0      0     0     0
#> 5e040c2397f959e7f0dfc6ec2855b812      0      0      0      0      9     0     0
#> 3105a594bba2604ae5ed19f5e0495d2e      0      0      0     29      0     0     0
#> d75c040bdaa1cca96cd0655c9af7ab28      4      0      0      0      0     0     0
#> 9b909edaebabc90b9c1a703626e6bd0c      0      0      0      0      0     0     0
#> a6593da81035cf121dbe46346a4ea960      0      0      0      0      0     0     0
#> 60e4c88400a750a8122b97476313a67e      0      0      0      0      0     0     0
#> 7854a78e3b881cee1c63ea0a50ba29ce      0      0      0      0      0     0     0
#> 1796734b56d45bdc997ccc09ab8731d1      0      0      0     18     48     0    33
#> 8243cd669ae7327f4326be78429ab05a      0      4      0     27     31     0    33
#> ebb93dc32c63402289114859d476cc5b      0      0      0      0      0     0     0
#> 8a01af5a531e851e2c50fc01131cda33      5      0      9     17     36     0    19
#> 03d2cc3f7b2bc352752d39e3738d0af9      0      0      0      0      0     0     0
#> ec378af6e32960b11b5f7391da4ae4fb      0      0      0      0      0     0     0
#> b407f2de868c064f65eb5ad8f82a8d0d      0      0      2      0      0     0     0
#> 73a95b6ba34fcd260989d848fe15d377      0      0      0      2      0     0     0
#> 3abb68202d15781590a2d1f9c5b37e15      0      0      6     25     39     0     0
#> 77923ef8e9ed01c9b6a8ac85eb028163      0      0      0      0      0     0     0
#> 9981b8fd20f400e052527f84f5db1b66      0      0      0      0      0     0     0
#> 597bdec49adf7675912f9fba4e262e70      0      0      0      0      0     0     0
#> 0fd3126118b13fb187bc55476a67c95a      0      0      0      6      0     0     0
#> 2a0c511ecc2513c42fe624e8e99d07da      0      0      0      0      0     0     0
#> d2912fea2dcf2e5b97e56c6a7d24295d      0      0      0      0      0     4     0
#> 2976af6aaf4261b6c274fcc3229a4d65      0      0      0      0      0     0     4
#> f097c45ce419da576aa111c98a4c7272      0      0      8      0      0     0     0
#> 3abc8260d34008f263e232ea25426c3f      0      0      0      0      0     0     0
#> d7e5b6432dc178aaeed9c81e0ec57afb      2      0      0      4      5     0     5
#> dfe71ce78ddfdbc7b453d17d628365f2     12      2     18     28     58     0    51
#> 289b18a86084273ce3b4038bf37e2840      0      0      0      0      0     0     0
#> 844a67032c52dc49019bf5dff248bf9c      0      0      0      4      0     0     0
#> 38cb3e3ce5ac5aed7dec50247ecfb267      0      0      0      0     36     0     0
#> 74e88e71e9929b1ddfcafccfab3b7bba      4      0      0      0      0     0     0
#> 66bc0737fc7fa1463ac7189ab7073998      0      0      0     16      0     0     0
#> b5a7aa588c9cae6ae8ca2f7d87ac53e8      7      0      0      0      0     0     0
#> db0540920457b21b170947aac9137771      0      0      0      0      0     0     0
#> ebbaf260510ab16ac6b8c2b9866d05b9      0      0     28     26     33     0     0
#> c1f461b67419d4a80a61db0867da9ed6      0      0      0     28     49     0    16
#> 0c6eba6452635eecfb865a1049299176      3      0     19      5      0     0    22
#> 7ecabb10f29a0c44eb53e9f9852fe22a      0      0      0      0      0     0     0
#> a33cb3037e7825920073b4dddaaa41e6      0      0      0      0     41     0     0
#> d625812c15c074c2e2b18099f97412c9      0      0      0      0      0     0     0
#> 17ef77f12e2a55717e910a5f6a625448      0      0     15     18     11     0    12
#> 7baf638fdd3f71166aafd22ecf410348      0      0      0      0      0     0     0
#> 9436cda211f60f722865cd23ef75d4ed      0      0      0      0     16     0     0
#> c8ab0b20e740eb6b35e68d58f69736e4      0      0     45      0     36     0    87
#> 5e43f2d8be71399fc8d6736754cb42ac      0      0      0     21      0     0     0
#> 19ad90ee1081a003406b02179ae63709      0      0      0      0      0     0     0
#> 0a25cdbd0660e68970c563c43a18f54a      0      0     30     26     44     0    23
#> 0aa6b210e902d8ab9d6f267407258ebb      0      0      0      0     33     0     0
#> a864c36d0d898cc5abd94d6d6f0bd5be      0      0      9      0     49     0    17
#> dea0dcb305ff0f27e3c36c7133cf933d      0      0      0      0      9     0     0
#> 463977baa69c7724fdc9407ad162714c      0      0      0      0     11     0     0
#> fcb72326c29342bcfb052953466eafc6      0      0     18     38     25     0    46
#> 5a554353f9ed728f0ab0631f29c8be2e      0      0      0      0      0     0     0
#> a9dd27b48e443415a743bf9902f515c4      0      0      0      6      0     0     0
#> 77a40bcf979c0099c52f19a3f1fe5217      0      0      0      0      0     0     0
#> b4007ac525a9b409d70426e24421cff0      0      0      0     26     31     0     0
#> 406ba18a175c4c79efb25f48a6beafe4      0      0      0      0      8     0     0
#> fe1c166755bc722cefd2e2a3b2b70672      0      0      0      0     14     0     0
#> 8afd274849ac07c0ca66f8f5f143df3c      0      0      0     29      0     0    16
#> 197301ed9d8728b5c061b99a17973161     20      0      0      0     62     0    15
#> 14779f9dc3505598ad3d799e4ee676c4      0      0      0      0      0     0     0
#> 4085c246640eee043d523ce2d799023a      0      0     37     40     19     0    58
#> bba09a403e67c1381de01c86511f1057      0      0      0      0     12     0    12
#> 392c59fee4e6b91d7654f2d3e3c65be8      0      0     28     38     38     0    29
#> d2b05ae73ad5ed296d7add9e2bb2cce4      0      0      0      0     13     0     8
#> 8306a9473dabfd5918ea869dbfcd1291      0      0     16     20     16     0    20
#> 06064c372bb624503c2905d13f7691f5      0      0     57    108    176     0   102
#> 3265692341e42dbd5cf4ed1938b906e3      0      0      0      0      0     0     0
#> da96a0797b9881a123795838b9f01140      0      0      0      0      0     0     0
#> 0bdddf10e577657becec88e42001fd16      0      0      0      0      0     0     0
#> 168999a950ebe17e6d2ea1f90e557d91      0      0      0     17     37     0     0
#> a266ffe9d38f626d106fa615da4618ef      0      0      0      0      0     0    16
#> 7706405d5bce7715f01d8b54188def42      0      0      0      0     48     0     0
#> 4b0cfc8899b466bd9992f7b22a5a06be      0      0      0      0     23     0     0
#> b81b2b4a7ca9eb2c23b170648dfaf22d     28      0      0     19     21     0    17
#> 75cf40872201f150cb5b191676f64b78      0      0     35     43     74     0    48
#> c8e5794c58de9326084a62d3282d257e      0      0      0      0      0     0     0
#> b0706a571f169c0bac7b4864052fcaa7      0      0      0      0      7     0     0
#> 0388f969e3527b2938acc676d6757bdb      0      0      0      0      0     0     0
#> 2dfe751390aedfbd9f62dbb77d1f10b4      0      0      0     13      0     0     0
#> 2548dc9697e1f9afe04a9012f40e2f0a      0      0      0      0      0     0     0
#> 2f019615d5141a4ab00b1aab98882a33      0      0     27      0    134     0    43
#> 629e4f0cf1826904ff74907b2507f58b      0      0      0      0      6     0     9
#> 774861f348cec8451102730c21a790b0      0      0      0      0     24     0     0
#> 11c0d6de28ab026a9ae17abe4efdd6bf      0      0      0      0      0     0     0
#> c11381dda0841c0d08cbfcef3a760aa3      0      0      0      0      0     0     5
#> 3c9619d4f9db4ac6ad7151a0c9d86ad3      0      0      0      0      0     0     0
#> 7e2c72204a8c45b1b0ac39e0e6047098      0      0     22      0     44     0    26
#> 8254e0f9f07709b62dd8483f3b25fad0      0      0     34      0    109     0    51
#> e78135f2886d427e678d12c446bf070e      0      0      0      4      0     0     0
#> abec938a8d53e7471f3f476d50694be9      0      0      0      0      0     0     0
#> 4963567d2979bb761291aa98e0b8595b      0      0     15     34     32     0    43
#> d4de9d1453102564f35262faa3952266      0      0      0     18      0     0     0
#> efa4a78c4b993e9d115841c8a02c795e      0      0      0      0      0     0     0
#> fcacd7424bf3543322ea3e9842d14e86      0      0      0     28      0     0     0
#> f941e15ea7f7560191d57abc4d2bedd8      0      0      0      4      0     0     0
#> ac294f5242de3c9b26e75502c48f6383      0      0      0      0      0     0     0
#> 7b6a8e3f8bfa1bcd1b5fc55db048d005      0      0      0      0      0     0     0
#> c9e776ab599fe75db286675b698c9067     23     33     34      0     66    35    33
#> ebe5cb945f6d6d849c07424591f80753      0      0      0      0      0     0     0
#> f116387f77802409103a559e4205df06      0      0     18      0     17     0    30
#> 93f9939b261cddc9a57071b42771c1d4      0      0     59      0     59     0    43
#> 971a5f48b3418449ca7cf3e190a5e831      0      0      0      0     67     0     0
#> 730de992ca06a95e7193a93146a97cf2      0      0      0      0      0     0    35
#> 8fa8b4ac69784c26901e036d13ab2332      0      0      0      0      0     0     0
#> c363c9efaa85023b488b732b97e47270      0      0      0      0      0     0    18
#> 2e2181cc893ca1e1e677e888bf8af441      0      0      0      0      6     0     0
#> 86fa662da6869634e7f90f85e9ee4933      0      0      0      0      0     0     0
#> 9a02177f2c235a6eb62984d5aab151ff      0      0      0     43      0     0     0
#> d6854be498e22b65da3c00572ee178f0      0      0      0      0     33     0     0
#> 015f364f150b230450da0a123e7bb8ab      0      0      0      0      0     0    11
#> af2c1b22758cd6653e8e841497d096e7      0      0      0      0      0     0     0
#> e7f2fab8c0e7c365768a43ef71406605      0      0      0      0     24    13     0
#> fb57cf192bce1876471f2f6181bf4b55      0      0      0      0      0     0     0
#> 6c740ce69c9f45bea684dc9e36054006      0      0      0     38      0     0     0
#> 3b27c7bf279535d8e39664d8dc1d0daf      0      0      0      0      0     0     0
#> ceea3c29ad012c2918fb6aa21a44e8cb      0      0      0      0      0     0     0
#> fa99ecdb27c1cc44fd866ee1da250c4c      0      0      0      0     36     0     0
#> 387cc83573db0e444e81c89a8a4f8f9c      0      0     53    115     93     0    75
#> c7a8d1af5f7f472671f700775de3e0a2      0      0      0      0     19     0     8
#> 98a6a84550e07ea252cfd8a0f51042d4     23      0      0     44      0    40     0
#> 08b809a68126c38db786cafe9710b013      0      0      0      0     41     0     0
#> 8db4cc0a5d2f03fce8efd2b8374ff359      0      0      0      0     10     0     0
#> 3e7c7d38423c93b43a1d9865cf0bc4a6     10      0     23     35    139     0    41
#> 6f02050760c030ed8cd1ed836697fa6b      0      0      0     21     47     0    21
#> 8fc1865cfcae32f72c1b1570299e2ef3      0      0     24     18     38     0    20
#> 38cc138fac95b006929361093d840a00      0      0     21     35     31     0    28
#> f5f8eb63e22362509e2c7c78271305f2      0      0      0     17      0     0     0
#> 0bb238aaee3ccb7dc8dfa86dbc1ba9f7      0      0      0      0     12     0     0
#> ca4546caf07edf5200cf79cea04aebcd      0      0      0      7      3     0     0
#> 893d12123eb7f623489a4f025a15ae7f      0      0      0      0      0     0     0
#> 192c5a817f48368339e519e0a0a907da      0      0      0      0      0     0     0
#> 3467185ca92b7bec2523214605adabf0      0      0      0      0      0     0     0
#> a2aae80d8fb922fe1074d243b59ce237     12      0     44     94     43     6    48
#> f11233f16d4177d38ab1c0bfc649e056      0      0      0      0      0     0     0
#> b219ca7793ab9f5e441f2a68bc76152c      0      0      0      0     21     0    24
#> e6e1848b51e705e55b3805703f088311      0      0      0      0      0     0     0
#> bd074b9bf365d11468bd76268d63b15b      0      0     15     60     63     0    41
#> f381b492d6651c52fd09ce9351e7eb2d      0      0      0      0      0     0     0
#> d9b31e1c30facc24fa778b1d65f5c457      6      0      0     91      0     0    14
#> 14a67c8c8cc31d203eac8109443fab5f      0      0      0      0     24     0    13
#> 026ad0ea53fb202ebb770acf65705d6e      0      4      0     17     22     0    10
#> 19ba9fea04969f975d6dc8d3da20d5e6      9      0     64      0      0     0   142
#> 0fc4426ed7dccc48ce27cfff8083433d      0     13     52    530    253    50    52
#> 45ad90e179842705fa984962f131fd61      0      0      0      0      0     0     0
#> d284afab226d22fd0596221d17db23d3      0     20      0     39     70     0    61
#> 38194655ba17dbfd7fcc1ed42a57c169     21      0     40     40     70    14    52
#> e30be5268bbebe597424bdf171cc1607      0      0      0      0      0     0     0
#> cb7ccdb4956134b452a3cb78d1751c57      0      0      0     12      0     0     0
#> b585436677ac88dd4ac6239655f2a268      0      0      0      0      5     0     0
#> f5d770eff5c9b9b0c859db57d175eb34      0      0      0      0     22     0     0
#> 1b45ffbc88c788dcee2239e68901cf05      0      0      0      0      0     0     0
#> 3fd951cf7cc3a0740a1f20f8a9c6d2c2      0      0      0      9      7     0    10
#> 2027d3e90ccc4003f1954322e99db6a9      0      0     18      0     38     0    22
#> 6157cb92771830955413f808e83fc420      0      0     14     19     35     0    13
#> 6349506cac32056a214b04239bdc3a33     52      0     48     85    173    33   111
#> 2254a8d37cf537b1e35097d70b6ac91b     26     19    113    114    780     0   192
#> 38f3bc8ef8836574de0fcbeafd348b65      0      0      0     29      0     0     0
#> fec23c70ae0c70ecceea1487cd978ee5      0      0      0      0     66     0     0
#> 503d123a6c78c6fb280a6472bd61b51e      0      0      0      3      0     0     5
#> 9f33de8570c639e727465f1ee4d750f6      0      0      0      0      0     0     0
#> d68e1899fa65515135d9bd8d356b7e3d      7      0      0      0      0     0     0
#> e15b88ccac7a4bef64215e915f6e0d89      0      0      0     14     48     0     0
#> 767042b4b19d8bfaaa426dad7dbce378      0      0      0      0      0     0     0
#> 4a3c4abee502d6185abd750522566bf3      0      0      0      0      4     0     0
#> c254aa5c00b7f27d756f7bbe96b3d4ad      0      0      0      0      0     0    27
#> fe8d8af791eb6e5d4335bb7604c74e9d      0      0     31      0      0     0    17
#> d81800d7e31b69877ce1eab9858b8f3a      0      0      0      0      0     0     0
#> c55c7cb80c5700ff9aded76192543b5a      0      0      0      0      0     0     0
#> 91dee7b3b90bbf2e53b2472b9bbde0e9      0      0      0      0      0     0     0
#> 6b28df61cb5dee74b76b1a13aaa0027d      0      0      0      0      0     0     0
#> 9a830fa26205ad91a02b5bcaebb0221a      3      0      0     38     26     5    58
#> a54b8142d2a09d058b81c377e8e286b8      0      0     18      0     25     0     0
#> 69ef3a03a4d733dbffd99fdaa1d2e62a      0      0      0      2      2     0     0
#> 083cbb131dd305a88b7d338bc371c35f      0      0      7     11     25     0     9
#> ef1a92caffb4d7cae286fdc77b0a93a6      0      0      0      0      0     0     0
#> 2cddd5dca166486f50a0828b28001e1e      0      0      0      0      0     0     0
#> bb899d1dcfd89cb200a76cb8951904f8      0      2     23      0     26     0    16
#> 7817f6851e9df03a3d992cf915310dba      0      0      0      0      0     0     0
#> 54eaf6571de743a86ccb531bee8aad7e      0      0      0      0      6     0     0
#>                                  4.1RO 2.16RI 1.2RI 1.2RO 2.10RO 4.18RI 2.10RI
#> 374a5eb6496d14e8e9d6f3d8c7e34d9d     0      0     0     0      0      0      0
#> 4c47d5cf81df1ac6e2a0560029a81578     0      0     0     0      0      0      0
#> 8fb6741a6685fad0374f85ad3e16156c     0      0     0     0      0      0      0
#> 310f3de009c95de5c938b8e811dfc96f     0      0     0     0      0      0      0
#> 7d3b39e8c0f0fdcb81b0c898a98c38f9     0      0     0     0      0      0      0
#> fc211549300b0954dbb4a4bf57e8a605     0      0     0     0      0      0      0
#> 53fad3d8c52f4b5a023ff92aee9c0f1a     0      0     0     0      0      0      0
#> 5e37342d650cde6cfc5bbe2726b12e3e     0     17     0     0     10      0      0
#> 56025f74cc8e54d4a054442d786ef096     0      0     0     0      0      0      0
#> d9e51d98f1c86a1c426552d252c944e4     0      0     0     0      0      0      0
#> bc905dcc609da83d3620f1986e63a9f7     0      0     0     0      0      0      0
#> 90009df763e86e6eca7baa54f01551e6     0      0     0     0      0      0      0
#> 6bbd12b392461134369d835715cd2765     0      0     0     0      0      0      0
#> 83c88a72c2dab171c38aadd98aa2d389     0      6     0     0      0      0      0
#> 2b93708de9f1c24ff6140814a599a3a8     0      0     0     0      0      0      0
#> 671a53a2acbf0c8c7e240e446daf0000     0      0     0     0      0      0      0
#> bfe35644221421300c0c17e00e465c05    16     40     0     0     45     66     26
#> d12bdf20f24b92aa8ec94c6c3f60497f     0      0     0     0      0      0      0
#> 3bdb72ce060bd89f107c5c2bc95399ed     0      0     0    94      0      0      0
#> e6a82be4b75aac8b928c1dfbe9ee9287     0      0     0     0      0      0      0
#> e3a63ba0a0ec40fbf270a57c180e0d3c     0     16     0     0      0      0      0
#> 2877cf981b128a18f788ba5437c1afa6     0      0     0     0      0      0      0
#> c015bc50696fe3e3d94ace01d9dcd691     0      0     0     0      0      0      0
#> eb2eeddf829080f8d9f25c6be5d507cd     0      0     0     0      0      0      0
#> b1fcda6d8df4c25c8e6d7124c076b1fe     0     17     0     0      0      0      0
#> de8bb5c39fb121fd048dfe0d639e6e8a     0      0     0     0      0      0      0
#> a234ac03223c275785dedf6df3b2aeff     0      0     0     0      0      0      0
#> 9ea1641eeed1f51fc228c653076d1d08     0      0     0     0      0      0      0
#> 694b926113eb8291fea166446890c8fe     0      9     0     0      0      0     14
#> b2af163540d17007833dc819704afb28     0      0     0     0      0      0      0
#> 6fc0ba1f8ff8f9259c9d492e59771755     0      0     0     0      0      0      0
#> 7b8c1dba059eac7ec7f3cefa203f6b2d     0      0     0     0      0      0      0
#> 0b06028193ef1ce19f5339a8ba7d4448     0      0     0     0      0      0      0
#> 56b4ec05ef3bbded515ec3dd2d573197     0      0     0     0      0      0      0
#> fc9c84f8767611251f0eb93af4ca10db     0      0     0     0      0      0      0
#> d3616bc0d8925a452274c484c73406c3     0      0    87     0      0      0      0
#> 7986b8b93f3afce920bd19a5d2229c09     0      0     0     0      0      0      0
#> bc78ee7297b6912452be0865d113f7ad     0      0     0     0      0      0      0
#> c472668b46b0d0f575108aa4170dd81f     0      0     0     0      0      0      0
#> 2658a78c5ff06760389ce2e3a548b4a2     0     37     0     0      0      0      0
#> 600c9e81fd1d1ddd5a37a5aa0f9f796d     0      0     0     0      0      0      0
#> f8d0bb97f093d8e8630ea62ac8a1bf13     0      0     0     0      0      0      0
#> 6c685b83d63d6d56a9dcaba34215ae6b     0      0     0     0      0      0      0
#> e0e7e6208280a5846a8b66f00667f8ca     0      0     0     0      0      0      0
#> 6d1aa1c1f93c423a856a5dd3b62e5c70     0      0     0     0      0      0      0
#> e553e766ea0e08f1f9fa2debe7946af6     0      0     0     0      0      0      0
#> 0f1ef43f6b23a567aba80599766d70bc     0      0     0     0      0      0      0
#> 3c4ccaa45e8d148c9a17b5d0f3720cd2     0     86    47   112    117    123    108
#> d9df6265d046a558d87d397b8a177d54     0      0     0     0      0      0      0
#> 4bbad82a39655a26a2de61b88c5af7fa     0      0     0     0      0      0      0
#> f5a6865451a3b138f4915c5ab445a0a7     0      0     0     0      0      0      0
#> 4a86fc550fa1ac37aa480d37edea875f     0      0     0     0      0      0      0
#> 02fc21e72f66adea757b355edea4a369     0      0     0     0      0      0      0
#> 25cf10f4f15dab00925cef52b9f08442     0      0     0     0      0      0      0
#> fa8c8b98d197c252368c160b43cefb98     0      0     0     0      0      0      0
#> c0b07d29aed1a875641aede53bc0a20c     0      0     0     0      0      0      0
#> 341433c1c8cc65ce47dcf32179d08fd8     0      0     0     0      0      0      0
#> 8bf175f329d16e4aa732cf2b32279df3     0      0     0     0      0      0      0
#> 062f7552b747dccab9586dff7275b5d9     0      0     0     0      0      0      0
#> 059e124ad6ba8f98c3d21db28a075ea8    63    146    63   113     98     35     70
#> 3d7a81faea3f783158352c9f739c8d17     0      0     0   138      0      0      0
#> 90be2fb2c82a98002d33ed7531e877e3     0      0     0     0      0      0      0
#> b0e263cdd1c189fc1757d8189db14be2     0     49    63    49      0     54     45
#> b962ce332fea338b05293d1a5a29cbca     0      0     0     0      0      0      0
#> 7a310a8f8e4bcbfb717263bd5376994d     0      0     0     0      0      0      0
#> 57136cbdae8fe066c2bb4a661003719b     0      0     0     0      0      0      0
#> 69fe46cab1cc9407470a2397bfc84103     0      0     0     0      0      0      0
#> b88efba97658a6936bb053985d026ace     0      0     0     0      0      0      0
#> 4d691fbbe2dbfb2815649fb72c757977     0      0     0     0      0      0      0
#> 49192cc8a222d97b965b839ff5d0f333     0      0     0     0      0      0      0
#> 859873fcd088c294a92a13340744eb6f     0      0     0     0      0      0      0
#> 0a8c44ca4d744e313e6dc5b239fca562     0      0     0     0      0     18      0
#> 09f94345824c82e8fe0d824cc4478fb5     0      0     0     0      0      0      0
#> 9dfbdb741fb0dc25d2206b8680bbf579     0      0     0     0      0      0      0
#> dc2bea463c720874cccd189ccdcf0e3e     0      0     0     0      0      0      0
#> 6273eb52991d528cf5281215f79b40a3     0      0     0     0      0      0     22
#> dc4bcbe74986e44ec3a040b866dfda9e     0      0     0     0      0      0      0
#> 96ec4781ee61f4083792f252f1b6cd3a     0      0     0     0      0      0      0
#> 9657d52acca51e701e7a3b7cedcfcc12     0      0     0     0      0      0      0
#> f65f183293d4710041057ab98f2c94ed     0      0     0     0      0      0      0
#> 8f6b3783b31a2dd8640de3df92dbb9c4     0      0     0     0      2      0      0
#> 0d90eb4842c9ad9126ed33c538c768ac     0      0     0     0      0      0      0
#> 8d3ee46aaaf7728594014be61c649347    19      0     0    37      0      0      0
#> e011a95fda1f1f2bafc4e5b83b52d989     0      0     0     0      0      0      0
#> ef8ea65ba4c2e2da3ba9556a87d2385d     0      0     0     0      0      0      0
#> dbbb92eb3a8a3feaf96092d46f2b6a45     0      0     0     0      0      0      0
#> 06cb03cddd8a12dbfaf27e8984b395c4     0    104    77   130     50    118     49
#> cdf134bb501c54d997945dabd44e35ff     0      0     0     0      0      0      0
#> f3c4561cc01a3b45ee3f134d709c6b94     0      0     0   114      0      0      0
#> b121442fd1eefe78ce4f4602aac1842e     0      0     0     0      0      0      0
#> ba8707aa88a2a62939d9338189b2a733     0      0     0     0      0      0      0
#> ecf7289628d49f9244469337ee1622cb     0     26    34    45     51     31      0
#> 9f06b400b93f859b11874b0bd7f09295     0      0     0     0      0      0      0
#> f775be5cf7c0dd2177006a44913503d5     0      0     0     0      0      0      0
#> e87ceb10f2becb2a6c4cf0692b1923ee     0     38     0    29     25     26      0
#> 036574bb8e05f3b9a11ab7a4090d6950     0      0     0     0      0      0      0
#> 0267ba9db97e00ae01fc2cccbb6263e0     0      0     0     0      0      0      0
#> be4f701c1dbfbda16c5150af29430141     0     81    61     0      0      0      0
#> 655a66e6566798ade91f73b2942c0b63     0      0     0     0      0      0      0
#> cd79f2696509d0ad906f5addba264079     0      0     0     0      0      0      0
#> 07ea3284bda2bf7daf1251ad2744e4c4     0     36     0     0      0      0      0
#> e3870da9de4a2b18a7ae9e9cae208e3e     0      0     0     0      0      0      0
#> e00bd59c1ea0ca9b1cc5f976801d4e1b     0      0     0     0      0      0      0
#> b39b381d0f533c9923c3b898ffab5614     0      0     0     0      0      0      0
#> e2959f0552fa6df5707d9ec8b57816d6     0      0     0     0      0      0      0
#> 53157acd1642f973842d01e1cd287421     0      0     0     0      0      0      0
#> 5991790d1c0aa8b73b2d7f906449a1f9     0      0     0     0      0      0      0
#> adee3552c211d9e33aaa28e82532bdea     0      0     0     0      0      0      0
#> 0089040d041888e3ae4aa4e0010a0e36     0      0     0     0      0      0      0
#> 8e876e13d051992c4c8616cf2db7d73a     0      0     0     0      0      0      0
#> 179abbc5c1f1549ffd0d7ae6963e3783     0      0     0     0      0      0      0
#> 340d095f9a9fb445582b60e1282d985d     0      0     0     0      0      0      0
#> 3d7e5432d8eaa3fa56ab7138605ea43d     0      0     0     0      0      0      0
#> d4b233740db907d69035915e8657fa8d     0      0     0     0      0      0      0
#> b22119df3b1de76745abaac585b646f8     0      0     0     0      0      0      0
#> 64761a24fc66e7c93145217dd153b1ee     0      0     0     0      0      0      0
#> eb6d7f9e0601207fdcfc4a88947d513c     0      0     0    76      0      0      0
#> 0d20271e3e2dbc68e1de868599a9a1ae     0      0     0     0      0      0      0
#> 369ea22772cb774950c830d0af7e8235     0      0     0     0      0      0      0
#> cae63066cfe088a49fcac50f119c9e2d     0      0     0     0      0      0      0
#> d9c57ec20ed90e2cebf6cacc4f9a5511     0      0     0     0      0      0      0
#> f5828cc3078b4c219a44c67e96e5d455     0      0     0     0      0      0      0
#> 73847ac2778d4cba27fe9758ca129d6c     0      0     0     0      0      0      2
#> d78003f13712da48d0bb36fd45116b91     0      0     0    77      0      0      0
#> eb40853d986ad72c82658d939f291ff2     0      0     0     0      0      0      0
#> 5d675a3518222fc99c8cba34feaeac81     0      0     0   111     52     70      0
#> 3ba390ae11c4e76375af082eefe7e41f     0      0     0     0      0      0      0
#> 3aa66daa4e560cf62c1833697bcd3837     0      0     0     0      0      0      0
#> 6f4388510be5ea70356f2257d8c57ddf     0      0     0     0      0      0      0
#> b33390d034148edf07d8a742778ddcbd     0      0     0     0      0      0      0
#> 09225ec7044394bbd3b6811a6799e481     0      0     0     0      0      0      0
#> fe6fbf9b86e0f487911a856009fa4a1e     0      0     0     0      0      0      0
#> 0aee6329e0b5b81449980151bf40205a     0      0     0     0      0      0      0
#> 0ef66a78236f90916afc2305f3ff6c80     0      0     0     0      0      0      0
#> fb1eb41359005c2e3e308fb668f07220     0      0     0     0      0      0      0
#> cb59c57800227006aead644df42b8986     0      0     0     0      0      0      0
#> 138da9dc244c3292536cacb8df79467b     0      0     0     0      0      0      0
#> 72f92361ec0696cb8df985b93d277a25     0     31    22    21     18     35     24
#> 4d5ab121763879fd59e8f00d35473540     0      0     0     0      0      0      0
#> 20032b2b586fee22805eccebe3059170    49      0     0     0      0      0      0
#> 710324ccc8cafa6fc9b38baf8299d8aa     0      0     0     0      0      0      0
#> 523e6b60340f4a215877b1c336822e50     0      0     0     0      0      0      0
#> 4501f98c83484fd902be11383c8d032b     0      0     0     0      0      0      0
#> d4cbf23757aad57c91b058cbb05566e4     0      0     0     0      0      0      0
#> de7964d5c22112693ed5796bb9f04b0b     0      0     0     0      0      0      0
#> f4a283be42989b920b83ea73ec72b25b     0      0     0     0      0      0      0
#> 2f3111649d0850ab86799f8bbadc28f8     0      0     0     0      0      0      0
#> b8acdd45cf33ab75fe4c0851a1787dc5     0      0     0     0      0      0      0
#> 5bb7fec39744a6d4e4d61596dabc6886     0      0     0     0      0      0      0
#> e8eddfb37c4215c86813163a27f202a5     0      0     0     0      0      0      0
#> 826e188c911c0846fd927300127c80c0     0     42     0     0      0      0      0
#> ff1fe599c44074ee6bdf919a119df93c     0      0     0     0      0      0      0
#> aff00055d3876f624d67fb80bdb59934     0      0     0     0      0      0      0
#> 1beded194f3ae4b57715b45e50adef46     0      0     0     0      0      0      0
#> 908d91feba7c1f1d8cc8a7145173d942     0      0     0     0      0      0      0
#> e7095743695f6d4337d01ceb9e29e224     0      0     0     0      0      0      0
#> 4cd420962855afdf0bc69bb268b846c7    52     58     0     0     25     78     36
#> 430ba5bea0fc08f65b69497bfbeb5c4b     0      0    12     0      0      0      0
#> e12cb2f0d2ce4401413176007dcd4158     0      0     0     0      0      0      0
#> ec615213f6c7febb65c9f2e45e438111     0      0     0     0      0      0      0
#> 412dd52cd54d0244ba8fc42ad140acda     0      0     0     0      0      0      0
#> 89192d9f779341d313f52aad0dc6a06f     0      0     0     0      0      0      0
#> ed543f92c61a06fc4e703655ce910a2f     0      0     0     0      0      0      0
#> fea1a70b0ef0f72e3bf27d8b295331a8     0      0    41     0      0      0      0
#> 6afe805b36f5aa29e7f2096aac42ace1     0      0     0     0      0      0      0
#> d9cc5a46831f24ec9192ae82dbc38b0a     0      0     0     0      0      0      0
#> 0d3eb38eb502c4ad13bd5a7060295411     0      0     0     0      0      0      0
#> 69a0ff5538a5490e291c75ce70325066     0      0     0     0      0      0      0
#> f91d545c815effa3c5afa83a8ad3397d     0      0     0     0      0      0      0
#> e6b032dea8e13c34af958fc48af09b58     0      0     0     0      0      0      0
#> 84a1b89045ab6285328f366224f3ec87     0      0     0     0      5      0      0
#> 4a6f0076e4ac71474b3c364294ad3b33     0      0     0     0      0      0      0
#> 7bdd39e4e3ccc7000862b487510e67d6     0      0     0     0      0      0      0
#> 998aedf12cdf054f5e504a2e21d20498     0      0     0     0     26      0      0
#> 567c2640108006c8eca081ad4180e346     0      0     0     0      0      0      0
#> e50a6625b9176b6856764b213013150f     0      0     0     0      0      0    106
#> 8d5a5ab7c784b0325a4a86a16fbba80e     0      0     0    22      0      0      0
#> f58ae28a309fa3d5898d2eb46d9eb9e0     0      0     0     0      0      0      0
#> 7ac65a8fd21111b68e7c2d6205112a83     0      0     0     0      0      0      0
#> d7add64c8700f5a8cadf40c8c7f5ff67     0      0     0     0      0      0      0
#> 6bad9d981d3cf0f90262547180dc6371     0      0     0     0      0      0      0
#> 149c0e63c83113e1246c0e041de67d9c     0      0     0     0      0      0      0
#> 6f196fae54ea90b5c11627c503c98940     0      0     0     0      0      0      0
#> b67c995d842e5b4b92fa0ce052603ef1     0      0     0     0      0      0      0
#> f44ffae67c7297e2c696b75fbd4e0f4c     0      0     0     0      0      0      0
#> 33997912cd952feaec18fbd11b31acac     0      0     0     0      0      0     22
#> 0a035423c4bf71a17566549f9b69448d     0      0     0     0      0      0      0
#> f65a6f1f1e5a4bc5c733ba364c0096d3     0      0     0     0      0      0      0
#> f15389a9f245465c33db3de72aeb17ba     0      0     0     0      0      0      0
#> cd2b676d3c0785a60359c0d416701ee4     0      0     0     0      0      0      0
#> 70ba2a378ebcd9a8f57b024a561cc6b1     0      0     0     0      0      0      0
#> 0273707384afdd41e9d022b539094211     0      0     0     0      0      0      0
#> f9ac98e876ecdb504202388d7ccdda5f     0      0     0     0      0      0      0
#> bdd9b87b782b8f935ecac62a9824e6a1     0      0     0     0      0      0      0
#> 1f12b97495668f1fc378a349546d1806     0      0     0     0      0      0      0
#> 83f1b4d776f6e08a6ba6f86634a1e5a9     0      0     0     0      0      0      0
#> 49f614420449f6c7ecc8c9f4e137ad05     0      0     0     0      0      0      0
#> aba5aa00a1b262f27517eca64ee3333f     0      0     0     0      0      0      0
#> 8971980484eeabeef00145b93b699e30     0      0     0     0      0      0      0
#> c1347ae54b0fc1dd385a1ccf3b8061ce     0      0     0     0      0      0      0
#> e9fe19783ebcc7e43a3a07ca8225db62     0      0     0    14      0      0      0
#> 7c523a5f0a7c736efd82afef84137901     0     16     0    44      0      0      0
#> e13c07ec336c89db5183cc6e9e57ee6e     0      0     0     0      0      0      0
#> a017ce431ccadea609b501d46e64e55b     0      0     0     0      0      0      0
#> fadf38d03d4b944f5782a759591ec3ae     0      0     0     0      0      0      0
#> 41bfda1c7f32205246e73a1e3e3ec4f4     0     19     0     0      0      0      0
#> 9e4f3e697b70758a7a996f3dd4533650     0      0     0     0      0      0      0
#> 39c4f6caa7510ed108f3ebcf72f79edf     0      0     0     0      0      0      0
#> 4aa083fe7c1b0639dbc48ec649e68a31     0      0     0     0      0      0      0
#> fd2ac2e14b72b0e22a052e3ab4580471     0      0     0     0      0      0      6
#> 4ee10294799bcfeab0c8ee567d5d9608     0      0     0     0      0      0      0
#> 46869f12c31de7a05781e930c25a1ccd     0      0     8     0      0      0      0
#> 5e2548d6a970645d0ca943a72d207489     0      0     0     0      0      0      0
#> 0538e1beed69f098c3f602130ee11a1c     0      0     0     0      0      0      0
#> 991c563278c4e8c91ffebce3560fe51b     0      0     0     0      0      0      0
#> ce05f58164e61517198f7ab0e3304cb5     0     14     0    12      0      0      7
#> fdbd0fed21bae91c05b2da268a025a89     0      0     0     0      0      0      0
#> 7412d11595d77f26792071cd28a26760     0      0     0     0      0      0      0
#> 4b087e13957eb38128f63e1c99df23c3     0      0     0     0      0      0      0
#> 310a6f0e6841db2ce75005cac942a972     0      0     0     0      0      0      0
#> 3df0fada803ac582802b47d4f4efe0c1     0      0     0     0      0      0      0
#> 3cd4e34c02ffdd1446df65bb54278c57     0      0     0     0      0      4      0
#> b14a1dea93554257791f16c2e94906b0     0      0     0     0      0      0      2
#> 902a4435fc3af101e8a9b5f51550fa4d     0      0     0     0      0      0      0
#> 1995185b31349529f58964c2e5a3becd     0      0     0     0      0      0      0
#> 69b247c565afdeacee0240482b1ce536     0      0     0     0      0      0      0
#> 95ff1b9528acca018e9c4697a8625504     0     34     0     0      0      0      0
#> 132d26d3ce6ac2ab25ccb4dd171879d7     0      0     0     0      0      0      0
#> bf1a6e4f94ee8ba60cf79089139641b6     0      0     0     0      0      0      0
#> 704189aa09f64b992b3a02ee001dc706     0      0     0     0      0      0      0
#> dea5ea27839ee67a31d18ba0bf0f7f65     0      0     0     0      0      0     10
#> 3840a4fff0250232f8570cae7c3514ba     0      0     0     0      0      0      0
#> 31ecabc032997aaf41dd0d8def3ec93f     0      0     0     0      0      6      0
#> 024a0e0348ed9abe9ae1c48b1b7ba09e     0      0     0     0      0      0      0
#> fc1d9940419420113ce3fbfacc8d703a     0      0     0     0      0     65     42
#> 0770fbf28cd307a08bf557af7ffecfd3     0      0     0    19      0      0      0
#> 883cbe8e4d1d98f468355e8ea974d222     0      0     0     0      0      0      0
#> cb6d6561383fb68f43f92a88a607749d     0      0     0     0      0      0      0
#> e07e0f1a9fbc696a9fe4f9b7739ae45b     0      0    11     0      0      0      0
#> 73d8d8f999acfa13cba2277f7526897e     0      0     0     0      0      0      0
#> 5218274c354d41965abbf4b0751d6f3a     0      0     0     0      2      0      0
#> 5c31570f706752fa4f7851aca5b5c291     0      0     0     0      0      0      0
#> 6a5e71c8b86a52ab1dd5e61164103940     0     31    21    22      9     30     20
#> 25641b0b803c7d8493ebc02cb104bf3d     0      0     0     0      0      0      0
#> 8844d70dd1fdb96f2a7faa65bccf78c6     0      0     0     0      0      0      0
#> 371b855afa8bed030e86ec78241b01a3     0      4     0     0      0      0      0
#> d0c180a47378d414a1d088ac14db236c     0      0     0     0      0      0      0
#> a072a588b3940fab3b227e4901c40f7f     0      0     0     0      0      0      0
#> 977a3523051e6b23da6221dbde795b34     0      0     0     0      0      0      0
#> 1c3adc87d2e93b08e576953508eb680d     0      0     0     0      0      0      0
#> 74eb0ddf379b00c5032b18fb9fbf1f32     0     25    20    15      0     27     16
#> c335b34a15db7d532f511e6c952abe62     0      0     6     4      0      0      5
#> 0c78285397fb853a0e39069213510c54     0      0     0     0      0      0      0
#> 779dec6de9c6a8b7df28fc8ec3c75bee     0      0     0     0      0      0      0
#> 5526624b473bb1f0541dc6fb95c98678     0      0     0     0      0      0      0
#> 649e7fa7b9ce5c5472d34675ed643639     0      6     0     0      0     10      0
#> a9c5addeea542953d1bae4d30d599b02     0     12     0     0      7      0      0
#> 68291fb3b2558d438da4810961c8d53a     0      0     0     0      0      0      0
#> 9c5d92bfefc7a190b93d5b9a6a77697b     0     16     0    10     13     13     23
#> 3ddab42d701754625cab0c236c5b6e6d     0      5     0     0      0      4      0
#> 287e840f9ecaef80d7a4762c8d6bf401     0      0     0     0      8      0      0
#> 885f411fa54c5576148678f665bb053f     0      0    12     0      0      0      0
#> b5709a57df9cbc8a741bb78a36f29c06     0     89    40    95     54     72     82
#> 26698cc60893017ef0b1ab173688d3de     0     34     0    15      0     26      0
#> db64f0028b6518839ab7129c4b7ef72d     0    147    98   115     68    182    174
#> f46fb43262da299ec7f0084298cafa50     0      0     0     0      0      0      0
#> c45e7087de669be270881fc8ace0fad8     0     45    69    44     42     33     75
#> 854a20d5019388e81b1e78419eb017fd     0     30     0     0      0      0      0
#> 9a2325ea015a8c5944047e3349cc1565     8     21     0    30     24     16     34
#> c5305ae3c6adc90814d70ab48f4f071a     0      0     0     0      0      0      0
#> e230123ec3ce8654942ee1a000a80010     0     48    36    62     33     41     74
#> 14e3c60cb885509667670151bc4a1df0     0      7     0     0      0      0      0
#> 7fd7b7bece92bbbfcceb2dcf1e835a8e     0      0     0     0      0      0      0
#> c098d1fb769dd975c2332dfb70da489e     0      0     0     0      0      0      0
#> 95dcdf1eb06b252e67dcdbdaffc505e8     0      0     0     0      0      0      0
#> c5399d0f9eda8f0815fd28a257952496     0      0     0     0      0      0      0
#> c7388cc314e992819b8f97ea35280982     0      0     0     0      0      0      0
#> ae8facc45d372f7108bd9499f74d1bfc     0      0     0     0      0      0      0
#> 16caf6f5bc653fb6ba668eb821076a21     0      0     0     0      0      0      0
#> e7e283b72ee15bcf0877c5eeb6137ebc     0      0    11     0      0      0      0
#> 1b3a779b384c73b33207dccd14b14fe6     0      0     0     0      0      0      0
#> c8d1afdf7b395ffc22193c5dc2fd5dc3     0      0     3     8      4      0      0
#> 0e3ffc2a68378086e2480831acd1a025     0      0     0     0      0      0      0
#> 0ddf7339b608a340252b33f63fcc0019     0     18     0    19      0     14     12
#> 5fd4b2f7dee3fb0bde6f856ae266dbcb     0      0     0     0      0      0      0
#> e2f4786bc66410c07bd8b96b5d0e2d92     0      0     0     0      0      0      0
#> c30aae8f5216c9536492d0beae74950f     0     63    53   119     43      0    106
#> 9087ecee0806a80861788a17dd7e6809     0      0    11     0      0      0      0
#> c29bc08418318360592b2b2717376a0a     0      0     0    20      0      0      0
#> f791c7db1a78d35e84803b5e03c01965     0      0     0     0      0      0      0
#> 44b44fc4055a58c01f6d1e9de5e0d1d9    57     89     0   137      0      0      0
#> 342df9401103d274851d3bd80704508d     0      0     0     0      0      0      0
#> 0590e6672bb2c55c0ff54999a2b6dcfe     0      0     0     0      0      0      0
#> d4e9147608bf883fdbbf9d23b17028d1     0      0     0     0      0      0      0
#> db24794eb8667aed224401386223f671     0     35    49    12      0     38     92
#> fe41e6fd9770030a21418b288930710f     0      0     0     0      0      0      0
#> 699f18f3b155bff97de941f2cfbbbadc     0      0     0     0      0      0     86
#> 13870cf319997750e3b96bbb0b2f1aea     0      0     0     0      0      0      0
#> 7ca4fb1b1740406241acb862445cf311     0      0    31     0     16      0     44
#> 9c0569fb95451d009fd3754303eaed12     0      0     0     0      0      0      0
#> 1ac26a3c12288d8fbadd6aa61366b1b7     0      0     0     0      0      0      0
#> 89f3ba99552605501b8acb816d288afe     0     37     0    33      0      0      0
#> b4051b5227746f04c1feb9d8979299ff     0      0     0     0      0      0      0
#> 087d358aaf77f81aaf50a623828865bb     0      0     0     0      0      0      0
#> e451c4b3ff1c3dabeee296529b3a7a07     0      4     0     0      0      0      0
#> 012871784f77c9dca7d446ddef2048d1     0      0     0     0      0      4      0
#> 95c6ce2788783bf2f83343deba23b8bd     0      0     0     0      0      0      0
#> 213227f244d97553a9a14c4e3fd04c98     0      0     0     0      0      0      2
#> c0e5e55c5cfbd576e7c085f1e4b69638     0      0     0     0      0      0      0
#> 51e98451d7e3f46127481ec5490f1c91     0      0     0    30      0      0      0
#> c75285a6c61d8fce20d823769c25c12a     0      0     0     0      0      0      0
#> 177966f7619eb24fdbce2266f5bed39a     0      0     0     0      0      0      0
#> b89347f0534679f890628364d2c7b4fe     0     36    37     0      0      0     17
#> 8dff1cb12c0beb5bc7b2546480a6a6c7     0      0     0     0      0     23      0
#> 84f2777cfcb0fbed6281917fcc765022     0      0     0     0      9     15     19
#> a5e4d106a2ab947984fdef187996a0b9     0      0     0     0      0      0      0
#> b6cbf6077b64acbb718971a1ffc97a27     0      0     0     0      0      0      0
#> 27e29495fda424728fd6a80b83a16772     0     22    21     0      0     29      0
#> 63c44ab3725dd0b803d43a29db9fd7fd     0      0     0     0      0      0      0
#> 9c931ab5bf18ef9b10996d59d292e6f5     0     26    15    18     16     45     21
#> 12fd8efb82c3d4909543f6db9853134e     0      4     0     0      0      8      0
#> a8a60c20f7ab7d0f2dadd9f8858c19be     0      0     0     0      0      0      0
#> d57e09eb1645f1abba6074d3bd7d4c52     0      0     0     0      0      0      0
#> 7e0f2b982ca2d35b6ac2cf251d78d3c4     0      0     0     6      0      0      0
#> 3e2e8e938a68a0318bcf7145ad5ed4a9     0      0     0     0     10      0      0
#> f2117b49bf57d9451bbef6de5b762483     0     15     0     0      0      0      0
#> be24dfe003725a0b0f28f4136106faf8     0      0     0     0      0      0      0
#> e63301d7d0af55d5c43db051dc51522f     0      0     0    11      0     22      0
#> a6e2e5501507562f2174f65e606b2d43    15      0     9    11     11      0      0
#> a06aa93cb50b9a39941a11c7242dde96     0      0     0     0      0      0      2
#> 6a0a8abb36aa3bfdee14edbe97e9d781     0      0     0     0      0      0      0
#> 4d3def5ffad0c392381591440bfc3668     0      0     0     0      0      0      0
#> ce437b81e7cf5ce986f7240c26464b0a     0      0     0     0      0      6      0
#> e57788f9d2fae76ed9ca566ad4d475a8     0      0     0     0      0      0      0
#> 8c525d0f63137bd068b75c57f5d04c09     0      0     0     0      0      0      0
#> 6d670aa40dafa152c2d59cf2be9bef74     0      0     0     0      0      0      0
#> e08f8e9e9ffb10a702f2ef617f70eecf     0      0     0     0      0      0      0
#> 7def1c807ca0767ccdfe5a6001dc51f1     0      0     0     0      0      0      0
#> 3b9a17275147adf1dc8c42eb00a06a36     0      0     0     0      0      0      0
#> 24be0b7e4c85ac53b09fddbca916ab96     0      0     0     0      0      0      0
#> f521a47d3bb4acfec6b261f880b4eea8     0     12     0     0      6      6     11
#> 91c24c9b6ffbedcc257d164176060fe0     0      0     0     0      0      0      0
#> 9b45f210647c5fe71b005eb9e25299a5     0      0     0     0      0      0      0
#> 974077d25fa6b32a7af22670b142a8f1     0      0     0     0      0      0      0
#> 203aa2535e54a11e4c48af4ab946686b     0      0     0     0      0      0      0
#> 9064afe569448558eb08e5bf5e458bce     0      0     0     0      0      0      0
#> fad432eeeff7f1d923e6d60aafe18dcd     0      3     0     0      0      0      0
#> a29db0e5f4a04c5c6937745ccccb167b     0     11     0    10      5     18      0
#> f179c10f8e86873c90ceb14b974a114e     0      0     0     0      0      0      0
#> a23e8588a2240d7af6f65519158c902c     0      0     0     0      0      0      0
#> 52034e0da6d5e7962369437956984e18     2      0     0     0      0      0      0
#> 522d20c27cf4237b60f36ffdca056c19     0      0     5     0      0      0      4
#> 957615aefb421e6ad4728e610b2bdfa6     0      0     0     0      0      0      0
#> 5ccc74e70b8a82a6eebea3385738ece3     0      0     0     0      0      0      0
#> e5924fd11a48419c6ec0e2a7cb4e43ee     0      0     0     0      0      0      0
#> eef7be4b2531cb7e909b58e91af5640d     0      0     0     0      0      0      0
#> 4416cbc2026260fa4e2ed0a8feff1443     0      0     0     0      0      0      0
#> 03f6091f5aacb69f06fc925acd1c9506     0      0     0     0      0      0      0
#> f655e5ddf5074d5e71bfa583926dec30     0     10     0     7      0      9      8
#> 9f733d39e588c7310499dd444cdd7812     0      2     0     0      0      0      0
#> 2b413de5e05a1385bab229f50879628f     0      0     0     0      0      0      0
#> 178b596f2fe5cb6cefc5513917988d84     0      3     0     0      0      0      0
#> aa94b014c726cf40983045bd606241c1     0      0     0     0      0      0      0
#> 55310b3aac4868c4e33b0a9492d7609a     0      0     0     0      2      0      0
#> c5d8d828df0e1baf5c2db05a6b8c4476     0      0     0     0      0      0      0
#> 6ceb26f60606233a7728d3457af84380     0      0     0     0      0      0      0
#> 1ff83868d1c9ecf52d9b776e0b3fd896     2      0     0     0      0      0      0
#> 3edfccfb9f81b340854e84ee5102a4c1     0      0     0     3      0      0      3
#> 963dc552c551872800f2fcbc5663e5fd     0      0     0     0      0      0      0
#> 3fda8c20e5da9fbc3fb7c5c864559019     0      9     0     0      0      0      0
#> 6fd88fa7e883bdead3c5ba89da853e8d     0      0     0     0      0      0      0
#> 35826fe826b7a9f202e8bef7b39d9a47     0      0     0    15      0      0      0
#> cc514f614ba146f489bf62b5a56074c5     0     24    41     0     30     12     57
#> 26260de3858081ce55e6f220b3ba25df     0     15     0     0      0      0     11
#> d6727a85782ffd2a069dccc23adaf521     0     64    78    47      0     49     97
#> eb1a74c2ca8831ef945775de0827ddda     0      0     0     0      0      0      0
#> 4ed3218615bbc1ed4da82d9806e4786b     0      0     0     0      0      0      0
#> 8a3367e36ed4e416d5a77a7946cf6fc5     0      0    12     5      0      6      0
#> 56f45eba0496ed29fb7171041a2a950d     0     13     6    18      8      0      6
#> 9d2782f19e60ec0cf72e1dc46b1d1766     0      0     0     0      0      0      0
#> 1559f8b64dc2bdc8cdc26ff8a5ffa4cb     0      0     0     0      0      0      0
#> 6d01f384a1cce140edea17dccd0f116d     0      7     0     0      0      0      0
#> 5950949b71082741302b712e97309b17     0      0     0     0      0      0      0
#> f6f238f95cc97995c9bf619b0c620b33     0      0     0     0     20      0      0
#> c7a112ffadf8009efd318625da44b18e     0      0     0     0      0      0      0
#> e79f4f30b0dacfcc50f3c1b12804ca3f     0      0     0     0      0     11      0
#> 297c825dc268f9a7b86b402edb6961b3     0      0     0     0      0      0      0
#> a1ed56ad4e86e0b4adc11baf2c4ab94b     0      0     0     0      0      0      0
#> 8afccb8f7f50c5603b471fb0c808c9e6     0      0     0     0      0      0      0
#> ce41c973bd4965a400989ec128890d55     0      7     0     0      0      0      0
#> 5d99ba689ad34e4393a905e8a725b72b     0      0     0     0      0      0      0
#> bfaa2d35e6dd13b9b32de811c6f92592     0      0     0     0      0      0      0
#> 0c544d5f3a49a010971a96d8c39e4905     0      0     0     0      0      0      0
#> bab3f08e1d75f42987012cae451639ff     6    175    72    77     85     69     89
#> fa97a81b8c33569e2aa1fc9a3137b435     0      0     0     7      0      0      8
#> a4d08886c3c6c62efd5678339a07828b     0      0     0     0      0      0     14
#> 11ed0972601847ec168473a6bc4eeba7     0      0     0     0      0      0      0
#> a6e4d88aa4b40402b1f94760d86b729c     0      0     0     0      0      0      0
#> 86deaff47d33946ede1b7c25c6bbe9f9     0      0     0     0      0      0      0
#> a1d3f1f8c151eaeb5cf12179861e0c87     0      0     0     0      0      0      0
#> 3980891c22ccbb387b65f9e437b2f94f     0      0     0     0      0      0      0
#> fc8ec321e0c2fc2a6de9abe98256b3f5     0     67     0     0     32      0      0
#> 84c11e6ffb29ebc2192f4b049b7dbddb     0      0     0     0      0      0      0
#> 82104585f0b2617c13eb69f0bbdf7d5a     0      0     0     0      0      0      0
#> 61cd4918428d3b36910b9cb2fc45946a     0      0     2     0      0      0      0
#> 0ab9a976935482bce5f8af4cc6685147     0      0     0     0      0      0      0
#> 87be3450216b45c53851fd173582a2d0     0      0     0     0      0      0      0
#> 93a51bb503fadf90121f518db980e158     0      0     0     0      0      0      0
#> b70cd843161996e075a853900e699d0f     0      0     0     0      0      0     11
#> 154b8be9eea158f93d27d314d8e0e2c8     0      0     0    11      0      0      0
#> 70691c2737cf7dc7c194de819f0e6cda     0     17     6     0      6     53      0
#> 2ae1d80ead92f5991f58da14d000acd6     0      0     0     0      0      0      0
#> 44fa299492152518fec34b9fa4554648     0      0     0     0      0      0      0
#> e15c8f9df87f363b0dbbf49de0bed249     0      0     0     0      0      0      0
#> 73689cad6eb4e86d9d0bb48c0da24463     0     17     0     0      0      0      0
#> 234e15e82d64f2e1d61802fa8ea1998a     0     40     0    39     36      0     47
#> 88c76ddb1b24acffe65d8128e2f48244     0      0     0     0      0     24      0
#> 5544139a16716cc6ef48b2c16c9613e9     0      0     0     0      0      0      0
#> 5b7fe68483e701059142cb531397b547     0      0     0     0      0      0      0
#> 0a0875a0acf4a28d1bc1c4d1d11aaa64     0      0     0     0      0      0      0
#> 6b47ec9d1c70f07d16834e9a4e54aa70     0      0     0     0      0      0      7
#> f8752b4824d239cfb19f3216cc34941b     0      0     0     0      0      0      0
#> ec8a0768cd2fb426d7e6bfc2c4f50e5b     0      0     0     0      0      0      0
#> 8a6b88f9515bf2dfac09e6a69045a4fe     0      0     0     0      0      0      0
#> f5151c67af19c8471743ef3e2914f6b1     0      0     0     0      0      0      0
#> c3842619fcbde52781c512c2561fca0b     0      0     0     0      0      0      0
#> b3b00faa9dee9fa6b3b33744f51da00e     0      0     0     0      0      0      0
#> c4f04e28c6267ac6c5b05d403e8f91c3     0      0     0     0      0      5      0
#> 7fc06f107a6ddf6486ef2d43fcf23626     0      0     0     0      0      0      0
#> 8d1dffcd377b338d2a0388637507b592     0      0     0     0      0      0      0
#> 15410d762889d51818986bfcdb5d8965     0      0     0     0      0      0      0
#> 48f68e44e258de63575b8bf6cf565d52     0      7     0     0      0      0      0
#> 4712befe3685661e0fc0968859c369d5     0      0     0     0      0      0      0
#> be9a3727a0422aea147100370d046fa8     0     41    24    56     21     34     20
#> 3c5942781761830f7b5406f551574424     0      0     0     0      2      0      0
#> 291d9f4d3d3c5cc3b0a1804770ed1e77     0      0     0     0      0      0      0
#> 784764c519a64dc9a147025a23222e93     0      0     0     5      0      0      0
#> b9e955ea5254dd68b54775aad2a67829     0      0     0     0      0      0      0
#> 8197f1594e8a17f3fc989fea64f54b8f     0      0     0     0      0      0      0
#> bf4a4f2383a16d2ecf27bf8e162b9ba4     0      0     0     0      0      0      0
#> 4e7d22d793d30bdedbd32a5b73b564b8     0      0    12     0      9     15     23
#> e88959b739c8cceb94b94b7887dd1df3     0     15     0     5      0      0      0
#> 5403445694339e2ae5ecac06d9f3ff47     0      0     0     0      0      0      0
#> ba3943fb075bbf63fb699638292826a1     0      0     0     0      0      0      0
#> 67b2d836db093cead5a8286142669613     0      0     0     0      0      0      0
#> 6c3a0ea164d3eb098dd2baac5519b38f     0      0     0     0      0      0      0
#> c643a776c36bcc00dfb5d790e1ea3dfe     0      0     0     0      0      0      0
#> d53cc81a9d57a64d194c876b76bd66e0     0      0     0     0      0      0      0
#> 7dc00d40971bb454b39fdd9b2dd9aa25     0      0     0     0      0      0      0
#> b4d0cec8543e940be1b5096370fdfb80     0      0     0     0      0      0      0
#> 6c54ceb94bc5e8fbd65224cffc60990a    15      0     0    11      7      0      6
#> b12ee9ca5a2320957213a114788cf2b7     0      0     0     0      0      0      0
#> 5e040c2397f959e7f0dfc6ec2855b812     0     12     0     0      0     15      9
#> 3105a594bba2604ae5ed19f5e0495d2e     0      0     0     0      5      0      0
#> d75c040bdaa1cca96cd0655c9af7ab28     0      0     0     0      0      0      0
#> 9b909edaebabc90b9c1a703626e6bd0c     0      3     0     0      0      0      0
#> a6593da81035cf121dbe46346a4ea960     0      0     0     0      0      0      0
#> 60e4c88400a750a8122b97476313a67e     0      0     0     0      0      0      0
#> 7854a78e3b881cee1c63ea0a50ba29ce     0      0     0     0      0      0      0
#> 1796734b56d45bdc997ccc09ab8731d1     0      0     0    13      0      0      0
#> 8243cd669ae7327f4326be78429ab05a     0      0     0     0      4      0     27
#> ebb93dc32c63402289114859d476cc5b     0      0     0     0      0      0      0
#> 8a01af5a531e851e2c50fc01131cda33     0     17     0    13      9     11     11
#> 03d2cc3f7b2bc352752d39e3738d0af9     0      0     0     0      0      0      0
#> ec378af6e32960b11b5f7391da4ae4fb     0      0     0     0      0      0      0
#> b407f2de868c064f65eb5ad8f82a8d0d     0      0     0     0      0      0      0
#> 73a95b6ba34fcd260989d848fe15d377     0      0     0     0      0      0      0
#> 3abb68202d15781590a2d1f9c5b37e15     0      0     0     0      0      0      0
#> 77923ef8e9ed01c9b6a8ac85eb028163     0     15     0     0      0      6      0
#> 9981b8fd20f400e052527f84f5db1b66     0      0     0     0      0      0      0
#> 597bdec49adf7675912f9fba4e262e70     0      0     0     0      0      0      0
#> 0fd3126118b13fb187bc55476a67c95a     0      0     0     0      0      0      0
#> 2a0c511ecc2513c42fe624e8e99d07da     0      3     0     0      0      0      0
#> d2912fea2dcf2e5b97e56c6a7d24295d     0      0     0     0      0      0      0
#> 2976af6aaf4261b6c274fcc3229a4d65     0      0     0     0      0      0      0
#> f097c45ce419da576aa111c98a4c7272     0      0     0     4      0      0      0
#> 3abc8260d34008f263e232ea25426c3f     0      0     0     0      0      0      0
#> d7e5b6432dc178aaeed9c81e0ec57afb     0      2     0     0      0      0      0
#> dfe71ce78ddfdbc7b453d17d628365f2     0     17     7     0     13     26      9
#> 289b18a86084273ce3b4038bf37e2840     0      0     0     0      0      0      0
#> 844a67032c52dc49019bf5dff248bf9c     0      0     0     0      0      0      0
#> 38cb3e3ce5ac5aed7dec50247ecfb267     0      0     0     0      0      0      0
#> 74e88e71e9929b1ddfcafccfab3b7bba     0      0     0     0      0      0      0
#> 66bc0737fc7fa1463ac7189ab7073998     0      0     0     0      0      0     32
#> b5a7aa588c9cae6ae8ca2f7d87ac53e8     0      0     0     0      0      3      0
#> db0540920457b21b170947aac9137771     0      0     0     0      0      0      0
#> ebbaf260510ab16ac6b8c2b9866d05b9    10     27     0     0      7      0      0
#> c1f461b67419d4a80a61db0867da9ed6     3      0     0     0      0      0      0
#> 0c6eba6452635eecfb865a1049299176     0      5     2     0      0      4      0
#> 7ecabb10f29a0c44eb53e9f9852fe22a     0      0     0     0      0      0      0
#> a33cb3037e7825920073b4dddaaa41e6     0      0     0     0      0      0      0
#> d625812c15c074c2e2b18099f97412c9     0      0     0     0      0      0      0
#> 17ef77f12e2a55717e910a5f6a625448     0      0     0     2      0      8      0
#> 7baf638fdd3f71166aafd22ecf410348     0      0     0     0      0      0      0
#> 9436cda211f60f722865cd23ef75d4ed     0      0     0     0      0      0      0
#> c8ab0b20e740eb6b35e68d58f69736e4     0     66    46    23      0     60     45
#> 5e43f2d8be71399fc8d6736754cb42ac     0      0     0     0      0      0      0
#> 19ad90ee1081a003406b02179ae63709     0      0     0     0      0      0      0
#> 0a25cdbd0660e68970c563c43a18f54a     0      0    23    31      0      0      0
#> 0aa6b210e902d8ab9d6f267407258ebb     0     16     0     0      0     19      0
#> a864c36d0d898cc5abd94d6d6f0bd5be     0     14     0     0      0      0      0
#> dea0dcb305ff0f27e3c36c7133cf933d     0      0     0     0      0      0      0
#> 463977baa69c7724fdc9407ad162714c     0      0     0     0      0      0      0
#> fcb72326c29342bcfb052953466eafc6     0     11    19    30      0     39     23
#> 5a554353f9ed728f0ab0631f29c8be2e     0      0     0     0      0      0      0
#> a9dd27b48e443415a743bf9902f515c4     0      0     0     0      0      0      0
#> 77a40bcf979c0099c52f19a3f1fe5217     0      0     0     0      0      0      0
#> b4007ac525a9b409d70426e24421cff0     0      0     0     0      0      0      0
#> 406ba18a175c4c79efb25f48a6beafe4     0      0     0     0      0      0      0
#> fe1c166755bc722cefd2e2a3b2b70672     0      0     0     0      0      0      0
#> 8afd274849ac07c0ca66f8f5f143df3c     0     15     0    61      0      0      0
#> 197301ed9d8728b5c061b99a17973161     0     50     0     0     55      0     29
#> 14779f9dc3505598ad3d799e4ee676c4     0      0     0     0      0      0      0
#> 4085c246640eee043d523ce2d799023a    21      0     0    65      0      3      0
#> bba09a403e67c1381de01c86511f1057     0      6     0     0      7      0      0
#> 392c59fee4e6b91d7654f2d3e3c65be8     0     19     0     0      0      0      0
#> d2b05ae73ad5ed296d7add9e2bb2cce4     0      0     0     0      0      9      0
#> 8306a9473dabfd5918ea869dbfcd1291     0     11    12    10     12      0      0
#> 06064c372bb624503c2905d13f7691f5    37     54    63   324     49     99     31
#> 3265692341e42dbd5cf4ed1938b906e3     0      0     0     0      0      0      0
#> da96a0797b9881a123795838b9f01140     0     11     0     0      0      0      0
#> 0bdddf10e577657becec88e42001fd16     0      0     0     0      0      0      0
#> 168999a950ebe17e6d2ea1f90e557d91     0      0     0     0      0      0      0
#> a266ffe9d38f626d106fa615da4618ef    34      0     0     0      0      0      0
#> 7706405d5bce7715f01d8b54188def42     0      7     0     0      0     14      0
#> 4b0cfc8899b466bd9992f7b22a5a06be     0      0     0     0      0      9      0
#> b81b2b4a7ca9eb2c23b170648dfaf22d     0      9     0    11      0      0     16
#> 75cf40872201f150cb5b191676f64b78     0     51     0    27      0     34     43
#> c8e5794c58de9326084a62d3282d257e     0      0     0     0      0      0      0
#> b0706a571f169c0bac7b4864052fcaa7     0      0     0     0      0      0      0
#> 0388f969e3527b2938acc676d6757bdb     0      0     0     0      0      0      0
#> 2dfe751390aedfbd9f62dbb77d1f10b4     0      0     0     0      0      0      0
#> 2548dc9697e1f9afe04a9012f40e2f0a     0      0     0     0      0      0      0
#> 2f019615d5141a4ab00b1aab98882a33     0      7     0     0      0     23      0
#> 629e4f0cf1826904ff74907b2507f58b     0      0     0     0      0      0      0
#> 774861f348cec8451102730c21a790b0     0      0     0     0      0      0      0
#> 11c0d6de28ab026a9ae17abe4efdd6bf     0      0     0     0      0      8      0
#> c11381dda0841c0d08cbfcef3a760aa3     0      0     0     0      0      0      0
#> 3c9619d4f9db4ac6ad7151a0c9d86ad3     0      0     0     0      0     17      0
#> 7e2c72204a8c45b1b0ac39e0e6047098     0      0     0     0      0      0      0
#> 8254e0f9f07709b62dd8483f3b25fad0    43     21     0    35     21     35     15
#> e78135f2886d427e678d12c446bf070e     0      0     0     0      0      0      0
#> abec938a8d53e7471f3f476d50694be9     0     17     0     0      0      0      0
#> 4963567d2979bb761291aa98e0b8595b     0      0     0     0     13      0      0
#> d4de9d1453102564f35262faa3952266     0      0     0     0      0      0      0
#> efa4a78c4b993e9d115841c8a02c795e     0      0     0     0      0     17      0
#> fcacd7424bf3543322ea3e9842d14e86     0      0     0     0      0      0      0
#> f941e15ea7f7560191d57abc4d2bedd8     0     22     0    13      8      0      0
#> ac294f5242de3c9b26e75502c48f6383     0      0     0     0      0      0      0
#> 7b6a8e3f8bfa1bcd1b5fc55db048d005     0      0     0     0      0      0      0
#> c9e776ab599fe75db286675b698c9067    16    136    17    21     68      0     37
#> ebe5cb945f6d6d849c07424591f80753     0      0     0     0      0      0      0
#> f116387f77802409103a559e4205df06     0      0     0    23      0     24      0
#> 93f9939b261cddc9a57071b42771c1d4     0      0     0     0      0      0      0
#> 971a5f48b3418449ca7cf3e190a5e831     0      0    19     0     19      0     21
#> 730de992ca06a95e7193a93146a97cf2     0      0     0     0     24     44     37
#> 8fa8b4ac69784c26901e036d13ab2332     0      0     0     0      0      0      0
#> c363c9efaa85023b488b732b97e47270     0      0     0     0      0     36      0
#> 2e2181cc893ca1e1e677e888bf8af441     0      0     0     0      0      0      0
#> 86fa662da6869634e7f90f85e9ee4933     0      0     0     0      0      0      0
#> 9a02177f2c235a6eb62984d5aab151ff    26      0     0    50      0      0      0
#> d6854be498e22b65da3c00572ee178f0     0     37     0     0      0      0      0
#> 015f364f150b230450da0a123e7bb8ab     0      0     0    20      0      0      0
#> af2c1b22758cd6653e8e841497d096e7     0      0     0     0      0      0      0
#> e7f2fab8c0e7c365768a43ef71406605    16      0    10    28      0      0      0
#> fb57cf192bce1876471f2f6181bf4b55    10      0     0    18      0      0      0
#> 6c740ce69c9f45bea684dc9e36054006    15      0     0    33      0      0      0
#> 3b27c7bf279535d8e39664d8dc1d0daf     0      0     0     0      0      0      0
#> ceea3c29ad012c2918fb6aa21a44e8cb     0      0     0     0      0      0      0
#> fa99ecdb27c1cc44fd866ee1da250c4c     0      0     0     9      0      0      0
#> 387cc83573db0e444e81c89a8a4f8f9c     0     29    22    43     27     31     30
#> c7a8d1af5f7f472671f700775de3e0a2     0      0     0     0      0     12      0
#> 98a6a84550e07ea252cfd8a0f51042d4     0     48     0     0      0      0      0
#> 08b809a68126c38db786cafe9710b013     0      0     0     0      0      0      0
#> 8db4cc0a5d2f03fce8efd2b8374ff359     0      0     0     0      0      0      0
#> 3e7c7d38423c93b43a1d9865cf0bc4a6     0     33     8    17      0     27      0
#> 6f02050760c030ed8cd1ed836697fa6b     0      0     0    16      0      0     26
#> 8fc1865cfcae32f72c1b1570299e2ef3     0     19     0     0      6      0      0
#> 38cc138fac95b006929361093d840a00     0      0     0    23     10      0      0
#> f5f8eb63e22362509e2c7c78271305f2     0      0     0    52      0      0      0
#> 0bb238aaee3ccb7dc8dfa86dbc1ba9f7     0      0     0     0      0      0      0
#> ca4546caf07edf5200cf79cea04aebcd     0      0     0     2      0      0      0
#> 893d12123eb7f623489a4f025a15ae7f     0      0     0     0      0      0      0
#> 192c5a817f48368339e519e0a0a907da     0      0     0     0      0      0      0
#> 3467185ca92b7bec2523214605adabf0     0      0     0     0      0      0      0
#> a2aae80d8fb922fe1074d243b59ce237     0     41     0    15     11     20      0
#> f11233f16d4177d38ab1c0bfc649e056     0      0     0     0      0      5      0
#> b219ca7793ab9f5e441f2a68bc76152c     0      0     0    28      0      0      0
#> e6e1848b51e705e55b3805703f088311     0      0     0     0      0      0      0
#> bd074b9bf365d11468bd76268d63b15b     0     31     0    23     15     31     15
#> f381b492d6651c52fd09ce9351e7eb2d     0      0     0     0      0      0      0
#> d9b31e1c30facc24fa778b1d65f5c457     0     30    17   289     18     17    272
#> 14a67c8c8cc31d203eac8109443fab5f     0      0     0     0      0      0      0
#> 026ad0ea53fb202ebb770acf65705d6e     9      0     0     0      0     14      0
#> 19ba9fea04969f975d6dc8d3da20d5e6     0      0     0     0      0      8      0
#> 0fc4426ed7dccc48ce27cfff8083433d   228     13     0   343     27     30     63
#> 45ad90e179842705fa984962f131fd61     0      0     0     0      0      0      0
#> d284afab226d22fd0596221d17db23d3    10      0     0    26      0      0      0
#> 38194655ba17dbfd7fcc1ed42a57c169     0     21    22    11     19     62     33
#> e30be5268bbebe597424bdf171cc1607     0      0     0     0      0      0      0
#> cb7ccdb4956134b452a3cb78d1751c57     0      0     0     0      0      0      0
#> b585436677ac88dd4ac6239655f2a268     0      0     0     0      0      0      0
#> f5d770eff5c9b9b0c859db57d175eb34     0      0     0     0      0      0      0
#> 1b45ffbc88c788dcee2239e68901cf05     5      0     0     0      0      0      0
#> 3fd951cf7cc3a0740a1f20f8a9c6d2c2     5      0     0     4      0      0      0
#> 2027d3e90ccc4003f1954322e99db6a9     0      0     0     0      0      0      0
#> 6157cb92771830955413f808e83fc420     0      0     0    22      0      0      0
#> 6349506cac32056a214b04239bdc3a33    35     13    28   109     23     56     27
#> 2254a8d37cf537b1e35097d70b6ac91b    34     89    97   240     59    213     69
#> 38f3bc8ef8836574de0fcbeafd348b65     0      0     0     0      0      0      0
#> fec23c70ae0c70ecceea1487cd978ee5     0      0     0     0      0      0     10
#> 503d123a6c78c6fb280a6472bd61b51e     0      7     0     0      0      0      0
#> 9f33de8570c639e727465f1ee4d750f6     0      0     0     0      0      0      0
#> d68e1899fa65515135d9bd8d356b7e3d     0      0     0     0      0      0      0
#> e15b88ccac7a4bef64215e915f6e0d89     0      0     0     0      0     16     10
#> 767042b4b19d8bfaaa426dad7dbce378     0      0     0     0      0      0      0
#> 4a3c4abee502d6185abd750522566bf3     0      0     3     0      0      0      8
#> c254aa5c00b7f27d756f7bbe96b3d4ad     0      0     0     0      0      0      0
#> fe8d8af791eb6e5d4335bb7604c74e9d     0      0     0     0      0     14      0
#> d81800d7e31b69877ce1eab9858b8f3a     0      0     0     0      0      0      0
#> c55c7cb80c5700ff9aded76192543b5a     0      0     0     0      0      0      0
#> 91dee7b3b90bbf2e53b2472b9bbde0e9     0      0     0     0      0      0      0
#> 6b28df61cb5dee74b76b1a13aaa0027d     0      0     0     0      0      0      0
#> 9a830fa26205ad91a02b5bcaebb0221a     0      0     0     0      0      0      0
#> a54b8142d2a09d058b81c377e8e286b8     0      0     0     0      0      0      0
#> 69ef3a03a4d733dbffd99fdaa1d2e62a     0      0     0     0      0      0      0
#> 083cbb131dd305a88b7d338bc371c35f     0      0     0     9      0      0      0
#> ef1a92caffb4d7cae286fdc77b0a93a6     0      0     0     0      0      0      7
#> 2cddd5dca166486f50a0828b28001e1e     0      0     0     0      0      0      0
#> bb899d1dcfd89cb200a76cb8951904f8     0     22     0     0      6      0      0
#> 7817f6851e9df03a3d992cf915310dba     0      0     0     0      0      0      0
#> 54eaf6571de743a86ccb531bee8aad7e     0      0     0     0      0      0      0
#>                                  3.7RO 3.7RI 7.9RI 4.12RO 7.1RI 7.23RI 7.23RO
#> 374a5eb6496d14e8e9d6f3d8c7e34d9d     0     0     0      0     0      0      0
#> 4c47d5cf81df1ac6e2a0560029a81578     0     0     0      0     0      0      0
#> 8fb6741a6685fad0374f85ad3e16156c     0     0     0      0     0      0      0
#> 310f3de009c95de5c938b8e811dfc96f     0     0    43      0     0      0      0
#> 7d3b39e8c0f0fdcb81b0c898a98c38f9     0     0     0      0     0      0      0
#> fc211549300b0954dbb4a4bf57e8a605     0     0     0      0     0      0      0
#> 53fad3d8c52f4b5a023ff92aee9c0f1a     0     0    23      0     0      0      0
#> 5e37342d650cde6cfc5bbe2726b12e3e     0     0     0      0     0      0      0
#> 56025f74cc8e54d4a054442d786ef096     0     0     0     51     0      0      0
#> d9e51d98f1c86a1c426552d252c944e4     0     0    15      0     0      0      0
#> bc905dcc609da83d3620f1986e63a9f7     0     0     0      0     0      0      0
#> 90009df763e86e6eca7baa54f01551e6     0     0     0      0     0      0      0
#> 6bbd12b392461134369d835715cd2765     0     0     0      0     0      0      0
#> 83c88a72c2dab171c38aadd98aa2d389     0    13     0      0     0      0      0
#> 2b93708de9f1c24ff6140814a599a3a8     0     3     0      0     0      0      0
#> 671a53a2acbf0c8c7e240e446daf0000     0     0     0      0     0      0      0
#> bfe35644221421300c0c17e00e465c05     0    41   106     71    35     75     51
#> d12bdf20f24b92aa8ec94c6c3f60497f     0     0     0      0     0      0      0
#> 3bdb72ce060bd89f107c5c2bc95399ed     0     0     0      0     0      0      0
#> e6a82be4b75aac8b928c1dfbe9ee9287    37     0     0      0     0      0      0
#> e3a63ba0a0ec40fbf270a57c180e0d3c     0     0     0      0     0      0      0
#> 2877cf981b128a18f788ba5437c1afa6     0     0     0      0     0      0      0
#> c015bc50696fe3e3d94ace01d9dcd691     0     0     0      0     5      0      0
#> eb2eeddf829080f8d9f25c6be5d507cd     0     6     0      0     0     13      8
#> b1fcda6d8df4c25c8e6d7124c076b1fe     0     0     0      0     0      0      0
#> de8bb5c39fb121fd048dfe0d639e6e8a     0     0     0      0     0      0      0
#> a234ac03223c275785dedf6df3b2aeff     0     0     0      0     0      0      0
#> 9ea1641eeed1f51fc228c653076d1d08     0     0     0     70     0      0      0
#> 694b926113eb8291fea166446890c8fe     0    10     0      0     0      0      0
#> b2af163540d17007833dc819704afb28     0     0     0      0     0      0      0
#> 6fc0ba1f8ff8f9259c9d492e59771755     0     0     0      0     0      0      0
#> 7b8c1dba059eac7ec7f3cefa203f6b2d    69     0    52      0    61     63      0
#> 0b06028193ef1ce19f5339a8ba7d4448     0     0     0      0     0     22      0
#> 56b4ec05ef3bbded515ec3dd2d573197     0     0    70      0     0      0      0
#> fc9c84f8767611251f0eb93af4ca10db     0     0     0      0     0      0      0
#> d3616bc0d8925a452274c484c73406c3     0     0     0      0     0      0      0
#> 7986b8b93f3afce920bd19a5d2229c09     0     0     0      0     0      0      0
#> bc78ee7297b6912452be0865d113f7ad     0     0     0      0     0      0      0
#> c472668b46b0d0f575108aa4170dd81f     0    35     0      0     0      0      0
#> 2658a78c5ff06760389ce2e3a548b4a2     0     0    59      0    49      0      0
#> 600c9e81fd1d1ddd5a37a5aa0f9f796d     0     0    44      0     0      0      0
#> f8d0bb97f093d8e8630ea62ac8a1bf13     0     0     0      0     0      0      0
#> 6c685b83d63d6d56a9dcaba34215ae6b     0     0    18      0    21      0      0
#> e0e7e6208280a5846a8b66f00667f8ca     0     0   130     59    71      0      0
#> 6d1aa1c1f93c423a856a5dd3b62e5c70     0     0     0      0     0      0      0
#> e553e766ea0e08f1f9fa2debe7946af6     0     0     0      0     0      0      0
#> 0f1ef43f6b23a567aba80599766d70bc     0     0     0      0     0      0      0
#> 3c4ccaa45e8d148c9a17b5d0f3720cd2    98     0    73    107   136    205     51
#> d9df6265d046a558d87d397b8a177d54     0     4     0      0     0      0      0
#> 4bbad82a39655a26a2de61b88c5af7fa     0     0     0      0     0      0      0
#> f5a6865451a3b138f4915c5ab445a0a7     0    12     0      0     0      0      0
#> 4a86fc550fa1ac37aa480d37edea875f     0     0     0      0     0      0      0
#> 02fc21e72f66adea757b355edea4a369     0     0     0      0     0      0      0
#> 25cf10f4f15dab00925cef52b9f08442     0     0     0      0     0      0      0
#> fa8c8b98d197c252368c160b43cefb98     0     0    37      0     0      0      0
#> c0b07d29aed1a875641aede53bc0a20c     0     0     0      0     0      0      0
#> 341433c1c8cc65ce47dcf32179d08fd8     0     0     0      0     0      0      0
#> 8bf175f329d16e4aa732cf2b32279df3     0     0     0      0     0      0      0
#> 062f7552b747dccab9586dff7275b5d9     0     0     0      0     0      0      0
#> 059e124ad6ba8f98c3d21db28a075ea8   139    14   140    107    72     80      0
#> 3d7a81faea3f783158352c9f739c8d17     0     0     0      0     0      0      0
#> 90be2fb2c82a98002d33ed7531e877e3     0     0     0      0     0      0      0
#> b0e263cdd1c189fc1757d8189db14be2    30    90   232     32    67     70      0
#> b962ce332fea338b05293d1a5a29cbca     0     0    38      0     0      0      0
#> 7a310a8f8e4bcbfb717263bd5376994d     0     0     0      0     0      0      0
#> 57136cbdae8fe066c2bb4a661003719b     0     0     0      0     0      0      0
#> 69fe46cab1cc9407470a2397bfc84103     0     0     0      0     0      0      0
#> b88efba97658a6936bb053985d026ace     0     0     0      0     0      0      0
#> 4d691fbbe2dbfb2815649fb72c757977     0     0     0      0     0      0      0
#> 49192cc8a222d97b965b839ff5d0f333     0     0     0      0     0      0      0
#> 859873fcd088c294a92a13340744eb6f    10     0     0      0     0      0      0
#> 0a8c44ca4d744e313e6dc5b239fca562     0     0     0      0     0      9      0
#> 09f94345824c82e8fe0d824cc4478fb5     0     0     0      0     0      0      5
#> 9dfbdb741fb0dc25d2206b8680bbf579     0     0     0      0     0      0      0
#> dc2bea463c720874cccd189ccdcf0e3e     0     0     0      0     0      0      0
#> 6273eb52991d528cf5281215f79b40a3     0     0     0      0     0      0      0
#> dc4bcbe74986e44ec3a040b866dfda9e     0     0     0      0     0      0      0
#> 96ec4781ee61f4083792f252f1b6cd3a     0     0     0      0     0      0      0
#> 9657d52acca51e701e7a3b7cedcfcc12     0     0     0      0     0      0      0
#> f65f183293d4710041057ab98f2c94ed     7     0     0      0     0      0      0
#> 8f6b3783b31a2dd8640de3df92dbb9c4     0     0     0      0     0      0      0
#> 0d90eb4842c9ad9126ed33c538c768ac     0     0     0      0     0      0      0
#> 8d3ee46aaaf7728594014be61c649347     0     0     0     32     0      0      0
#> e011a95fda1f1f2bafc4e5b83b52d989     0     0     0      0     0      0      0
#> ef8ea65ba4c2e2da3ba9556a87d2385d     0     0     0      0     0      0      0
#> dbbb92eb3a8a3feaf96092d46f2b6a45     0     0     0      0    86      0      0
#> 06cb03cddd8a12dbfaf27e8984b395c4    65    58   254     56    86    128    107
#> cdf134bb501c54d997945dabd44e35ff     0     0     0      0     0      0      0
#> f3c4561cc01a3b45ee3f134d709c6b94     0     0     0    108     0      0      0
#> b121442fd1eefe78ce4f4602aac1842e     0    71     0      0    79      0      0
#> ba8707aa88a2a62939d9338189b2a733     0     0    17     20    11      0      0
#> ecf7289628d49f9244469337ee1622cb    25     0    44     46     0     32     18
#> 9f06b400b93f859b11874b0bd7f09295     0     0     0      0     0      0      0
#> f775be5cf7c0dd2177006a44913503d5     0     0     0      0     0      0      0
#> e87ceb10f2becb2a6c4cf0692b1923ee    38     0    37      0    39     26     13
#> 036574bb8e05f3b9a11ab7a4090d6950     0     0     0      0     0      0      0
#> 0267ba9db97e00ae01fc2cccbb6263e0     0     0     0      0     0      0      0
#> be4f701c1dbfbda16c5150af29430141     0    46    43      0    35     42     32
#> 655a66e6566798ade91f73b2942c0b63     0     0     0     18     0      0      0
#> cd79f2696509d0ad906f5addba264079     0     0     0      0     0      0      0
#> 07ea3284bda2bf7daf1251ad2744e4c4     0     0     0      0     0      0      0
#> e3870da9de4a2b18a7ae9e9cae208e3e     0     0     0      0     0      0      0
#> e00bd59c1ea0ca9b1cc5f976801d4e1b     0     0     0      0     0      0      0
#> b39b381d0f533c9923c3b898ffab5614     0     0     0      0     0      0      0
#> e2959f0552fa6df5707d9ec8b57816d6     0     0     9      0     0      0      0
#> 53157acd1642f973842d01e1cd287421     0     0     0      0     0      0      0
#> 5991790d1c0aa8b73b2d7f906449a1f9     0     0     0      0     0      0      0
#> adee3552c211d9e33aaa28e82532bdea     0     0     0      0     0      0      0
#> 0089040d041888e3ae4aa4e0010a0e36     0     0     0      0     0      0      0
#> 8e876e13d051992c4c8616cf2db7d73a     0     0     0      0     0      0      0
#> 179abbc5c1f1549ffd0d7ae6963e3783     0     0     0      0     0      0      0
#> 340d095f9a9fb445582b60e1282d985d     0     0    21      0     0      0      0
#> 3d7e5432d8eaa3fa56ab7138605ea43d     0     0     0      0     0      0      0
#> d4b233740db907d69035915e8657fa8d     0     0    16      0    10      0      0
#> b22119df3b1de76745abaac585b646f8     0     0     0      0     0      0      0
#> 64761a24fc66e7c93145217dd153b1ee     0     0     0      0     0      0      0
#> eb6d7f9e0601207fdcfc4a88947d513c     0     0     0      0     0      0      0
#> 0d20271e3e2dbc68e1de868599a9a1ae     0     0     0      0     0      0      0
#> 369ea22772cb774950c830d0af7e8235     0     0     0      0     0      0      0
#> cae63066cfe088a49fcac50f119c9e2d     0     0     0      0     0      0      0
#> d9c57ec20ed90e2cebf6cacc4f9a5511     0     0     0      0     0      0      0
#> f5828cc3078b4c219a44c67e96e5d455     0     0     0      0     0      0      0
#> 73847ac2778d4cba27fe9758ca129d6c     0     0     0      0     0      0      0
#> d78003f13712da48d0bb36fd45116b91    52     0     0      0     0      0      0
#> eb40853d986ad72c82658d939f291ff2     0     0     0      0     0     29      0
#> 5d675a3518222fc99c8cba34feaeac81     0     0    83      0    78     72      0
#> 3ba390ae11c4e76375af082eefe7e41f     0     0     0      0     0      0      0
#> 3aa66daa4e560cf62c1833697bcd3837     0     0     0      0     0      0      0
#> 6f4388510be5ea70356f2257d8c57ddf    64     0     0      0     0      0      0
#> b33390d034148edf07d8a742778ddcbd     0     0     0      0     0      0      0
#> 09225ec7044394bbd3b6811a6799e481     0     0     0      0     0      0      0
#> fe6fbf9b86e0f487911a856009fa4a1e     0    10     0      0     0     37      0
#> 0aee6329e0b5b81449980151bf40205a     0     0     0      0     0      0      0
#> 0ef66a78236f90916afc2305f3ff6c80     0     0    27     14    20     15      0
#> fb1eb41359005c2e3e308fb668f07220     0     0     0      0     0      0      0
#> cb59c57800227006aead644df42b8986     0     0     0      0     0      0      0
#> 138da9dc244c3292536cacb8df79467b     0     0     0      0    17      0      0
#> 72f92361ec0696cb8df985b93d277a25    23    44     0      0     0     50      0
#> 4d5ab121763879fd59e8f00d35473540     0     0     0      0     0      0      0
#> 20032b2b586fee22805eccebe3059170    31     0     0      0     0     31      0
#> 710324ccc8cafa6fc9b38baf8299d8aa     0     0     0      0     0      0      0
#> 523e6b60340f4a215877b1c336822e50     0     0     0      0     0      0      0
#> 4501f98c83484fd902be11383c8d032b     0    22     0      0     0      0      0
#> d4cbf23757aad57c91b058cbb05566e4     0     0     0     14     0      0      0
#> de7964d5c22112693ed5796bb9f04b0b     0     0     0      0     0      0      0
#> f4a283be42989b920b83ea73ec72b25b     0     0     0      0    49      0      0
#> 2f3111649d0850ab86799f8bbadc28f8     0     0    15      0    11      0      0
#> b8acdd45cf33ab75fe4c0851a1787dc5     0     0     0      0     0      0      0
#> 5bb7fec39744a6d4e4d61596dabc6886     0     0     0      0     0      0      0
#> e8eddfb37c4215c86813163a27f202a5     0     0     0      0     0      0      0
#> 826e188c911c0846fd927300127c80c0     0     0     0      0    25     23      0
#> ff1fe599c44074ee6bdf919a119df93c     0     0     5      0     0      0      0
#> aff00055d3876f624d67fb80bdb59934     0     0     0      0     0      0      0
#> 1beded194f3ae4b57715b45e50adef46     0     0     0      0     0      0      0
#> 908d91feba7c1f1d8cc8a7145173d942     0     0     0      0     0      0      0
#> e7095743695f6d4337d01ceb9e29e224     0     0     0      0     0      0      0
#> 4cd420962855afdf0bc69bb268b846c7    32    57   122     59    67     60      0
#> 430ba5bea0fc08f65b69497bfbeb5c4b     0     0     0      0     0      0      0
#> e12cb2f0d2ce4401413176007dcd4158     0     0    69     26    59      0     25
#> ec615213f6c7febb65c9f2e45e438111     0     0     0      0     0      0      0
#> 412dd52cd54d0244ba8fc42ad140acda    17     0     0      0     0      0      0
#> 89192d9f779341d313f52aad0dc6a06f     0     0     0     14     0      0      0
#> ed543f92c61a06fc4e703655ce910a2f     0     0     0      0     0      0      0
#> fea1a70b0ef0f72e3bf27d8b295331a8     0     0     0      0     0      0      0
#> 6afe805b36f5aa29e7f2096aac42ace1     0     0     0      0     0      0      0
#> d9cc5a46831f24ec9192ae82dbc38b0a     0     0     0      0     0      0      0
#> 0d3eb38eb502c4ad13bd5a7060295411     0     0     0      0     0      0      0
#> 69a0ff5538a5490e291c75ce70325066    17     0     0      0     0      0      0
#> f91d545c815effa3c5afa83a8ad3397d     0     0     0      0     0      0      0
#> e6b032dea8e13c34af958fc48af09b58     0     0     0      0     0      0      0
#> 84a1b89045ab6285328f366224f3ec87     0     0     0      0     0      0      0
#> 4a6f0076e4ac71474b3c364294ad3b33     0     0     0      0     0     24      0
#> 7bdd39e4e3ccc7000862b487510e67d6     0     0     0      0     0      0      0
#> 998aedf12cdf054f5e504a2e21d20498    37     0     0      0     0      0      0
#> 567c2640108006c8eca081ad4180e346    46     0     0      0     0      0      0
#> e50a6625b9176b6856764b213013150f     0     0     0      0     0    102      0
#> 8d5a5ab7c784b0325a4a86a16fbba80e    72     0     0      0     0      0      0
#> f58ae28a309fa3d5898d2eb46d9eb9e0     0     0     0      0     0      0      0
#> 7ac65a8fd21111b68e7c2d6205112a83     0     0     0      0     0      0      0
#> d7add64c8700f5a8cadf40c8c7f5ff67     0     0     0      0     0     11      0
#> 6bad9d981d3cf0f90262547180dc6371     0     0     0      0     0      0      0
#> 149c0e63c83113e1246c0e041de67d9c     0     0     0      0     0      0      0
#> 6f196fae54ea90b5c11627c503c98940     0     0     0      0     0      0      0
#> b67c995d842e5b4b92fa0ce052603ef1     0     0     0      0     0      0      0
#> f44ffae67c7297e2c696b75fbd4e0f4c     0     0    22      0     0      0      0
#> 33997912cd952feaec18fbd11b31acac     0     0     0      0     0     44      0
#> 0a035423c4bf71a17566549f9b69448d     0     0     0      0     0      0      0
#> f65a6f1f1e5a4bc5c733ba364c0096d3     0     0     0      0     0     18      0
#> f15389a9f245465c33db3de72aeb17ba    28     0     0      0     0      0      0
#> cd2b676d3c0785a60359c0d416701ee4     0     0     0      0     0      0      0
#> 70ba2a378ebcd9a8f57b024a561cc6b1     0     0     0      0     0      0      0
#> 0273707384afdd41e9d022b539094211    85    56     0      0     0      0      0
#> f9ac98e876ecdb504202388d7ccdda5f     0     0     0      0     0      0      0
#> bdd9b87b782b8f935ecac62a9824e6a1     0     0     0      0     0      0      0
#> 1f12b97495668f1fc378a349546d1806     0    24     0      0     0      0      0
#> 83f1b4d776f6e08a6ba6f86634a1e5a9     0     0     0      0     0      0      0
#> 49f614420449f6c7ecc8c9f4e137ad05     0     0     0      0     0      0      0
#> aba5aa00a1b262f27517eca64ee3333f     0     0     0      0     0      0      0
#> 8971980484eeabeef00145b93b699e30     0     0     0      0     0      0      0
#> c1347ae54b0fc1dd385a1ccf3b8061ce    22     0     0      0     0      0      0
#> e9fe19783ebcc7e43a3a07ca8225db62     0     0     0      0     0      0      0
#> 7c523a5f0a7c736efd82afef84137901     0     0     0      0     0     36      0
#> e13c07ec336c89db5183cc6e9e57ee6e    18     0     0      0     0      0      0
#> a017ce431ccadea609b501d46e64e55b     0     0     0     30     0      0      0
#> fadf38d03d4b944f5782a759591ec3ae     0     0     0      0     0      0      0
#> 41bfda1c7f32205246e73a1e3e3ec4f4     0     0     0      0     0      0      0
#> 9e4f3e697b70758a7a996f3dd4533650     0     0     0      0     0      0      0
#> 39c4f6caa7510ed108f3ebcf72f79edf     0     0     0      0     0      0      0
#> 4aa083fe7c1b0639dbc48ec649e68a31     0     0     0      0     0      0      0
#> fd2ac2e14b72b0e22a052e3ab4580471     0     0    10      0     9      0      0
#> 4ee10294799bcfeab0c8ee567d5d9608     0     0     0      0     0      0      0
#> 46869f12c31de7a05781e930c25a1ccd     0     0     0      0     0      0      0
#> 5e2548d6a970645d0ca943a72d207489     0     0     0     52     0      0      0
#> 0538e1beed69f098c3f602130ee11a1c     0     0     0      0     0      0      0
#> 991c563278c4e8c91ffebce3560fe51b     0     0     0      0     0      0      0
#> ce05f58164e61517198f7ab0e3304cb5     0     0    21     11    26     10      0
#> fdbd0fed21bae91c05b2da268a025a89     0     0     0      0     0      0      0
#> 7412d11595d77f26792071cd28a26760     0     0     0      0    19      0      0
#> 4b087e13957eb38128f63e1c99df23c3     0     0     0      0     0      0      0
#> 310a6f0e6841db2ce75005cac942a972     0     0     0      0     0      0      0
#> 3df0fada803ac582802b47d4f4efe0c1     0     0     0      0     0      0      0
#> 3cd4e34c02ffdd1446df65bb54278c57     0     0     0      0     0      0      0
#> b14a1dea93554257791f16c2e94906b0     0     0     0      0     0      0      0
#> 902a4435fc3af101e8a9b5f51550fa4d     0     0     0      0     0      0      0
#> 1995185b31349529f58964c2e5a3becd     0     0     0      0     0      0      0
#> 69b247c565afdeacee0240482b1ce536     0     0     0      0     0      0      0
#> 95ff1b9528acca018e9c4697a8625504     0     0     0      0     0     46      0
#> 132d26d3ce6ac2ab25ccb4dd171879d7     0     0   127      0     0      0      0
#> bf1a6e4f94ee8ba60cf79089139641b6     0     0     0      0     0      0      0
#> 704189aa09f64b992b3a02ee001dc706     0     0     0      0     0      0      0
#> dea5ea27839ee67a31d18ba0bf0f7f65     0     0     0      0     0     25      0
#> 3840a4fff0250232f8570cae7c3514ba     0     0     0      0     0      0      0
#> 31ecabc032997aaf41dd0d8def3ec93f     0     0     0      0     0      0      0
#> 024a0e0348ed9abe9ae1c48b1b7ba09e     0     0     0      0     0      0      0
#> fc1d9940419420113ce3fbfacc8d703a    68     0    85     79   205     15      0
#> 0770fbf28cd307a08bf557af7ffecfd3     0     0     0      0     0      0      0
#> 883cbe8e4d1d98f468355e8ea974d222     0     0     0      0     0      0      0
#> cb6d6561383fb68f43f92a88a607749d     0     0     0      0     0      0      0
#> e07e0f1a9fbc696a9fe4f9b7739ae45b     0     0    17      0     0      0      0
#> 73d8d8f999acfa13cba2277f7526897e     0     0     0      0     0      0      0
#> 5218274c354d41965abbf4b0751d6f3a     0     0     0      0     0      0      0
#> 5c31570f706752fa4f7851aca5b5c291     0     0    10      0     0      0      0
#> 6a5e71c8b86a52ab1dd5e61164103940     0    25    67      0    39     61      0
#> 25641b0b803c7d8493ebc02cb104bf3d     0     0     0      0     0      0      0
#> 8844d70dd1fdb96f2a7faa65bccf78c6     0     0     0      0     7      0      0
#> 371b855afa8bed030e86ec78241b01a3     0     2     2      0     0      0      0
#> d0c180a47378d414a1d088ac14db236c     0     0     0      0     0      6      0
#> a072a588b3940fab3b227e4901c40f7f     0     0    20      0     0      0      0
#> 977a3523051e6b23da6221dbde795b34     0     0     0      0     0      0      0
#> 1c3adc87d2e93b08e576953508eb680d     0     0     8      0     0      0      0
#> 74eb0ddf379b00c5032b18fb9fbf1f32     0    26    29      7    19     44     12
#> c335b34a15db7d532f511e6c952abe62     0     5    31      2    26      0      0
#> 0c78285397fb853a0e39069213510c54     0     0     0      0     0      0      0
#> 779dec6de9c6a8b7df28fc8ec3c75bee     0     0     0      0     0      0      0
#> 5526624b473bb1f0541dc6fb95c98678     0     4     0      0     0      0      0
#> 649e7fa7b9ce5c5472d34675ed643639     0     0     0      0     0      8      0
#> a9c5addeea542953d1bae4d30d599b02     6     0    19     20    24     16      0
#> 68291fb3b2558d438da4810961c8d53a     0     0     0      0     0      0      0
#> 9c5d92bfefc7a190b93d5b9a6a77697b     9     0    23      9    22     24      0
#> 3ddab42d701754625cab0c236c5b6e6d     0     6     0      0     0      9      0
#> 287e840f9ecaef80d7a4762c8d6bf401     0    15     0      0     0     24      0
#> 885f411fa54c5576148678f665bb053f     0     0     0      0     0      0      0
#> b5709a57df9cbc8a741bb78a36f29c06    98   132   102     81    89    162     35
#> 26698cc60893017ef0b1ab173688d3de    35    25    56      0    38     14      0
#> db64f0028b6518839ab7129c4b7ef72d   104   213   173     67   147    271     75
#> f46fb43262da299ec7f0084298cafa50     0     0     0      0     0      0      0
#> c45e7087de669be270881fc8ace0fad8    35    50    55     41    57     94     25
#> 854a20d5019388e81b1e78419eb017fd     0    28     0      0    15     41      0
#> 9a2325ea015a8c5944047e3349cc1565     0    19    88     19    81     38     12
#> c5305ae3c6adc90814d70ab48f4f071a     0     0     0      0     0      0      0
#> e230123ec3ce8654942ee1a000a80010     0    48   102      0    82    106     22
#> 14e3c60cb885509667670151bc4a1df0     0     0     0      0    17      0      0
#> 7fd7b7bece92bbbfcceb2dcf1e835a8e     0     0     0      0     0      8      0
#> c098d1fb769dd975c2332dfb70da489e     0     0     0      0     0      0      0
#> 95dcdf1eb06b252e67dcdbdaffc505e8     0     0     8      0     0      0      0
#> c5399d0f9eda8f0815fd28a257952496     0     0    62      0    26      0      0
#> c7388cc314e992819b8f97ea35280982     0     0    35      0     0      0      0
#> ae8facc45d372f7108bd9499f74d1bfc    15     0    26      0    15     16      0
#> 16caf6f5bc653fb6ba668eb821076a21     0     0     0      0     0      0      0
#> e7e283b72ee15bcf0877c5eeb6137ebc     0     0     0      0     0      0      0
#> 1b3a779b384c73b33207dccd14b14fe6     0     0     0      0     0      0      0
#> c8d1afdf7b395ffc22193c5dc2fd5dc3     0     0     7     21     0      0      0
#> 0e3ffc2a68378086e2480831acd1a025     0     0     0      0    23      0      0
#> 0ddf7339b608a340252b33f63fcc0019    34    17    31     16    16     42      0
#> 5fd4b2f7dee3fb0bde6f856ae266dbcb     0     0     0      0     0      0      0
#> e2f4786bc66410c07bd8b96b5d0e2d92     0     0    34      0    40     19      0
#> c30aae8f5216c9536492d0beae74950f     0    96   315     74   221    260     17
#> 9087ecee0806a80861788a17dd7e6809     0     0     0      0     0      0      0
#> c29bc08418318360592b2b2717376a0a     0     0     0      0     0      0      0
#> f791c7db1a78d35e84803b5e03c01965     0     0     0      0     0      0      0
#> 44b44fc4055a58c01f6d1e9de5e0d1d9   106     0   110    149     0     41     58
#> 342df9401103d274851d3bd80704508d     0     0     0      0     0      0      0
#> 0590e6672bb2c55c0ff54999a2b6dcfe     0     0     0      0     0      0      0
#> d4e9147608bf883fdbbf9d23b17028d1     0     0    14      0     0      0      0
#> db24794eb8667aed224401386223f671    10    87    47     10    22     25      0
#> fe41e6fd9770030a21418b288930710f     0     0     0      0     0      0      0
#> 699f18f3b155bff97de941f2cfbbbadc     0     0    62     72    39     14      0
#> 13870cf319997750e3b96bbb0b2f1aea     0     0     0      0     0      9      0
#> 7ca4fb1b1740406241acb862445cf311     0     0    17      0    15     42     13
#> 9c0569fb95451d009fd3754303eaed12    32     0     0      0     0      0      0
#> 1ac26a3c12288d8fbadd6aa61366b1b7     0     0     0      0    36      0      0
#> 89f3ba99552605501b8acb816d288afe     0    29     0      0     0      0      0
#> b4051b5227746f04c1feb9d8979299ff     0     0     4      0     0      0      0
#> 087d358aaf77f81aaf50a623828865bb     0     0     0      0     0      0      0
#> e451c4b3ff1c3dabeee296529b3a7a07     0     0    15      0     0      0      6
#> 012871784f77c9dca7d446ddef2048d1     0     0     0      0     0      0      0
#> 95c6ce2788783bf2f83343deba23b8bd     0     0     0      0     5      0      0
#> 213227f244d97553a9a14c4e3fd04c98     0     0     0      0     0      0      0
#> c0e5e55c5cfbd576e7c085f1e4b69638     0     0     0      0     0      0      0
#> 51e98451d7e3f46127481ec5490f1c91     0     0     0     65     0      0      0
#> c75285a6c61d8fce20d823769c25c12a     0     0     0      0     0      0      0
#> 177966f7619eb24fdbce2266f5bed39a     0     0    16      0     0      0      0
#> b89347f0534679f890628364d2c7b4fe     0     0     0      0     0      0      0
#> 8dff1cb12c0beb5bc7b2546480a6a6c7     0    28     0      0     0     42      0
#> 84f2777cfcb0fbed6281917fcc765022     0     0    38      0     0     12      0
#> a5e4d106a2ab947984fdef187996a0b9     0     0     0      0     0      0      0
#> b6cbf6077b64acbb718971a1ffc97a27     0     0     0      0     0      0      0
#> 27e29495fda424728fd6a80b83a16772     0     0    46      0    50      0     14
#> 63c44ab3725dd0b803d43a29db9fd7fd    17     0    31      0     0      0      9
#> 9c931ab5bf18ef9b10996d59d292e6f5     0     0    40      0    44     46     24
#> 12fd8efb82c3d4909543f6db9853134e     0     0     8      5     9      4      6
#> a8a60c20f7ab7d0f2dadd9f8858c19be     0     0     0      0     0      0      0
#> d57e09eb1645f1abba6074d3bd7d4c52     0     0     0      0     0      0      0
#> 7e0f2b982ca2d35b6ac2cf251d78d3c4     0     0     0      0     0      0      0
#> 3e2e8e938a68a0318bcf7145ad5ed4a9    19     0    20      0     0      0      0
#> f2117b49bf57d9451bbef6de5b762483     0     0     0      6     0     15      0
#> be24dfe003725a0b0f28f4136106faf8     0     0     0      0    24      0      0
#> e63301d7d0af55d5c43db051dc51522f     0     0    29     45    18      0      0
#> a6e2e5501507562f2174f65e606b2d43    30     6    44      0     0     47     20
#> a06aa93cb50b9a39941a11c7242dde96    11     0     0      0     0      0      0
#> 6a0a8abb36aa3bfdee14edbe97e9d781     0     0     0      0     0      0      0
#> 4d3def5ffad0c392381591440bfc3668     0     0     0      0     0      0      0
#> ce437b81e7cf5ce986f7240c26464b0a    29     0    15      0     0      0      0
#> e57788f9d2fae76ed9ca566ad4d475a8     0     0     5      0     0      0      0
#> 8c525d0f63137bd068b75c57f5d04c09     0     0     0      0     0      0      0
#> 6d670aa40dafa152c2d59cf2be9bef74     0     0     0      0     0      4      0
#> e08f8e9e9ffb10a702f2ef617f70eecf     0     0     0      0     0      0      0
#> 7def1c807ca0767ccdfe5a6001dc51f1     0     0     0      0     0      0      0
#> 3b9a17275147adf1dc8c42eb00a06a36    20     0     0      0     0      8      0
#> 24be0b7e4c85ac53b09fddbca916ab96    14     0     0      0     0      0      0
#> f521a47d3bb4acfec6b261f880b4eea8     0     0     0      0     0      0      4
#> 91c24c9b6ffbedcc257d164176060fe0     0     0     0      0     0      0      0
#> 9b45f210647c5fe71b005eb9e25299a5     0     0     0      0     6      7      0
#> 974077d25fa6b32a7af22670b142a8f1     0     0     2      0     0      0      0
#> 203aa2535e54a11e4c48af4ab946686b     0     0    10      0     0      0      0
#> 9064afe569448558eb08e5bf5e458bce     0     0     0      0     4      0      0
#> fad432eeeff7f1d923e6d60aafe18dcd     0     0     0      0     0      0      0
#> a29db0e5f4a04c5c6937745ccccb167b     0     0    48     12    22     18      0
#> f179c10f8e86873c90ceb14b974a114e     0     0     0      0     0      0      0
#> a23e8588a2240d7af6f65519158c902c     0     0     0      0     0      0      0
#> 52034e0da6d5e7962369437956984e18     0     0     0      0     0      0      0
#> 522d20c27cf4237b60f36ffdca056c19     0    13     0      0     0      0      0
#> 957615aefb421e6ad4728e610b2bdfa6     0     0     4      0     0      0      2
#> 5ccc74e70b8a82a6eebea3385738ece3     0     5     0      0     0      0      0
#> e5924fd11a48419c6ec0e2a7cb4e43ee     0     0     0      0     0      0      0
#> eef7be4b2531cb7e909b58e91af5640d     0     0     6      0     0      0      0
#> 4416cbc2026260fa4e2ed0a8feff1443     0     0     0      0     0      0      0
#> 03f6091f5aacb69f06fc925acd1c9506     0     0     0      0     0      4      0
#> f655e5ddf5074d5e71bfa583926dec30     0     8    21      0    14      0      0
#> 9f733d39e588c7310499dd444cdd7812     0     0     0      0     0      0      0
#> 2b413de5e05a1385bab229f50879628f     0     0     0      0     0      0      0
#> 178b596f2fe5cb6cefc5513917988d84     0     0     7      0     0      0      0
#> aa94b014c726cf40983045bd606241c1     0     0     3      0     0      0      0
#> 55310b3aac4868c4e33b0a9492d7609a     0     0     0      0     0      0      0
#> c5d8d828df0e1baf5c2db05a6b8c4476     0     0     9      0     0      0      0
#> 6ceb26f60606233a7728d3457af84380     0     0     0      0     0      0      0
#> 1ff83868d1c9ecf52d9b776e0b3fd896     0     0    11      0     0      0      0
#> 3edfccfb9f81b340854e84ee5102a4c1     6     0     0     10     0      6      0
#> 963dc552c551872800f2fcbc5663e5fd     0     0     0      0    23      0      0
#> 3fda8c20e5da9fbc3fb7c5c864559019     0     0     0      0     0      8      0
#> 6fd88fa7e883bdead3c5ba89da853e8d     0     0     0      0     0      8      0
#> 35826fe826b7a9f202e8bef7b39d9a47    14     0     0      0     0      0      0
#> cc514f614ba146f489bf62b5a56074c5    35    80    74     25    48     37      0
#> 26260de3858081ce55e6f220b3ba25df     0     0     0      0     0     26      0
#> d6727a85782ffd2a069dccc23adaf521     0    68   191     56   135     94      0
#> eb1a74c2ca8831ef945775de0827ddda     0     0     0      2     0      0      0
#> 4ed3218615bbc1ed4da82d9806e4786b    11     0     0      0     0      0      0
#> 8a3367e36ed4e416d5a77a7946cf6fc5     0     0     9      0     0      0      0
#> 56f45eba0496ed29fb7171041a2a950d    11     0    30      5    26     24      0
#> 9d2782f19e60ec0cf72e1dc46b1d1766     4     0     0      0     0      0      0
#> 1559f8b64dc2bdc8cdc26ff8a5ffa4cb     0     0     0      0     0      0      0
#> 6d01f384a1cce140edea17dccd0f116d     0     0     0      0     0      0      0
#> 5950949b71082741302b712e97309b17     0     0     0      0     0      0      0
#> f6f238f95cc97995c9bf619b0c620b33     0     0     0      0     0      0      0
#> c7a112ffadf8009efd318625da44b18e     0     0     0      0     0      0      0
#> e79f4f30b0dacfcc50f3c1b12804ca3f     0     0    13      0     0      0      0
#> 297c825dc268f9a7b86b402edb6961b3     4     0     0      0     0      0      8
#> a1ed56ad4e86e0b4adc11baf2c4ab94b     0     0     0      0     0      0      0
#> 8afccb8f7f50c5603b471fb0c808c9e6     0     0     0      0     0      0      0
#> ce41c973bd4965a400989ec128890d55     4     0     0      0     0     11      4
#> 5d99ba689ad34e4393a905e8a725b72b     0     0     0      0     0      0      0
#> bfaa2d35e6dd13b9b32de811c6f92592     0     0     0      0     4      0      0
#> 0c544d5f3a49a010971a96d8c39e4905     0     0     0      0     0      0      0
#> bab3f08e1d75f42987012cae451639ff   482    65   186     65   145    142     29
#> fa97a81b8c33569e2aa1fc9a3137b435     0     0     0      0     7     16      0
#> a4d08886c3c6c62efd5678339a07828b     0     0     0      0     0      0      0
#> 11ed0972601847ec168473a6bc4eeba7     0     0     0     17    32      0      0
#> a6e4d88aa4b40402b1f94760d86b729c     0     0     0      0     0      8      0
#> 86deaff47d33946ede1b7c25c6bbe9f9     0     0     0      0     0      0      0
#> a1d3f1f8c151eaeb5cf12179861e0c87     0     0     0      0     0      0      0
#> 3980891c22ccbb387b65f9e437b2f94f     0     0     0      0     0      0      0
#> fc8ec321e0c2fc2a6de9abe98256b3f5     0   179   133     70   137      0      0
#> 84c11e6ffb29ebc2192f4b049b7dbddb     0     0     0      0     0      0      0
#> 82104585f0b2617c13eb69f0bbdf7d5a     9     0     0      0     0      0      0
#> 61cd4918428d3b36910b9cb2fc45946a     0     0     0      0     0      0      0
#> 0ab9a976935482bce5f8af4cc6685147    41     0    11      0    12      0      0
#> 87be3450216b45c53851fd173582a2d0     0     0     0      0     0      0      0
#> 93a51bb503fadf90121f518db980e158     0     0     0      0     0      5      0
#> b70cd843161996e075a853900e699d0f     0     0     0      0     0     28      0
#> 154b8be9eea158f93d27d314d8e0e2c8     0     0     0      0     0     34      3
#> 70691c2737cf7dc7c194de819f0e6cda     0    18     8      0     0     11     18
#> 2ae1d80ead92f5991f58da14d000acd6     0     0     0      0     0      0      0
#> 44fa299492152518fec34b9fa4554648     0     0     0     11     0      0      0
#> e15c8f9df87f363b0dbbf49de0bed249     0     0     0      0     0      0      0
#> 73689cad6eb4e86d9d0bb48c0da24463     0     0     0      0     0      0      0
#> 234e15e82d64f2e1d61802fa8ea1998a     0     0    35      0    30     62      0
#> 88c76ddb1b24acffe65d8128e2f48244     0     0     0      0     0     29     14
#> 5544139a16716cc6ef48b2c16c9613e9     0     0     0      0     0      0      0
#> 5b7fe68483e701059142cb531397b547     0     0     0      0     0      0      0
#> 0a0875a0acf4a28d1bc1c4d1d11aaa64     0     0     0      0     0      0      0
#> 6b47ec9d1c70f07d16834e9a4e54aa70     0     0     0      0     0      0      0
#> f8752b4824d239cfb19f3216cc34941b     0     0     0      0     0      0      0
#> ec8a0768cd2fb426d7e6bfc2c4f50e5b     0     0     0      0     0      0      0
#> 8a6b88f9515bf2dfac09e6a69045a4fe     0     0     7      0     0      0      0
#> f5151c67af19c8471743ef3e2914f6b1     0     0     6      0     0      0      0
#> c3842619fcbde52781c512c2561fca0b    19     0     0      0     0      0      0
#> b3b00faa9dee9fa6b3b33744f51da00e     0     0     0      0     0      0      0
#> c4f04e28c6267ac6c5b05d403e8f91c3     0     0     0      0     0      0      0
#> 7fc06f107a6ddf6486ef2d43fcf23626     0     0     0      0     0      0      0
#> 8d1dffcd377b338d2a0388637507b592     0     0     7      0     0      0      0
#> 15410d762889d51818986bfcdb5d8965     0     0     0      0     0      0      0
#> 48f68e44e258de63575b8bf6cf565d52     0     0    23      4    16      0      0
#> 4712befe3685661e0fc0968859c369d5     0     0    18     20     0      0      0
#> be9a3727a0422aea147100370d046fa8    21    26    84     41   100     67      3
#> 3c5942781761830f7b5406f551574424     0     0     0      0     0      0      0
#> 291d9f4d3d3c5cc3b0a1804770ed1e77     0     0     0      0     0      0      0
#> 784764c519a64dc9a147025a23222e93     0     0     9      7    24     14      0
#> b9e955ea5254dd68b54775aad2a67829     0     0     0      0     0      0      0
#> 8197f1594e8a17f3fc989fea64f54b8f     0     0     0      0     0      0      0
#> bf4a4f2383a16d2ecf27bf8e162b9ba4     0     0     0      0     0      0      0
#> 4e7d22d793d30bdedbd32a5b73b564b8     9    28    14      0     0     24      4
#> e88959b739c8cceb94b94b7887dd1df3     0     0     0      0     0      0      2
#> 5403445694339e2ae5ecac06d9f3ff47     0     0     6      0     0      0      0
#> ba3943fb075bbf63fb699638292826a1     0     0     0     10     0      0      0
#> 67b2d836db093cead5a8286142669613     0     0     0      0     0      0      0
#> 6c3a0ea164d3eb098dd2baac5519b38f     0     0     0      0     0      0      0
#> c643a776c36bcc00dfb5d790e1ea3dfe     0     8    21      0     0      0      0
#> d53cc81a9d57a64d194c876b76bd66e0     0     0     0      0     0      0      0
#> 7dc00d40971bb454b39fdd9b2dd9aa25     0     0     0      0     0      0      0
#> b4d0cec8543e940be1b5096370fdfb80     4     0     0      0     0      0      0
#> 6c54ceb94bc5e8fbd65224cffc60990a     0    15    24     20     0      0      0
#> b12ee9ca5a2320957213a114788cf2b7     0     0     0      0     0      0      0
#> 5e040c2397f959e7f0dfc6ec2855b812     0     0     0      0     0     19      0
#> 3105a594bba2604ae5ed19f5e0495d2e    13     0     0      0    13      0      0
#> d75c040bdaa1cca96cd0655c9af7ab28     0     0     0      0     0     30      0
#> 9b909edaebabc90b9c1a703626e6bd0c     0     0     0      0     0      0      0
#> a6593da81035cf121dbe46346a4ea960     0     0     0      0     0      0      0
#> 60e4c88400a750a8122b97476313a67e     0     0     0      0    28      0      0
#> 7854a78e3b881cee1c63ea0a50ba29ce     0     0     0      0     0      0      0
#> 1796734b56d45bdc997ccc09ab8731d1    12     0     0     10     0      0      0
#> 8243cd669ae7327f4326be78429ab05a     0     0    30      3     5      0      0
#> ebb93dc32c63402289114859d476cc5b     0     0     8      0     0      0      0
#> 8a01af5a531e851e2c50fc01131cda33     9    13    22      8    15     18      0
#> 03d2cc3f7b2bc352752d39e3738d0af9     0     0     0      0     0      0      0
#> ec378af6e32960b11b5f7391da4ae4fb     0     0     0      0     0      0      0
#> b407f2de868c064f65eb5ad8f82a8d0d     0     0     0      0     0      0      0
#> 73a95b6ba34fcd260989d848fe15d377     0     0     0      0     0      0      0
#> 3abb68202d15781590a2d1f9c5b37e15     0     0     0      0     0      0      0
#> 77923ef8e9ed01c9b6a8ac85eb028163     0     0    39      0     0      0      0
#> 9981b8fd20f400e052527f84f5db1b66     0     0     8      0     0      0      0
#> 597bdec49adf7675912f9fba4e262e70     0     0     3      0     0      0      0
#> 0fd3126118b13fb187bc55476a67c95a     0     4     0      0     0      0      0
#> 2a0c511ecc2513c42fe624e8e99d07da     0     0     0      0     0      0      0
#> d2912fea2dcf2e5b97e56c6a7d24295d     0     4     0      0     0      0      0
#> 2976af6aaf4261b6c274fcc3229a4d65     0     0     0      0     0      0      0
#> f097c45ce419da576aa111c98a4c7272     0     0     0      0     0      0      0
#> 3abc8260d34008f263e232ea25426c3f     0     0     3      0     0      0      0
#> d7e5b6432dc178aaeed9c81e0ec57afb     0     0     0      0     0      0      0
#> dfe71ce78ddfdbc7b453d17d628365f2    17    20    46     12    30     36      0
#> 289b18a86084273ce3b4038bf37e2840     0     0     0      0     0      0      0
#> 844a67032c52dc49019bf5dff248bf9c     0     0     3      0     0      0      0
#> 38cb3e3ce5ac5aed7dec50247ecfb267     0     0     0      0     0      0      0
#> 74e88e71e9929b1ddfcafccfab3b7bba     0     0     0      0     0      0      0
#> 66bc0737fc7fa1463ac7189ab7073998     0    42     0     17     0      0      0
#> b5a7aa588c9cae6ae8ca2f7d87ac53e8     0     0     0      0     0      0      0
#> db0540920457b21b170947aac9137771     0     0     0      0     0      0      0
#> ebbaf260510ab16ac6b8c2b9866d05b9     6     0    33      0     3      0      0
#> c1f461b67419d4a80a61db0867da9ed6     0     8    22      0    28     24      0
#> 0c6eba6452635eecfb865a1049299176     3     0     0     10     0      5      0
#> 7ecabb10f29a0c44eb53e9f9852fe22a     0     0     0      0     0      0      0
#> a33cb3037e7825920073b4dddaaa41e6     0     0    30      0    28     21      0
#> d625812c15c074c2e2b18099f97412c9     0     0     0      0     0      0      0
#> 17ef77f12e2a55717e910a5f6a625448     7     0    28      0     4      0      2
#> 7baf638fdd3f71166aafd22ecf410348     0     0     0      0     9      0      0
#> 9436cda211f60f722865cd23ef75d4ed     0     0     9      0     0      0      0
#> c8ab0b20e740eb6b35e68d58f69736e4     0    33   129     54    99     40     37
#> 5e43f2d8be71399fc8d6736754cb42ac    18     0     0      0     0      0      0
#> 19ad90ee1081a003406b02179ae63709     0     0     0      0     0      0      0
#> 0a25cdbd0660e68970c563c43a18f54a    30     0    38     20    21      0      0
#> 0aa6b210e902d8ab9d6f267407258ebb     0     0     0      0     0      0     18
#> a864c36d0d898cc5abd94d6d6f0bd5be     0     0    50      0    15      0      0
#> dea0dcb305ff0f27e3c36c7133cf933d     0     0     0      0     0      0      0
#> 463977baa69c7724fdc9407ad162714c     0     0     0      0     0      0      0
#> fcb72326c29342bcfb052953466eafc6    26     0    20     27     0     36      0
#> 5a554353f9ed728f0ab0631f29c8be2e     0     0     0      0     0      0      0
#> a9dd27b48e443415a743bf9902f515c4     0     0     0      0     0      0      0
#> 77a40bcf979c0099c52f19a3f1fe5217     0     0     0      0     0      0      0
#> b4007ac525a9b409d70426e24421cff0     0    19    19      0     0      7     18
#> 406ba18a175c4c79efb25f48a6beafe4     0     0     0      0     0      0      0
#> fe1c166755bc722cefd2e2a3b2b70672     0     0    18      0     0      0      0
#> 8afd274849ac07c0ca66f8f5f143df3c     0     0    15      0     0      0      0
#> 197301ed9d8728b5c061b99a17973161    24     0    69      0    39     15      0
#> 14779f9dc3505598ad3d799e4ee676c4     0     0     0     38     0      0      0
#> 4085c246640eee043d523ce2d799023a    32     0    28     35     0      0      0
#> bba09a403e67c1381de01c86511f1057     0     0    28      0    30     11      0
#> 392c59fee4e6b91d7654f2d3e3c65be8    21     0     0      0     0     18      0
#> d2b05ae73ad5ed296d7add9e2bb2cce4     0     0     0      0     0      8      0
#> 8306a9473dabfd5918ea869dbfcd1291     0    15     0      0     0     12      0
#> 06064c372bb624503c2905d13f7691f5   262    57   171    221   103    101     85
#> 3265692341e42dbd5cf4ed1938b906e3     0     0     0      0     0      0      0
#> da96a0797b9881a123795838b9f01140     0     0     0      0     0      0      0
#> 0bdddf10e577657becec88e42001fd16     0    33     0      0     0      0      0
#> 168999a950ebe17e6d2ea1f90e557d91     0     0    18     32     0      0      0
#> a266ffe9d38f626d106fa615da4618ef     0     0     0      0     0      0      0
#> 7706405d5bce7715f01d8b54188def42     0     0    23      0     0      0      0
#> 4b0cfc8899b466bd9992f7b22a5a06be     0     0    11      0     0      0      0
#> b81b2b4a7ca9eb2c23b170648dfaf22d     0    17     0      0     0     14      0
#> 75cf40872201f150cb5b191676f64b78    42    39    62     33    54     42     27
#> c8e5794c58de9326084a62d3282d257e     0    26     0      0     0      0      0
#> b0706a571f169c0bac7b4864052fcaa7     0     0     0      0     0      0      0
#> 0388f969e3527b2938acc676d6757bdb     0     0     0      0     0      0      0
#> 2dfe751390aedfbd9f62dbb77d1f10b4     0     0     0      0     0     10      0
#> 2548dc9697e1f9afe04a9012f40e2f0a     0     0     0      0     0      0      0
#> 2f019615d5141a4ab00b1aab98882a33     0     0   118     29    42      0     14
#> 629e4f0cf1826904ff74907b2507f58b     0     0     6      0     0      4      0
#> 774861f348cec8451102730c21a790b0     0     0    29      0     0      0      0
#> 11c0d6de28ab026a9ae17abe4efdd6bf     0     0    16      0     0      8      0
#> c11381dda0841c0d08cbfcef3a760aa3     0     0     0      0     0      0      0
#> 3c9619d4f9db4ac6ad7151a0c9d86ad3     0     0     0      0     0      0      0
#> 7e2c72204a8c45b1b0ac39e0e6047098     0     0    43     16    17      0      0
#> 8254e0f9f07709b62dd8483f3b25fad0    35     0    60     23    33     37     20
#> e78135f2886d427e678d12c446bf070e     0     0     0      0     0      0      0
#> abec938a8d53e7471f3f476d50694be9     0     0     0      0     0      0      0
#> 4963567d2979bb761291aa98e0b8595b    18     0     0      0    15     28     16
#> d4de9d1453102564f35262faa3952266     0     0     0      0     0      0      0
#> efa4a78c4b993e9d115841c8a02c795e     0     0     0      0     0      0      0
#> fcacd7424bf3543322ea3e9842d14e86     0     0    35     13    37     32      0
#> f941e15ea7f7560191d57abc4d2bedd8     0     7     0      9     0     19      0
#> ac294f5242de3c9b26e75502c48f6383     0     0     0      0     0      0      0
#> 7b6a8e3f8bfa1bcd1b5fc55db048d005     0     0     0      0     0      0      0
#> c9e776ab599fe75db286675b698c9067     0    47    55      0     0     16      0
#> ebe5cb945f6d6d849c07424591f80753     0     0     0      0     0      0      0
#> f116387f77802409103a559e4205df06     0     0    20      0     0     19     33
#> 93f9939b261cddc9a57071b42771c1d4     0     0    35      0     0      0      0
#> 971a5f48b3418449ca7cf3e190a5e831    27     0    90     42    59      0      0
#> 730de992ca06a95e7193a93146a97cf2     0    26     0      0     0      0     23
#> 8fa8b4ac69784c26901e036d13ab2332     0     0     0      0     0      0      0
#> c363c9efaa85023b488b732b97e47270     0     0     0     27     0     12     20
#> 2e2181cc893ca1e1e677e888bf8af441     0     0     0      0     0      0      0
#> 86fa662da6869634e7f90f85e9ee4933     0     0     0     23     0      0      0
#> 9a02177f2c235a6eb62984d5aab151ff     0     0     0      0     0      0      0
#> d6854be498e22b65da3c00572ee178f0     0     0     0      0     0      0      0
#> 015f364f150b230450da0a123e7bb8ab     0     0     0      0     0      0      0
#> af2c1b22758cd6653e8e841497d096e7     0     0    35      0    43      0      0
#> e7f2fab8c0e7c365768a43ef71406605    33     0    18     44     0      0      0
#> fb57cf192bce1876471f2f6181bf4b55     0     0     0      0     0      0      0
#> 6c740ce69c9f45bea684dc9e36054006     0     0     0      0     0      0      0
#> 3b27c7bf279535d8e39664d8dc1d0daf     0     0     0     17     0      0      0
#> ceea3c29ad012c2918fb6aa21a44e8cb     0     0     0      0     0      0      0
#> fa99ecdb27c1cc44fd866ee1da250c4c     0     0     0      0     0     23      0
#> 387cc83573db0e444e81c89a8a4f8f9c    53    33    63     29    93     89     22
#> c7a8d1af5f7f472671f700775de3e0a2     0     0     0      0     0      0      0
#> 98a6a84550e07ea252cfd8a0f51042d4   217    26    36      0     0     19      0
#> 08b809a68126c38db786cafe9710b013     0     0     7     25     0      0      0
#> 8db4cc0a5d2f03fce8efd2b8374ff359     0     0    10      0     0      5      0
#> 3e7c7d38423c93b43a1d9865cf0bc4a6     0    30    57      0    40      0     15
#> 6f02050760c030ed8cd1ed836697fa6b     0    29    30      0    20     32      0
#> 8fc1865cfcae32f72c1b1570299e2ef3     0     0     0      0     0     13      0
#> 38cc138fac95b006929361093d840a00    15     0     6     11    14     32      0
#> f5f8eb63e22362509e2c7c78271305f2    65    15     0      6     0      0      0
#> 0bb238aaee3ccb7dc8dfa86dbc1ba9f7     0     0     2      0     0      0      0
#> ca4546caf07edf5200cf79cea04aebcd     0     0     0      0     0      0      0
#> 893d12123eb7f623489a4f025a15ae7f    22     0     0      0     0      0      0
#> 192c5a817f48368339e519e0a0a907da     0     0     0      0     0      0      0
#> 3467185ca92b7bec2523214605adabf0     3     0     0      0     0      0      0
#> a2aae80d8fb922fe1074d243b59ce237    28     0    22     22    21     52      0
#> f11233f16d4177d38ab1c0bfc649e056     0     0     0      0     0      0      0
#> b219ca7793ab9f5e441f2a68bc76152c     0     0     6      0     0     20      0
#> e6e1848b51e705e55b3805703f088311    76     0     0      0     0      0      0
#> bd074b9bf365d11468bd76268d63b15b    22    22    30     11    18     40     20
#> f381b492d6651c52fd09ce9351e7eb2d     0     0     0      0     0      0      0
#> d9b31e1c30facc24fa778b1d65f5c457     0     0   811    251   797    666     17
#> 14a67c8c8cc31d203eac8109443fab5f     0     0    12      0     0      0      0
#> 026ad0ea53fb202ebb770acf65705d6e     9     0    24     32    19      0      0
#> 19ba9fea04969f975d6dc8d3da20d5e6     0     0    43      8    18      0      0
#> 0fc4426ed7dccc48ce27cfff8083433d   317    26   458    267   266    224     36
#> 45ad90e179842705fa984962f131fd61     0     0     0      0     0      0      0
#> d284afab226d22fd0596221d17db23d3     0    23    49     39    18     29     15
#> 38194655ba17dbfd7fcc1ed42a57c169     0    55    43     12    41     47     27
#> e30be5268bbebe597424bdf171cc1607     0     9     0      0     0      0      0
#> cb7ccdb4956134b452a3cb78d1751c57     0     0     0      0     0      0      0
#> b585436677ac88dd4ac6239655f2a268     0     0     0      0     0      0      0
#> f5d770eff5c9b9b0c859db57d175eb34     0     0    19      6     0      0      0
#> 1b45ffbc88c788dcee2239e68901cf05     0     0     0      0     0      0      0
#> 3fd951cf7cc3a0740a1f20f8a9c6d2c2     7     0    12     15     0      0      0
#> 2027d3e90ccc4003f1954322e99db6a9     0     0    31      0    23      0      0
#> 6157cb92771830955413f808e83fc420    13     0    23     30    13     12      5
#> 6349506cac32056a214b04239bdc3a33    60    42    97     50    66    167     43
#> 2254a8d37cf537b1e35097d70b6ac91b   155    67   206    118   127    161    102
#> 38f3bc8ef8836574de0fcbeafd348b65     0     0     0      0     0      0      0
#> fec23c70ae0c70ecceea1487cd978ee5     0     0    36      0    23      0      0
#> 503d123a6c78c6fb280a6472bd61b51e     0     0     9      0     4      0      0
#> 9f33de8570c639e727465f1ee4d750f6     0     0     0      0     0      0      0
#> d68e1899fa65515135d9bd8d356b7e3d     0     0     0      0     0      0      0
#> e15b88ccac7a4bef64215e915f6e0d89    15     0    23      0     7      2      6
#> 767042b4b19d8bfaaa426dad7dbce378     0    10     0      0     0      0      0
#> 4a3c4abee502d6185abd750522566bf3     0     0     0      5     0      0      0
#> c254aa5c00b7f27d756f7bbe96b3d4ad     0     0     0      0     0      0      0
#> fe8d8af791eb6e5d4335bb7604c74e9d     0     0     0      0     0      0      0
#> d81800d7e31b69877ce1eab9858b8f3a     0     0     0      0     0      0      0
#> c55c7cb80c5700ff9aded76192543b5a     0     0     0      0     0      0      0
#> 91dee7b3b90bbf2e53b2472b9bbde0e9     0     0     0      0     0      0      0
#> 6b28df61cb5dee74b76b1a13aaa0027d     0     0     0      0     0      6      0
#> 9a830fa26205ad91a02b5bcaebb0221a    10    16    29     15    44     32      0
#> a54b8142d2a09d058b81c377e8e286b8     0     0     0      0     0      0      0
#> 69ef3a03a4d733dbffd99fdaa1d2e62a     0     0     0      0     0      0      0
#> 083cbb131dd305a88b7d338bc371c35f    14     0    23      0    17      0      0
#> ef1a92caffb4d7cae286fdc77b0a93a6     0     0     0      0     0      0      0
#> 2cddd5dca166486f50a0828b28001e1e     0     0     0      0    11      0      0
#> bb899d1dcfd89cb200a76cb8951904f8     0     0    24      0     0     32      0
#> 7817f6851e9df03a3d992cf915310dba     0     0     0      0     0      0      0
#> 54eaf6571de743a86ccb531bee8aad7e     0     0     4      0     0      0      0
#>                                  4.18RO 2.13RI 5.14RI 5.14RO 6.13RO 1.8RI
#> 374a5eb6496d14e8e9d6f3d8c7e34d9d      0      0      0      0      0     0
#> 4c47d5cf81df1ac6e2a0560029a81578      0      0      0      0      0     0
#> 8fb6741a6685fad0374f85ad3e16156c      0      0      0      0      0     0
#> 310f3de009c95de5c938b8e811dfc96f      0      0     35      0      0     0
#> 7d3b39e8c0f0fdcb81b0c898a98c38f9      0      0      0      3      0     0
#> fc211549300b0954dbb4a4bf57e8a605      0      0      0      0      0     0
#> 53fad3d8c52f4b5a023ff92aee9c0f1a      0      0      0      0      0     0
#> 5e37342d650cde6cfc5bbe2726b12e3e      0      0      0      0      0     0
#> 56025f74cc8e54d4a054442d786ef096      0      0      0      0      0     0
#> d9e51d98f1c86a1c426552d252c944e4      0      0      0      0      0     0
#> bc905dcc609da83d3620f1986e63a9f7      0      0      0      0      0     0
#> 90009df763e86e6eca7baa54f01551e6      0      0      9      0      0     0
#> 6bbd12b392461134369d835715cd2765      0      0      0      0      0     0
#> 83c88a72c2dab171c38aadd98aa2d389      0      0     18      0      0     0
#> 2b93708de9f1c24ff6140814a599a3a8      0      0      0      0      0     0
#> 671a53a2acbf0c8c7e240e446daf0000      0      0      0      0      0     0
#> bfe35644221421300c0c17e00e465c05     31      0     57     60      0    46
#> d12bdf20f24b92aa8ec94c6c3f60497f      0      0      0      0      0     0
#> 3bdb72ce060bd89f107c5c2bc95399ed      0      0      0      0      0     0
#> e6a82be4b75aac8b928c1dfbe9ee9287      0      0     54      0      0     0
#> e3a63ba0a0ec40fbf270a57c180e0d3c      0      0      0      0      0     0
#> 2877cf981b128a18f788ba5437c1afa6      0      0      0     32      0    23
#> c015bc50696fe3e3d94ace01d9dcd691      0      0      0      0      0     0
#> eb2eeddf829080f8d9f25c6be5d507cd      0      0      0      0      0     0
#> b1fcda6d8df4c25c8e6d7124c076b1fe      0      0      0      0      0     0
#> de8bb5c39fb121fd048dfe0d639e6e8a      0      0      0      0      0     0
#> a234ac03223c275785dedf6df3b2aeff      0      0      0      0      0     0
#> 9ea1641eeed1f51fc228c653076d1d08      0      0      0      0      0     0
#> 694b926113eb8291fea166446890c8fe      0      0      0      0      0     0
#> b2af163540d17007833dc819704afb28      0      0      0      0      0     0
#> 6fc0ba1f8ff8f9259c9d492e59771755      0      0      0      0      0     0
#> 7b8c1dba059eac7ec7f3cefa203f6b2d      0     41      0      0      0    37
#> 0b06028193ef1ce19f5339a8ba7d4448      0      0      0      0      0     0
#> 56b4ec05ef3bbded515ec3dd2d573197      0      0      0      0      0     0
#> fc9c84f8767611251f0eb93af4ca10db      0      0      0      0      0     0
#> d3616bc0d8925a452274c484c73406c3      0      0     38    132      0     0
#> 7986b8b93f3afce920bd19a5d2229c09      0      0      0      0      0     0
#> bc78ee7297b6912452be0865d113f7ad      0      0      0      0      0     0
#> c472668b46b0d0f575108aa4170dd81f      0      0     26      0      0     0
#> 2658a78c5ff06760389ce2e3a548b4a2      0      0     48     54      0     0
#> 600c9e81fd1d1ddd5a37a5aa0f9f796d      0      0      0      0      0     0
#> f8d0bb97f093d8e8630ea62ac8a1bf13      0      0      0      0      0     0
#> 6c685b83d63d6d56a9dcaba34215ae6b      0      0      0      0      0     0
#> e0e7e6208280a5846a8b66f00667f8ca      0      0     62     87      0    38
#> 6d1aa1c1f93c423a856a5dd3b62e5c70      0      0      9      0      0     0
#> e553e766ea0e08f1f9fa2debe7946af6      0      0      0      0      0     0
#> 0f1ef43f6b23a567aba80599766d70bc      0      0     27      0      0     0
#> 3c4ccaa45e8d148c9a17b5d0f3720cd2    139    105      0     62      0   102
#> d9df6265d046a558d87d397b8a177d54      0      0      0      0      0     0
#> 4bbad82a39655a26a2de61b88c5af7fa      0      0    115      0      0     0
#> f5a6865451a3b138f4915c5ab445a0a7      0      8      0      0      0     0
#> 4a86fc550fa1ac37aa480d37edea875f      0     26      0      0      0     0
#> 02fc21e72f66adea757b355edea4a369      0      0      0      0      0     0
#> 25cf10f4f15dab00925cef52b9f08442      0      0      0      0      0     0
#> fa8c8b98d197c252368c160b43cefb98      0      0      0      0      0     0
#> c0b07d29aed1a875641aede53bc0a20c      0      0      0      0      0     0
#> 341433c1c8cc65ce47dcf32179d08fd8      0      0      0      0      0     0
#> 8bf175f329d16e4aa732cf2b32279df3      0      0      0      0      0     0
#> 062f7552b747dccab9586dff7275b5d9      0      0     14      0      0     0
#> 059e124ad6ba8f98c3d21db28a075ea8     66     67     89     70     50    82
#> 3d7a81faea3f783158352c9f739c8d17      0      0     68      0      0     0
#> 90be2fb2c82a98002d33ed7531e877e3      0      0      0      0      0     0
#> b0e263cdd1c189fc1757d8189db14be2      0      0     89    198      0    77
#> b962ce332fea338b05293d1a5a29cbca      0      0     19      0      0     0
#> 7a310a8f8e4bcbfb717263bd5376994d      0      0      0      7      0     0
#> 57136cbdae8fe066c2bb4a661003719b      0      0      0      0      0     0
#> 69fe46cab1cc9407470a2397bfc84103      0      0      0     11      0     0
#> b88efba97658a6936bb053985d026ace      0      0      0      0      0     0
#> 4d691fbbe2dbfb2815649fb72c757977      0      0      0      0      0     0
#> 49192cc8a222d97b965b839ff5d0f333      0      0      0      0      0     0
#> 859873fcd088c294a92a13340744eb6f      0      0      0      0      0     0
#> 0a8c44ca4d744e313e6dc5b239fca562      0      0      0      0      0     0
#> 09f94345824c82e8fe0d824cc4478fb5      0      0      0      0      0     0
#> 9dfbdb741fb0dc25d2206b8680bbf579      0      0      0      0      0     0
#> dc2bea463c720874cccd189ccdcf0e3e      0      0      0      0      0     0
#> 6273eb52991d528cf5281215f79b40a3      0      0      0      0      0     0
#> dc4bcbe74986e44ec3a040b866dfda9e      0      0      0      0      0     0
#> 96ec4781ee61f4083792f252f1b6cd3a      0      0      0      0      0     0
#> 9657d52acca51e701e7a3b7cedcfcc12      0      0      0      0      0    52
#> f65f183293d4710041057ab98f2c94ed      0      0      0      0      0     0
#> 8f6b3783b31a2dd8640de3df92dbb9c4      0      0      0      0      0     0
#> 0d90eb4842c9ad9126ed33c538c768ac      0      0      0      0      0     0
#> 8d3ee46aaaf7728594014be61c649347      0      0     15      0      0     0
#> e011a95fda1f1f2bafc4e5b83b52d989      0      0      7      0      0     0
#> ef8ea65ba4c2e2da3ba9556a87d2385d      0      0      0      0      0     0
#> dbbb92eb3a8a3feaf96092d46f2b6a45      0      0      0      0      0     0
#> 06cb03cddd8a12dbfaf27e8984b395c4     69     79     84    246      0    88
#> cdf134bb501c54d997945dabd44e35ff      0      0      0      0      0     0
#> f3c4561cc01a3b45ee3f134d709c6b94      0      0      0      0      0     0
#> b121442fd1eefe78ce4f4602aac1842e      0      0      0      0      0     0
#> ba8707aa88a2a62939d9338189b2a733      0      0      0      0      0     0
#> ecf7289628d49f9244469337ee1622cb     30      0     23     53      0     0
#> 9f06b400b93f859b11874b0bd7f09295      0      0      0      0      0     0
#> f775be5cf7c0dd2177006a44913503d5      0      0      0      0      0     0
#> e87ceb10f2becb2a6c4cf0692b1923ee      0     38     27     31      0    33
#> 036574bb8e05f3b9a11ab7a4090d6950      0      0      0      0      0     0
#> 0267ba9db97e00ae01fc2cccbb6263e0      0      0      0      0      0     0
#> be4f701c1dbfbda16c5150af29430141      0      0     47      0      0     0
#> 655a66e6566798ade91f73b2942c0b63      0      0      0      0      0     0
#> cd79f2696509d0ad906f5addba264079      0      0      0      0      0     0
#> 07ea3284bda2bf7daf1251ad2744e4c4      0      0      0      0      0     0
#> e3870da9de4a2b18a7ae9e9cae208e3e      0      0      0      0      0     0
#> e00bd59c1ea0ca9b1cc5f976801d4e1b      0      0     21      0      0     0
#> b39b381d0f533c9923c3b898ffab5614      0      0      0      0      0     0
#> e2959f0552fa6df5707d9ec8b57816d6      0      0      0      0      0     0
#> 53157acd1642f973842d01e1cd287421      0      0      0      0      0     0
#> 5991790d1c0aa8b73b2d7f906449a1f9      0      0     28      0      0     0
#> adee3552c211d9e33aaa28e82532bdea      0      0     11      0      0     0
#> 0089040d041888e3ae4aa4e0010a0e36      0      0      0      0      0     0
#> 8e876e13d051992c4c8616cf2db7d73a      0      0      0      0      0     0
#> 179abbc5c1f1549ffd0d7ae6963e3783      0      0      0      0      0     0
#> 340d095f9a9fb445582b60e1282d985d      0      0      0      0      0     0
#> 3d7e5432d8eaa3fa56ab7138605ea43d      0      0      0     16      0     0
#> d4b233740db907d69035915e8657fa8d      0      0      0      0      0     0
#> b22119df3b1de76745abaac585b646f8      0      0      0     10      0     0
#> 64761a24fc66e7c93145217dd153b1ee      0      0     23      0      0     0
#> eb6d7f9e0601207fdcfc4a88947d513c      0      0      0      0      0     0
#> 0d20271e3e2dbc68e1de868599a9a1ae      0      0      0      0      0     0
#> 369ea22772cb774950c830d0af7e8235      0      0      0      0      0     0
#> cae63066cfe088a49fcac50f119c9e2d      0      0      0      0      0     0
#> d9c57ec20ed90e2cebf6cacc4f9a5511      0      0      0      0      0     0
#> f5828cc3078b4c219a44c67e96e5d455      0      0      0      0      0     0
#> 73847ac2778d4cba27fe9758ca129d6c      0      0      0      0      0     0
#> d78003f13712da48d0bb36fd45116b91      0      0      0      0      0     0
#> eb40853d986ad72c82658d939f291ff2      0      0      0      0      0     0
#> 5d675a3518222fc99c8cba34feaeac81      0      0      0      0      0     0
#> 3ba390ae11c4e76375af082eefe7e41f      0      0      0      0      0     0
#> 3aa66daa4e560cf62c1833697bcd3837      0      0      0      0      0     0
#> 6f4388510be5ea70356f2257d8c57ddf      0      0      0      0      0     0
#> b33390d034148edf07d8a742778ddcbd      0      0      0      0      0     0
#> 09225ec7044394bbd3b6811a6799e481      0      0      0      0      0     0
#> fe6fbf9b86e0f487911a856009fa4a1e      0      0      0      0      0     0
#> 0aee6329e0b5b81449980151bf40205a      0      0      0     57      0     0
#> 0ef66a78236f90916afc2305f3ff6c80      0      0      0     24      0    33
#> fb1eb41359005c2e3e308fb668f07220      0      0      0      0      0    49
#> cb59c57800227006aead644df42b8986      0      0      0      0      0     0
#> 138da9dc244c3292536cacb8df79467b      0      0      0      0      0     0
#> 72f92361ec0696cb8df985b93d277a25      0     27      0      0      0     0
#> 4d5ab121763879fd59e8f00d35473540      0      0     16      0      0     0
#> 20032b2b586fee22805eccebe3059170      0      0     56      0      0     0
#> 710324ccc8cafa6fc9b38baf8299d8aa      0      0      0     59      0     0
#> 523e6b60340f4a215877b1c336822e50      0      0      0      0      0     0
#> 4501f98c83484fd902be11383c8d032b      0      0     11      0      0     0
#> d4cbf23757aad57c91b058cbb05566e4      0      0      0      0      0     0
#> de7964d5c22112693ed5796bb9f04b0b      0      0      0      0      0     0
#> f4a283be42989b920b83ea73ec72b25b      0      0      0      0      0     0
#> 2f3111649d0850ab86799f8bbadc28f8      0      0      0      0      0     0
#> b8acdd45cf33ab75fe4c0851a1787dc5      0      0      0      0      0     0
#> 5bb7fec39744a6d4e4d61596dabc6886      0      0      0      0      0    14
#> e8eddfb37c4215c86813163a27f202a5      0      0      0      0      0     0
#> 826e188c911c0846fd927300127c80c0      0     20     22     18      0     0
#> ff1fe599c44074ee6bdf919a119df93c      0      0      0      0      0     0
#> aff00055d3876f624d67fb80bdb59934      0      0     10      0      0     0
#> 1beded194f3ae4b57715b45e50adef46      0      0      0      0      0     0
#> 908d91feba7c1f1d8cc8a7145173d942      0      0      0      0      0     0
#> e7095743695f6d4337d01ceb9e29e224      0      0      8      4      0     0
#> 4cd420962855afdf0bc69bb268b846c7      0      0     56    121     31    64
#> 430ba5bea0fc08f65b69497bfbeb5c4b      0      0      0      0      0     0
#> e12cb2f0d2ce4401413176007dcd4158      9      0      0     80      0    29
#> ec615213f6c7febb65c9f2e45e438111      0      0      0      0      0     0
#> 412dd52cd54d0244ba8fc42ad140acda      0      0      0      0      0     0
#> 89192d9f779341d313f52aad0dc6a06f      0      0      0      0      0     0
#> ed543f92c61a06fc4e703655ce910a2f      0      0      0      0      0     0
#> fea1a70b0ef0f72e3bf27d8b295331a8      0      0      0      0      0     0
#> 6afe805b36f5aa29e7f2096aac42ace1      0      0      0      0      0    17
#> d9cc5a46831f24ec9192ae82dbc38b0a      0      0      0      0      0     0
#> 0d3eb38eb502c4ad13bd5a7060295411      0      0      0      0      0     0
#> 69a0ff5538a5490e291c75ce70325066      0      0      0      0      0     0
#> f91d545c815effa3c5afa83a8ad3397d      0      0      0      0      0     0
#> e6b032dea8e13c34af958fc48af09b58      0      0      0      0      0     0
#> 84a1b89045ab6285328f366224f3ec87      0      0      0      0      0     0
#> 4a6f0076e4ac71474b3c364294ad3b33      0      0      0      0      0     0
#> 7bdd39e4e3ccc7000862b487510e67d6      0      0      0      0      0    20
#> 998aedf12cdf054f5e504a2e21d20498      0      0      0      0      0    12
#> 567c2640108006c8eca081ad4180e346      0      0      0      0      0     0
#> e50a6625b9176b6856764b213013150f      0     51      0      0      0     0
#> 8d5a5ab7c784b0325a4a86a16fbba80e      0      0     29      0      0    32
#> f58ae28a309fa3d5898d2eb46d9eb9e0      0      0      0      0      0     0
#> 7ac65a8fd21111b68e7c2d6205112a83      0      0      0      0      0     0
#> d7add64c8700f5a8cadf40c8c7f5ff67      0      0      0      0      0     0
#> 6bad9d981d3cf0f90262547180dc6371      0      0      0      0      0     0
#> 149c0e63c83113e1246c0e041de67d9c      0      0      0      0      0     0
#> 6f196fae54ea90b5c11627c503c98940      0      0     29      0      0     0
#> b67c995d842e5b4b92fa0ce052603ef1      0      0      0      0      0     0
#> f44ffae67c7297e2c696b75fbd4e0f4c      0      0      0      0      0     0
#> 33997912cd952feaec18fbd11b31acac      0      0      0      0      0    33
#> 0a035423c4bf71a17566549f9b69448d      0      0      3      0      0     0
#> f65a6f1f1e5a4bc5c733ba364c0096d3      0      0      0      0      0     0
#> f15389a9f245465c33db3de72aeb17ba      0      0      0      0      0     0
#> cd2b676d3c0785a60359c0d416701ee4      0      0      0      0      0     0
#> 70ba2a378ebcd9a8f57b024a561cc6b1      0     16      0      0      0     0
#> 0273707384afdd41e9d022b539094211      0      0      0      0      0     0
#> f9ac98e876ecdb504202388d7ccdda5f      0      0      0      0      0     0
#> bdd9b87b782b8f935ecac62a9824e6a1      0      0      0      0      0     0
#> 1f12b97495668f1fc378a349546d1806      0      0      0      0      0     0
#> 83f1b4d776f6e08a6ba6f86634a1e5a9      0      0      0      0      0     0
#> 49f614420449f6c7ecc8c9f4e137ad05      0      0      0      0      0     0
#> aba5aa00a1b262f27517eca64ee3333f      0      0      0      0      0     0
#> 8971980484eeabeef00145b93b699e30      0      0      0      0      0     0
#> c1347ae54b0fc1dd385a1ccf3b8061ce      0      0      0     31      0     0
#> e9fe19783ebcc7e43a3a07ca8225db62      0      0      0      0      0     0
#> 7c523a5f0a7c736efd82afef84137901      0      0      0      0      0     0
#> e13c07ec336c89db5183cc6e9e57ee6e      0      0      0      0      0     0
#> a017ce431ccadea609b501d46e64e55b      0      0     35      0      0     0
#> fadf38d03d4b944f5782a759591ec3ae      0      0      0      0     38     0
#> 41bfda1c7f32205246e73a1e3e3ec4f4      0      0      0      0      0     0
#> 9e4f3e697b70758a7a996f3dd4533650      0      0      0      0      0     0
#> 39c4f6caa7510ed108f3ebcf72f79edf      0      0     30      0      0     0
#> 4aa083fe7c1b0639dbc48ec649e68a31      0      0      0      0      0     0
#> fd2ac2e14b72b0e22a052e3ab4580471      0      0      0     20      0    17
#> 4ee10294799bcfeab0c8ee567d5d9608      0      0      0      0      0     0
#> 46869f12c31de7a05781e930c25a1ccd      0      0      0      0      0     0
#> 5e2548d6a970645d0ca943a72d207489      0      0      0      0      0     0
#> 0538e1beed69f098c3f602130ee11a1c      0      0      0     83      0     0
#> 991c563278c4e8c91ffebce3560fe51b      0      0     20      0      0     0
#> ce05f58164e61517198f7ab0e3304cb5      0      0     10     24      0    12
#> fdbd0fed21bae91c05b2da268a025a89      0      0      0      0      0     0
#> 7412d11595d77f26792071cd28a26760      0      0      0      0      0     0
#> 4b087e13957eb38128f63e1c99df23c3      0      0      0      0      0     0
#> 310a6f0e6841db2ce75005cac942a972      0      0      0      0      0     0
#> 3df0fada803ac582802b47d4f4efe0c1      0      0      0      0      0     0
#> 3cd4e34c02ffdd1446df65bb54278c57      0      0      0      0      0     0
#> b14a1dea93554257791f16c2e94906b0      0      0      0      0      0     0
#> 902a4435fc3af101e8a9b5f51550fa4d      0      0      0      0      0     0
#> 1995185b31349529f58964c2e5a3becd      0      0      0     20      0     0
#> 69b247c565afdeacee0240482b1ce536      0      0      0      0      0     0
#> 95ff1b9528acca018e9c4697a8625504     26      0      0      0      0     0
#> 132d26d3ce6ac2ab25ccb4dd171879d7      0      0      0      0      0     0
#> bf1a6e4f94ee8ba60cf79089139641b6      0      0      0      0      0     0
#> 704189aa09f64b992b3a02ee001dc706      0      0      0      0      0     0
#> dea5ea27839ee67a31d18ba0bf0f7f65      0      0      0      0      0    34
#> 3840a4fff0250232f8570cae7c3514ba      0      0      0      0      0     0
#> 31ecabc032997aaf41dd0d8def3ec93f      0      0      0      0      0     5
#> 024a0e0348ed9abe9ae1c48b1b7ba09e      0      0      0      0      0     0
#> fc1d9940419420113ce3fbfacc8d703a      0     16      0     34      0    85
#> 0770fbf28cd307a08bf557af7ffecfd3      0      0      0      0      0     0
#> 883cbe8e4d1d98f468355e8ea974d222      0      0      0      0      0     0
#> cb6d6561383fb68f43f92a88a607749d      0      0      0      0      0     0
#> e07e0f1a9fbc696a9fe4f9b7739ae45b      0      0      0      0      0     0
#> 73d8d8f999acfa13cba2277f7526897e      0      0      0     29      0     0
#> 5218274c354d41965abbf4b0751d6f3a      0      0      0      0      0     0
#> 5c31570f706752fa4f7851aca5b5c291      0      0      0      7      0     5
#> 6a5e71c8b86a52ab1dd5e61164103940      6     28     50      0      0     0
#> 25641b0b803c7d8493ebc02cb104bf3d      0      0      0      0      0     0
#> 8844d70dd1fdb96f2a7faa65bccf78c6      0      0      0      0      0     0
#> 371b855afa8bed030e86ec78241b01a3      0      0      8      0      0     0
#> d0c180a47378d414a1d088ac14db236c      0      0      0      0      0     0
#> a072a588b3940fab3b227e4901c40f7f      0      0      0      0      0     0
#> 977a3523051e6b23da6221dbde795b34      0      0      0      0      0     0
#> 1c3adc87d2e93b08e576953508eb680d      0      0      0      5      0     0
#> 74eb0ddf379b00c5032b18fb9fbf1f32      0     41     50     18      0     0
#> c335b34a15db7d532f511e6c952abe62      0      0     49      0      0     0
#> 0c78285397fb853a0e39069213510c54      0      0      0      8      0     0
#> 779dec6de9c6a8b7df28fc8ec3c75bee      0      0      0      0      0     0
#> 5526624b473bb1f0541dc6fb95c98678      0      0      2      0      0     0
#> 649e7fa7b9ce5c5472d34675ed643639      3      0     35      0      0     0
#> a9c5addeea542953d1bae4d30d599b02      0     12      7     40      0     0
#> 68291fb3b2558d438da4810961c8d53a      0      0      0      0      0     0
#> 9c5d92bfefc7a190b93d5b9a6a77697b      3     18     35     26      0     0
#> 3ddab42d701754625cab0c236c5b6e6d      0      8     14      0      0     0
#> 287e840f9ecaef80d7a4762c8d6bf401      0      0      0      0      0     0
#> 885f411fa54c5576148678f665bb053f      0      0      0      0      0     0
#> b5709a57df9cbc8a741bb78a36f29c06     18    118    286    143      0    63
#> 26698cc60893017ef0b1ab173688d3de      0      0     26     34      0     0
#> db64f0028b6518839ab7129c4b7ef72d     23    194    424    203      0    65
#> f46fb43262da299ec7f0084298cafa50      0      0      0      0      0     0
#> c45e7087de669be270881fc8ace0fad8      0     85    201     44      0    16
#> 854a20d5019388e81b1e78419eb017fd      0      0     49      0      0     0
#> 9a2325ea015a8c5944047e3349cc1565      0     69     65     50      0    53
#> c5305ae3c6adc90814d70ab48f4f071a      0      0      0      0      0     0
#> e230123ec3ce8654942ee1a000a80010     27     60     44     60      0    33
#> 14e3c60cb885509667670151bc4a1df0      0      0     12      0      0     0
#> 7fd7b7bece92bbbfcceb2dcf1e835a8e      0      0      0      0      0     0
#> c098d1fb769dd975c2332dfb70da489e      0      0      0      0      0     0
#> 95dcdf1eb06b252e67dcdbdaffc505e8      0      0      0      0      0     0
#> c5399d0f9eda8f0815fd28a257952496      0      0     56     17      0     0
#> c7388cc314e992819b8f97ea35280982      0      0     16      0      0    32
#> ae8facc45d372f7108bd9499f74d1bfc      0     19     19      0      0     0
#> 16caf6f5bc653fb6ba668eb821076a21      0      0      0      0      0     0
#> e7e283b72ee15bcf0877c5eeb6137ebc      0      0      0      0      0     0
#> 1b3a779b384c73b33207dccd14b14fe6      0      0      0      0      0     0
#> c8d1afdf7b395ffc22193c5dc2fd5dc3      0      0      0      0      0     0
#> 0e3ffc2a68378086e2480831acd1a025      0      0      0      0      0     0
#> 0ddf7339b608a340252b33f63fcc0019      0      0     27     22      0    11
#> 5fd4b2f7dee3fb0bde6f856ae266dbcb      0      0     14      0      0     0
#> e2f4786bc66410c07bd8b96b5d0e2d92      0      0      0     65      0     0
#> c30aae8f5216c9536492d0beae74950f     38    126    165    114      0    88
#> 9087ecee0806a80861788a17dd7e6809      0      0     24      0      0     0
#> c29bc08418318360592b2b2717376a0a      0      0     12      0      0     0
#> f791c7db1a78d35e84803b5e03c01965      0      0     14      0      0     0
#> 44b44fc4055a58c01f6d1e9de5e0d1d9      0     67    146      0     45    72
#> 342df9401103d274851d3bd80704508d      0      0      0      0      0     0
#> 0590e6672bb2c55c0ff54999a2b6dcfe      0      0      0      0      0     0
#> d4e9147608bf883fdbbf9d23b17028d1      0      0      0     18      0    10
#> db24794eb8667aed224401386223f671      0     16    261     14      0    16
#> fe41e6fd9770030a21418b288930710f      0      0      0      0      0     0
#> 699f18f3b155bff97de941f2cfbbbadc     21    119     86      0      0     0
#> 13870cf319997750e3b96bbb0b2f1aea      0      0      0      0      0    26
#> 7ca4fb1b1740406241acb862445cf311      0      0     59     34      0     0
#> 9c0569fb95451d009fd3754303eaed12      0      0      0      0      0     0
#> 1ac26a3c12288d8fbadd6aa61366b1b7      0      0      0      0      0     0
#> 89f3ba99552605501b8acb816d288afe      0      0     37      0      0     0
#> b4051b5227746f04c1feb9d8979299ff      0      0      0      0      0     0
#> 087d358aaf77f81aaf50a623828865bb      0      0     10      0      0     0
#> e451c4b3ff1c3dabeee296529b3a7a07      0      0      5     12      0     6
#> 012871784f77c9dca7d446ddef2048d1      0      0      0      0      0     0
#> 95c6ce2788783bf2f83343deba23b8bd      0      0      0      0      0     0
#> 213227f244d97553a9a14c4e3fd04c98      0      0      0      0      0     0
#> c0e5e55c5cfbd576e7c085f1e4b69638      0      0      0      0      0     0
#> 51e98451d7e3f46127481ec5490f1c91      0      0      0      0      0     0
#> c75285a6c61d8fce20d823769c25c12a      0      0      0      0      0     0
#> 177966f7619eb24fdbce2266f5bed39a      0      0      0     10      0     0
#> b89347f0534679f890628364d2c7b4fe      0      0      0      0      0     0
#> 8dff1cb12c0beb5bc7b2546480a6a6c7      0      0     32     15      0     0
#> 84f2777cfcb0fbed6281917fcc765022      0      0      0     16      0    77
#> a5e4d106a2ab947984fdef187996a0b9      0      0      0      0      0     0
#> b6cbf6077b64acbb718971a1ffc97a27      0      0     52      0      0     0
#> 27e29495fda424728fd6a80b83a16772      0      0     17     27      0    31
#> 63c44ab3725dd0b803d43a29db9fd7fd      0      0      0     25      0     0
#> 9c931ab5bf18ef9b10996d59d292e6f5     19     22     42     70      0    43
#> 12fd8efb82c3d4909543f6db9853134e      0      0      0     20      0     0
#> a8a60c20f7ab7d0f2dadd9f8858c19be      0      0      3      0      0     0
#> d57e09eb1645f1abba6074d3bd7d4c52      0      0      0      0      0     0
#> 7e0f2b982ca2d35b6ac2cf251d78d3c4      0      0      0      0      0     0
#> 3e2e8e938a68a0318bcf7145ad5ed4a9      0      0     20      0      0     8
#> f2117b49bf57d9451bbef6de5b762483      0      0      0      0      0     0
#> be24dfe003725a0b0f28f4136106faf8      0      9      0      0      0     0
#> e63301d7d0af55d5c43db051dc51522f     10      0      0      0      0    10
#> a6e2e5501507562f2174f65e606b2d43      0      0     42     35      0    18
#> a06aa93cb50b9a39941a11c7242dde96      0      0     19      2      0     4
#> 6a0a8abb36aa3bfdee14edbe97e9d781      0      0      6      0      0     6
#> 4d3def5ffad0c392381591440bfc3668      0      0      0      0      0     0
#> ce437b81e7cf5ce986f7240c26464b0a      0      5     10      0      0     6
#> e57788f9d2fae76ed9ca566ad4d475a8      0      0      0      0      0     0
#> 8c525d0f63137bd068b75c57f5d04c09      0      0      3      0      0     0
#> 6d670aa40dafa152c2d59cf2be9bef74      0      0      0      0      0     0
#> e08f8e9e9ffb10a702f2ef617f70eecf      0      0      0      0      0     0
#> 7def1c807ca0767ccdfe5a6001dc51f1      0      0      0      0      0     0
#> 3b9a17275147adf1dc8c42eb00a06a36      0      0      0      0      0     0
#> 24be0b7e4c85ac53b09fddbca916ab96      0      0     10      0      0     0
#> f521a47d3bb4acfec6b261f880b4eea8      0      0      0      0      0     0
#> 91c24c9b6ffbedcc257d164176060fe0      0      0      0      0      0     0
#> 9b45f210647c5fe71b005eb9e25299a5      0      0      0      4      0     2
#> 974077d25fa6b32a7af22670b142a8f1      0      0      0      0      0     0
#> 203aa2535e54a11e4c48af4ab946686b      0      0      0      9      0     6
#> 9064afe569448558eb08e5bf5e458bce      0      0      0      0      0     0
#> fad432eeeff7f1d923e6d60aafe18dcd      0      0      0      0      0     0
#> a29db0e5f4a04c5c6937745ccccb167b      6      0      4     43      0    14
#> f179c10f8e86873c90ceb14b974a114e      0      0      0      0      0     2
#> a23e8588a2240d7af6f65519158c902c      0      0      0      0      0     0
#> 52034e0da6d5e7962369437956984e18      0      0      0      0      0     0
#> 522d20c27cf4237b60f36ffdca056c19      0      0      7      0      0     0
#> 957615aefb421e6ad4728e610b2bdfa6      0      0      0      0      0     0
#> 5ccc74e70b8a82a6eebea3385738ece3      0      0      0      0      0     0
#> e5924fd11a48419c6ec0e2a7cb4e43ee      0      0      0      0      0     0
#> eef7be4b2531cb7e909b58e91af5640d      3      0      0      4      0     6
#> 4416cbc2026260fa4e2ed0a8feff1443      0      0      5      0      0     0
#> 03f6091f5aacb69f06fc925acd1c9506      0      0      6      0      0     0
#> f655e5ddf5074d5e71bfa583926dec30      0      0     12      8      0     0
#> 9f733d39e588c7310499dd444cdd7812      0      0      0      0      0     0
#> 2b413de5e05a1385bab229f50879628f      0      0      0      0      0     0
#> 178b596f2fe5cb6cefc5513917988d84      0      0      0      7      0     0
#> aa94b014c726cf40983045bd606241c1      0      0      0      0      0     0
#> 55310b3aac4868c4e33b0a9492d7609a      0      0      0      0      0     0
#> c5d8d828df0e1baf5c2db05a6b8c4476      0      0      0      0      0     0
#> 6ceb26f60606233a7728d3457af84380      0      0      0      0      0     0
#> 1ff83868d1c9ecf52d9b776e0b3fd896      0      0      0      0      0     0
#> 3edfccfb9f81b340854e84ee5102a4c1      0      0      0      0      0     0
#> 963dc552c551872800f2fcbc5663e5fd      0      0      0     15      0    12
#> 3fda8c20e5da9fbc3fb7c5c864559019      0     10      0      0      0     0
#> 6fd88fa7e883bdead3c5ba89da853e8d      0      0     12      8      0     0
#> 35826fe826b7a9f202e8bef7b39d9a47      0     18      0     23      0     0
#> cc514f614ba146f489bf62b5a56074c5      0     21    172     33      0    20
#> 26260de3858081ce55e6f220b3ba25df      0     12     35      0      0     0
#> d6727a85782ffd2a069dccc23adaf521      0     61    155    132      0    68
#> eb1a74c2ca8831ef945775de0827ddda      0      0      0      0      0     0
#> 4ed3218615bbc1ed4da82d9806e4786b      0      0      0      0      0     0
#> 8a3367e36ed4e416d5a77a7946cf6fc5      0      0     34      8      0     0
#> 56f45eba0496ed29fb7171041a2a950d      6      0      0     17      0     7
#> 9d2782f19e60ec0cf72e1dc46b1d1766      0      0      0      0      0     0
#> 1559f8b64dc2bdc8cdc26ff8a5ffa4cb      0      0      0      0      0     0
#> 6d01f384a1cce140edea17dccd0f116d      0      0      0      0      0     6
#> 5950949b71082741302b712e97309b17      0      0      0      2      0     0
#> f6f238f95cc97995c9bf619b0c620b33      0      0      0      0      0     0
#> c7a112ffadf8009efd318625da44b18e      0      0      0      0      0     0
#> e79f4f30b0dacfcc50f3c1b12804ca3f      0      0      0      8      0     7
#> 297c825dc268f9a7b86b402edb6961b3      0      0      0      0      0     0
#> a1ed56ad4e86e0b4adc11baf2c4ab94b      0      0      0      0      0     0
#> 8afccb8f7f50c5603b471fb0c808c9e6      0      0      0      0      0     0
#> ce41c973bd4965a400989ec128890d55      0      5      5      0      0     0
#> 5d99ba689ad34e4393a905e8a725b72b      0      0     14      0      0     0
#> bfaa2d35e6dd13b9b32de811c6f92592      0      0      0      0      0     0
#> 0c544d5f3a49a010971a96d8c39e4905      0      0     11      2      0     0
#> bab3f08e1d75f42987012cae451639ff     57    238    165    191      0   111
#> fa97a81b8c33569e2aa1fc9a3137b435      0     10      0      0      0     9
#> a4d08886c3c6c62efd5678339a07828b      0      4      0      0      0     0
#> 11ed0972601847ec168473a6bc4eeba7      0      0      0      0      0     0
#> a6e4d88aa4b40402b1f94760d86b729c      0      0      0      0      0     0
#> 86deaff47d33946ede1b7c25c6bbe9f9      0      0      0      0      0     3
#> a1d3f1f8c151eaeb5cf12179861e0c87      0      0      0      0      0     0
#> 3980891c22ccbb387b65f9e437b2f94f      0      0      0      0      0     0
#> fc8ec321e0c2fc2a6de9abe98256b3f5     12      0      0     16      0     0
#> 84c11e6ffb29ebc2192f4b049b7dbddb      0      2      0      0      0     0
#> 82104585f0b2617c13eb69f0bbdf7d5a      0      0      0      0      0     0
#> 61cd4918428d3b36910b9cb2fc45946a      0      0      0      0      0     0
#> 0ab9a976935482bce5f8af4cc6685147      0      0      0      0      2     0
#> 87be3450216b45c53851fd173582a2d0      0      0      0      0      0     0
#> 93a51bb503fadf90121f518db980e158      0      0      0      0      0     0
#> b70cd843161996e075a853900e699d0f      0     16      0      0      0     0
#> 154b8be9eea158f93d27d314d8e0e2c8      0     20      4      0      0     0
#> 70691c2737cf7dc7c194de819f0e6cda      0      0      0      0      0     0
#> 2ae1d80ead92f5991f58da14d000acd6      0      0      0      0      0     7
#> 44fa299492152518fec34b9fa4554648      0      0      0      0      0    13
#> e15c8f9df87f363b0dbbf49de0bed249      0      0      0      0      0     0
#> 73689cad6eb4e86d9d0bb48c0da24463      0      0      0     15      0     0
#> 234e15e82d64f2e1d61802fa8ea1998a     38     12     43      0      0     0
#> 88c76ddb1b24acffe65d8128e2f48244      0      0      0     28      0     0
#> 5544139a16716cc6ef48b2c16c9613e9      0      0      0      0      0     0
#> 5b7fe68483e701059142cb531397b547      0      0      0      0      0     0
#> 0a0875a0acf4a28d1bc1c4d1d11aaa64      0      0      0      0      0     4
#> 6b47ec9d1c70f07d16834e9a4e54aa70      0      0      0      0      0     0
#> f8752b4824d239cfb19f3216cc34941b      0      0      0     43      0     0
#> ec8a0768cd2fb426d7e6bfc2c4f50e5b      0      0      0      4      0     0
#> 8a6b88f9515bf2dfac09e6a69045a4fe      0      0      0      0      0     0
#> f5151c67af19c8471743ef3e2914f6b1      0      0      0      0      0     4
#> c3842619fcbde52781c512c2561fca0b      0      0      0      0     10     0
#> b3b00faa9dee9fa6b3b33744f51da00e      0      0      0      0      0     0
#> c4f04e28c6267ac6c5b05d403e8f91c3      0      0      0      0      0     0
#> 7fc06f107a6ddf6486ef2d43fcf23626      0      0      0      0      0     0
#> 8d1dffcd377b338d2a0388637507b592      9      0     11      0      0     0
#> 15410d762889d51818986bfcdb5d8965      0      0      0      0      0     0
#> 48f68e44e258de63575b8bf6cf565d52      0     15     25      0      7    17
#> 4712befe3685661e0fc0968859c369d5      0      0      0      0      0     0
#> be9a3727a0422aea147100370d046fa8     13     27     72     81      0    40
#> 3c5942781761830f7b5406f551574424      0      0      0      0      0     7
#> 291d9f4d3d3c5cc3b0a1804770ed1e77      0      0      8      0      0     0
#> 784764c519a64dc9a147025a23222e93      0      0     10      2      0     5
#> b9e955ea5254dd68b54775aad2a67829      0      0      0      0      0     0
#> 8197f1594e8a17f3fc989fea64f54b8f      0      0      0      0      0     0
#> bf4a4f2383a16d2ecf27bf8e162b9ba4      0      2      0      0      0     3
#> 4e7d22d793d30bdedbd32a5b73b564b8      0     23     37      8      0     0
#> e88959b739c8cceb94b94b7887dd1df3      0      0     12      0      0     0
#> 5403445694339e2ae5ecac06d9f3ff47      0      0      0      0      0     0
#> ba3943fb075bbf63fb699638292826a1      0      0      0      0      0     0
#> 67b2d836db093cead5a8286142669613      0      0      0      0      0     0
#> 6c3a0ea164d3eb098dd2baac5519b38f      0      0      0      0      0    14
#> c643a776c36bcc00dfb5d790e1ea3dfe      0      0      3      0      0     0
#> d53cc81a9d57a64d194c876b76bd66e0      0      0      0     11      0     0
#> 7dc00d40971bb454b39fdd9b2dd9aa25      0      0      0      0      0     0
#> b4d0cec8543e940be1b5096370fdfb80      0      0      0      0      0     0
#> 6c54ceb94bc5e8fbd65224cffc60990a      0      0     41     16      5    12
#> b12ee9ca5a2320957213a114788cf2b7      0      0      0      0      0     0
#> 5e040c2397f959e7f0dfc6ec2855b812      0      0      0      0      0     0
#> 3105a594bba2604ae5ed19f5e0495d2e      0      0      8      0      0     0
#> d75c040bdaa1cca96cd0655c9af7ab28      0      0      0      0      0     0
#> 9b909edaebabc90b9c1a703626e6bd0c      0      0     15      0      0     0
#> a6593da81035cf121dbe46346a4ea960      0      0      0      0      0     0
#> 60e4c88400a750a8122b97476313a67e      0      0      0      0      0     0
#> 7854a78e3b881cee1c63ea0a50ba29ce      0      0      0      0      0     0
#> 1796734b56d45bdc997ccc09ab8731d1      0     14     18     22      0    17
#> 8243cd669ae7327f4326be78429ab05a      0     15     35      0      0     3
#> ebb93dc32c63402289114859d476cc5b      0      0      0      0      0     3
#> 8a01af5a531e851e2c50fc01131cda33      0      8     19     12      0     9
#> 03d2cc3f7b2bc352752d39e3738d0af9      0      0      0      2      0     0
#> ec378af6e32960b11b5f7391da4ae4fb      0      0      0      0      0     0
#> b407f2de868c064f65eb5ad8f82a8d0d      0      0      0      0      0     0
#> 73a95b6ba34fcd260989d848fe15d377      0      0      0      0      0     0
#> 3abb68202d15781590a2d1f9c5b37e15      0      0      0      0      7     0
#> 77923ef8e9ed01c9b6a8ac85eb028163      0      0      0      0      0     0
#> 9981b8fd20f400e052527f84f5db1b66      0      0     19     10      0     0
#> 597bdec49adf7675912f9fba4e262e70      0      0      0      0      0     0
#> 0fd3126118b13fb187bc55476a67c95a      0      3      0      3      0     0
#> 2a0c511ecc2513c42fe624e8e99d07da      0      0      0      0      0     0
#> d2912fea2dcf2e5b97e56c6a7d24295d      0      0      7      0      0     0
#> 2976af6aaf4261b6c274fcc3229a4d65      0      0      0      0      0     0
#> f097c45ce419da576aa111c98a4c7272      0      0      0      0      0     0
#> 3abc8260d34008f263e232ea25426c3f      0      0      2      0      0     0
#> d7e5b6432dc178aaeed9c81e0ec57afb      0      6      0      0      0     0
#> dfe71ce78ddfdbc7b453d17d628365f2      0     20     34     31      0    21
#> 289b18a86084273ce3b4038bf37e2840      0      0      6      0      0     0
#> 844a67032c52dc49019bf5dff248bf9c      0      0      0      0      0     0
#> 38cb3e3ce5ac5aed7dec50247ecfb267      3      0     30      9      0     0
#> 74e88e71e9929b1ddfcafccfab3b7bba      0      0      0      0      0     0
#> 66bc0737fc7fa1463ac7189ab7073998      0      5      0      0      0     0
#> b5a7aa588c9cae6ae8ca2f7d87ac53e8      0      0      0     31      0     0
#> db0540920457b21b170947aac9137771      0      0      0      0      0     0
#> ebbaf260510ab16ac6b8c2b9866d05b9      0      0      0      3      0     6
#> c1f461b67419d4a80a61db0867da9ed6      0     11      0     16      0     0
#> 0c6eba6452635eecfb865a1049299176      0      0     13      0      0     0
#> 7ecabb10f29a0c44eb53e9f9852fe22a      0      0      0      0      0     8
#> a33cb3037e7825920073b4dddaaa41e6      0      0      0      0      0     0
#> d625812c15c074c2e2b18099f97412c9      0      0      0      0      0     0
#> 17ef77f12e2a55717e910a5f6a625448      0     10      0     13      0     0
#> 7baf638fdd3f71166aafd22ecf410348      0      0      0      0      0     0
#> 9436cda211f60f722865cd23ef75d4ed      0      0      0     10      0     0
#> c8ab0b20e740eb6b35e68d58f69736e4     56     31     27     32      0     0
#> 5e43f2d8be71399fc8d6736754cb42ac      0      0      0      0      0     0
#> 19ad90ee1081a003406b02179ae63709      0      0      0      0      0     0
#> 0a25cdbd0660e68970c563c43a18f54a      0      0     19     38      0    24
#> 0aa6b210e902d8ab9d6f267407258ebb      0      0      0      0      0     0
#> a864c36d0d898cc5abd94d6d6f0bd5be      0      0      0     22      0     0
#> dea0dcb305ff0f27e3c36c7133cf933d      0      0      0      7      0     0
#> 463977baa69c7724fdc9407ad162714c      0      0      0      0      0     0
#> fcb72326c29342bcfb052953466eafc6     20     20     75     18      0     0
#> 5a554353f9ed728f0ab0631f29c8be2e      0      0      0      0      0     0
#> a9dd27b48e443415a743bf9902f515c4      0      0      0      0      0     0
#> 77a40bcf979c0099c52f19a3f1fe5217      0      0      0      0      0     0
#> b4007ac525a9b409d70426e24421cff0      0      0     13      0      0    24
#> 406ba18a175c4c79efb25f48a6beafe4      0      0      0      0      0     0
#> fe1c166755bc722cefd2e2a3b2b70672      0      0      0      0      0     0
#> 8afd274849ac07c0ca66f8f5f143df3c      0      0      0     17      0    20
#> 197301ed9d8728b5c061b99a17973161      0      0     18     23      0    35
#> 14779f9dc3505598ad3d799e4ee676c4      0      0      0      0      0     0
#> 4085c246640eee043d523ce2d799023a     30      0      0      0      0     0
#> bba09a403e67c1381de01c86511f1057      0      0      0      0      0     0
#> 392c59fee4e6b91d7654f2d3e3c65be8      0     27      0      0      0     0
#> d2b05ae73ad5ed296d7add9e2bb2cce4      0      0     12      7      0     0
#> 8306a9473dabfd5918ea869dbfcd1291      0     12     30      0      0     0
#> 06064c372bb624503c2905d13f7691f5     86     68     43     76      0    80
#> 3265692341e42dbd5cf4ed1938b906e3      0      0      0      0      0     0
#> da96a0797b9881a123795838b9f01140      0      0      0     24      0     0
#> 0bdddf10e577657becec88e42001fd16      0      0      0      0      0     0
#> 168999a950ebe17e6d2ea1f90e557d91      0      0      0      0      0     0
#> a266ffe9d38f626d106fa615da4618ef      0      0      0      0      0     0
#> 7706405d5bce7715f01d8b54188def42      0      0      0     20      0    14
#> 4b0cfc8899b466bd9992f7b22a5a06be      0      0      0     20      0     0
#> b81b2b4a7ca9eb2c23b170648dfaf22d      0     32     26     17      0     0
#> 75cf40872201f150cb5b191676f64b78      0      0     57     41      0     0
#> c8e5794c58de9326084a62d3282d257e      0      0      0      0      0     0
#> b0706a571f169c0bac7b4864052fcaa7      0      0      0      3      0     0
#> 0388f969e3527b2938acc676d6757bdb      0      0      0      0      0     0
#> 2dfe751390aedfbd9f62dbb77d1f10b4      0      0      0      0      0     0
#> 2548dc9697e1f9afe04a9012f40e2f0a      0      0      0      0      0    11
#> 2f019615d5141a4ab00b1aab98882a33      0      0      0     59      0    18
#> 629e4f0cf1826904ff74907b2507f58b      0      0      0      6      0     0
#> 774861f348cec8451102730c21a790b0      0      0      0     10      0    16
#> 11c0d6de28ab026a9ae17abe4efdd6bf      0      0      0      0      0     0
#> c11381dda0841c0d08cbfcef3a760aa3      0      0      0      0      0     0
#> 3c9619d4f9db4ac6ad7151a0c9d86ad3      0      0      0      0      0     0
#> 7e2c72204a8c45b1b0ac39e0e6047098      0      0      0     20      0    24
#> 8254e0f9f07709b62dd8483f3b25fad0      0      0     35     78     21     0
#> e78135f2886d427e678d12c446bf070e      0      0      0      0      0    11
#> abec938a8d53e7471f3f476d50694be9      0      0      0      0      0     0
#> 4963567d2979bb761291aa98e0b8595b      0     19     36     11      0     0
#> d4de9d1453102564f35262faa3952266      0      0      0      0      0     0
#> efa4a78c4b993e9d115841c8a02c795e      0      0      0      0      0     0
#> fcacd7424bf3543322ea3e9842d14e86      0      0     23      0      0     0
#> f941e15ea7f7560191d57abc4d2bedd8      7      0      0      0      0     0
#> ac294f5242de3c9b26e75502c48f6383      0      0      0      0      0     0
#> 7b6a8e3f8bfa1bcd1b5fc55db048d005      0      0     12      0      0     0
#> c9e776ab599fe75db286675b698c9067     18     56     15     25      0    26
#> ebe5cb945f6d6d849c07424591f80753      0      0      0      0      0     0
#> f116387f77802409103a559e4205df06      0      0      0      0      0     0
#> 93f9939b261cddc9a57071b42771c1d4      0      0      0     34      0    17
#> 971a5f48b3418449ca7cf3e190a5e831      0      0     10     30      0   108
#> 730de992ca06a95e7193a93146a97cf2     34      0     49      0      0    55
#> 8fa8b4ac69784c26901e036d13ab2332      0      0      0      0      0     0
#> c363c9efaa85023b488b732b97e47270      0      0      0      0      0     0
#> 2e2181cc893ca1e1e677e888bf8af441      0      0      0      0      0     0
#> 86fa662da6869634e7f90f85e9ee4933      0      0      0      0      0     0
#> 9a02177f2c235a6eb62984d5aab151ff      0      0      0      0      0     0
#> d6854be498e22b65da3c00572ee178f0      0     34      0      0      0     0
#> 015f364f150b230450da0a123e7bb8ab      0      0     11      0      0    20
#> af2c1b22758cd6653e8e841497d096e7      0      0      0      0      0    30
#> e7f2fab8c0e7c365768a43ef71406605      0      0     38      0      0     0
#> fb57cf192bce1876471f2f6181bf4b55      0      0      0      0      0     0
#> 6c740ce69c9f45bea684dc9e36054006      0      0      0      0      0     4
#> 3b27c7bf279535d8e39664d8dc1d0daf      0      0      0      0      0     0
#> ceea3c29ad012c2918fb6aa21a44e8cb      0      0      0      0      0     0
#> fa99ecdb27c1cc44fd866ee1da250c4c      0     15      0      0      0     0
#> 387cc83573db0e444e81c89a8a4f8f9c     17     75     97     63      0    46
#> c7a8d1af5f7f472671f700775de3e0a2      0      0     14      0      0     0
#> 98a6a84550e07ea252cfd8a0f51042d4      2      0     78      0      0    17
#> 08b809a68126c38db786cafe9710b013      0      0      0      4      0     0
#> 8db4cc0a5d2f03fce8efd2b8374ff359      0      0      0      0      0     0
#> 3e7c7d38423c93b43a1d9865cf0bc4a6      0      0     31     66      0    35
#> 6f02050760c030ed8cd1ed836697fa6b      0      0     18     15      0     0
#> 8fc1865cfcae32f72c1b1570299e2ef3      0      0     13     14      0    11
#> 38cc138fac95b006929361093d840a00      0     17      0      0      0     0
#> f5f8eb63e22362509e2c7c78271305f2      0      0     15      0      0     0
#> 0bb238aaee3ccb7dc8dfa86dbc1ba9f7      0      0      0      3      0     0
#> ca4546caf07edf5200cf79cea04aebcd      0      0      0      0      0     0
#> 893d12123eb7f623489a4f025a15ae7f      0      0      0      0      0     0
#> 192c5a817f48368339e519e0a0a907da      0      0      0      0      0     0
#> 3467185ca92b7bec2523214605adabf0      0      0      0      0      0     0
#> a2aae80d8fb922fe1074d243b59ce237     12     36     33     19      0     0
#> f11233f16d4177d38ab1c0bfc649e056      0      0      0      0      0     0
#> b219ca7793ab9f5e441f2a68bc76152c      0      0      0      0      0     0
#> e6e1848b51e705e55b3805703f088311      0      0      0      0      0     0
#> bd074b9bf365d11468bd76268d63b15b      0     24     37     22      0    13
#> f381b492d6651c52fd09ce9351e7eb2d      0      0      0      0      0     0
#> d9b31e1c30facc24fa778b1d65f5c457     40     12    144      0      0    14
#> 14a67c8c8cc31d203eac8109443fab5f      0      0      0      0      0     0
#> 026ad0ea53fb202ebb770acf65705d6e     12      0     20     17      0     0
#> 19ba9fea04969f975d6dc8d3da20d5e6      0      0      0      0      0     0
#> 0fc4426ed7dccc48ce27cfff8083433d    189     16    342    111    198    26
#> 45ad90e179842705fa984962f131fd61      0      0      0      0      0    15
#> d284afab226d22fd0596221d17db23d3      0     26     63     25      0     0
#> 38194655ba17dbfd7fcc1ed42a57c169     18     28     58     27      0    27
#> e30be5268bbebe597424bdf171cc1607      0      0     18      0      0     0
#> cb7ccdb4956134b452a3cb78d1751c57      0     10      0      0      0     0
#> b585436677ac88dd4ac6239655f2a268      0      0      0      0      0     0
#> f5d770eff5c9b9b0c859db57d175eb34      0      0      0      0      0    10
#> 1b45ffbc88c788dcee2239e68901cf05      0      0      0      0      0     0
#> 3fd951cf7cc3a0740a1f20f8a9c6d2c2      0      0      5     10      0     0
#> 2027d3e90ccc4003f1954322e99db6a9      0      0      0     18      0    20
#> 6157cb92771830955413f808e83fc420      0      0     23     25      0    17
#> 6349506cac32056a214b04239bdc3a33     30     49     99     82      0    58
#> 2254a8d37cf537b1e35097d70b6ac91b     83     21    156    346     15   177
#> 38f3bc8ef8836574de0fcbeafd348b65      0      0     52      0      0     0
#> fec23c70ae0c70ecceea1487cd978ee5      0      0      0     40      0    51
#> 503d123a6c78c6fb280a6472bd61b51e      0      0      0      0      0     0
#> 9f33de8570c639e727465f1ee4d750f6      0      0      2      0      0     0
#> d68e1899fa65515135d9bd8d356b7e3d      0      0      0      0      0     0
#> e15b88ccac7a4bef64215e915f6e0d89      0      0     13     30      0    13
#> 767042b4b19d8bfaaa426dad7dbce378      0      0      0      0      0     0
#> 4a3c4abee502d6185abd750522566bf3      0      0      5      0      0     0
#> c254aa5c00b7f27d756f7bbe96b3d4ad      0      0      0      0      0     0
#> fe8d8af791eb6e5d4335bb7604c74e9d      0      0      0      0      0     0
#> d81800d7e31b69877ce1eab9858b8f3a      0      0     30      0      0     0
#> c55c7cb80c5700ff9aded76192543b5a      0      0      5      0      0     0
#> 91dee7b3b90bbf2e53b2472b9bbde0e9      0      0      0     21      0     0
#> 6b28df61cb5dee74b76b1a13aaa0027d      0      0      0      0      0     0
#> 9a830fa26205ad91a02b5bcaebb0221a      0     13      6      0      0    25
#> a54b8142d2a09d058b81c377e8e286b8      0      0      0      0      0     0
#> 69ef3a03a4d733dbffd99fdaa1d2e62a      0      0      0      0      0     0
#> 083cbb131dd305a88b7d338bc371c35f      0      0      0     16      0     0
#> ef1a92caffb4d7cae286fdc77b0a93a6      0      0     26      0      0     0
#> 2cddd5dca166486f50a0828b28001e1e      0      0      0      0      0     0
#> bb899d1dcfd89cb200a76cb8951904f8      0      0      0      0      0    14
#> 7817f6851e9df03a3d992cf915310dba      0      0      0     11      0     0
#> 54eaf6571de743a86ccb531bee8aad7e      0      0      0      0      0     0
#>                                  3.13RI 5.8RO 3.13RO 5.8RI 2.11RI 2.3RI 6.16RO
#> 374a5eb6496d14e8e9d6f3d8c7e34d9d      0     0      0     0      0     0      0
#> 4c47d5cf81df1ac6e2a0560029a81578      0     0      0     0      0     0      0
#> 8fb6741a6685fad0374f85ad3e16156c      0     0      0     0      0     0      0
#> 310f3de009c95de5c938b8e811dfc96f      0     0      0     0      0     0      0
#> 7d3b39e8c0f0fdcb81b0c898a98c38f9      0     0      0     0      0     0      0
#> fc211549300b0954dbb4a4bf57e8a605      0     0      0    21      0     0      0
#> 53fad3d8c52f4b5a023ff92aee9c0f1a      0     0      0     0      0     0      0
#> 5e37342d650cde6cfc5bbe2726b12e3e      0     0      0     0      0     0      0
#> 56025f74cc8e54d4a054442d786ef096      0     0      0     0      0     0      0
#> d9e51d98f1c86a1c426552d252c944e4      0     0      0     0      0     0      0
#> bc905dcc609da83d3620f1986e63a9f7      0     0      0     0      0     0      0
#> 90009df763e86e6eca7baa54f01551e6      0     0      0     0      0     0      0
#> 6bbd12b392461134369d835715cd2765      0     0      0     0      0     0      0
#> 83c88a72c2dab171c38aadd98aa2d389      0     0      0     0      0     0     10
#> 2b93708de9f1c24ff6140814a599a3a8      0     0      0     0      0     0      0
#> 671a53a2acbf0c8c7e240e446daf0000      0     0      0     0      0     0      0
#> bfe35644221421300c0c17e00e465c05     39    18      0     0      0    45      0
#> d12bdf20f24b92aa8ec94c6c3f60497f      0     0      0     0      0     0      0
#> 3bdb72ce060bd89f107c5c2bc95399ed      0     0      0     0      0     0      0
#> e6a82be4b75aac8b928c1dfbe9ee9287      0     0      0     0      0     0     11
#> e3a63ba0a0ec40fbf270a57c180e0d3c      0     0      0     0      0     0      0
#> 2877cf981b128a18f788ba5437c1afa6      0     0      0     0      0     0      0
#> c015bc50696fe3e3d94ace01d9dcd691      0     0      0     0      0     0      0
#> eb2eeddf829080f8d9f25c6be5d507cd      0     0      0    22      0     0      0
#> b1fcda6d8df4c25c8e6d7124c076b1fe      0     0      0     0      0     0      0
#> de8bb5c39fb121fd048dfe0d639e6e8a      0     0     68     0      0     0      0
#> a234ac03223c275785dedf6df3b2aeff      0     0      0     0      0     0      0
#> 9ea1641eeed1f51fc228c653076d1d08      0    35      0     0      0     0      0
#> 694b926113eb8291fea166446890c8fe      0     0      0    39      0     0      0
#> b2af163540d17007833dc819704afb28      0     0      0    59      0     0      0
#> 6fc0ba1f8ff8f9259c9d492e59771755      0     0      0     0      0     0      0
#> 7b8c1dba059eac7ec7f3cefa203f6b2d      0     0      0     0      0     0      0
#> 0b06028193ef1ce19f5339a8ba7d4448      0     0      0     0      0     0      0
#> 56b4ec05ef3bbded515ec3dd2d573197      0     0      0     0      0     0      0
#> fc9c84f8767611251f0eb93af4ca10db      0     0      0     0      0     0      0
#> d3616bc0d8925a452274c484c73406c3      0     0      0     0      0     0      0
#> 7986b8b93f3afce920bd19a5d2229c09      0     0      0     0      0     0      0
#> bc78ee7297b6912452be0865d113f7ad      0     0      0     0      0     0      0
#> c472668b46b0d0f575108aa4170dd81f      0     0      0     0      0     0      0
#> 2658a78c5ff06760389ce2e3a548b4a2      0     0      0     0      0     0      0
#> 600c9e81fd1d1ddd5a37a5aa0f9f796d      0     0      0     0      0     0      0
#> f8d0bb97f093d8e8630ea62ac8a1bf13      0     0      0     0      0     0      0
#> 6c685b83d63d6d56a9dcaba34215ae6b     21     0      0     0      0     0      0
#> e0e7e6208280a5846a8b66f00667f8ca      0     0      0    37      0     0      0
#> 6d1aa1c1f93c423a856a5dd3b62e5c70      0     0      0     0      0    61      0
#> e553e766ea0e08f1f9fa2debe7946af6      0     0      0     0      0     0      0
#> 0f1ef43f6b23a567aba80599766d70bc      0     0     18     0      9     0      0
#> 3c4ccaa45e8d148c9a17b5d0f3720cd2    111     0      0   104    256    78     77
#> d9df6265d046a558d87d397b8a177d54      0     3      0     0      3     0      0
#> 4bbad82a39655a26a2de61b88c5af7fa      0     0      0     0      0     0      0
#> f5a6865451a3b138f4915c5ab445a0a7      0     0      0     0      0     0      9
#> 4a86fc550fa1ac37aa480d37edea875f      0     0      0     0      0     0      0
#> 02fc21e72f66adea757b355edea4a369      0     0      0     0      0     0      0
#> 25cf10f4f15dab00925cef52b9f08442      0     0      0     0      0     0      0
#> fa8c8b98d197c252368c160b43cefb98      0     0      0     0      0     0      0
#> c0b07d29aed1a875641aede53bc0a20c      0     0      0    11      0     0      0
#> 341433c1c8cc65ce47dcf32179d08fd8      0     0      0     0      0     0      0
#> 8bf175f329d16e4aa732cf2b32279df3      0     0     45     0      0     0      0
#> 062f7552b747dccab9586dff7275b5d9      0     0      0     0      0     0      0
#> 059e124ad6ba8f98c3d21db28a075ea8      0    99      0   148     37    81     85
#> 3d7a81faea3f783158352c9f739c8d17      0    68      0     0      0     0      0
#> 90be2fb2c82a98002d33ed7531e877e3      0     0      0     0      0     0      0
#> b0e263cdd1c189fc1757d8189db14be2     46    35      0    38      0    49     52
#> b962ce332fea338b05293d1a5a29cbca      0     0      0     0      0     0      0
#> 7a310a8f8e4bcbfb717263bd5376994d      0     0      0     0      0     0      0
#> 57136cbdae8fe066c2bb4a661003719b      0     0      0     0      0     0      0
#> 69fe46cab1cc9407470a2397bfc84103      0     0      0     0      0     0      0
#> b88efba97658a6936bb053985d026ace      0     0      0     0      0     0      0
#> 4d691fbbe2dbfb2815649fb72c757977      0     0      0     0     24     0      0
#> 49192cc8a222d97b965b839ff5d0f333      0     0      0     0      0     0      0
#> 859873fcd088c294a92a13340744eb6f      0     0      0     0      0     0      0
#> 0a8c44ca4d744e313e6dc5b239fca562      0     0      0    14      0     0      0
#> 09f94345824c82e8fe0d824cc4478fb5      0     0      0     0      0     0      0
#> 9dfbdb741fb0dc25d2206b8680bbf579      0     0      0     0      0     0      0
#> dc2bea463c720874cccd189ccdcf0e3e      0     0      0     0      0     0      0
#> 6273eb52991d528cf5281215f79b40a3      0     0      0     0      0     0      0
#> dc4bcbe74986e44ec3a040b866dfda9e      0     0      0     0      0     0      0
#> 96ec4781ee61f4083792f252f1b6cd3a      0     0      0    11      0     0      0
#> 9657d52acca51e701e7a3b7cedcfcc12      0     0      0     0      0     0      0
#> f65f183293d4710041057ab98f2c94ed      0     0      0     0      0     0      0
#> 8f6b3783b31a2dd8640de3df92dbb9c4      0     0      0     0      0     0      0
#> 0d90eb4842c9ad9126ed33c538c768ac      0     0      0     0      0     0      0
#> 8d3ee46aaaf7728594014be61c649347      0     6      0     0      0     0      0
#> e011a95fda1f1f2bafc4e5b83b52d989      0     0      0     0      0     0      0
#> ef8ea65ba4c2e2da3ba9556a87d2385d      0     0      0     0      0     0      0
#> dbbb92eb3a8a3feaf96092d46f2b6a45      0     0      0     0      0     0      0
#> 06cb03cddd8a12dbfaf27e8984b395c4     61    69     76   110    120    92     87
#> cdf134bb501c54d997945dabd44e35ff      0     0      0     0      0     0      0
#> f3c4561cc01a3b45ee3f134d709c6b94      0     0      0     0      0     0      0
#> b121442fd1eefe78ce4f4602aac1842e      0     0      0     0      0     0      0
#> ba8707aa88a2a62939d9338189b2a733      9     0      0     0      0     0      0
#> ecf7289628d49f9244469337ee1622cb      0    57      0    49      0    36     44
#> 9f06b400b93f859b11874b0bd7f09295      0     0      0     0      0     0      0
#> f775be5cf7c0dd2177006a44913503d5      0     0      0     0      0     0      0
#> e87ceb10f2becb2a6c4cf0692b1923ee      0    22      0    18     11    25      0
#> 036574bb8e05f3b9a11ab7a4090d6950      0     0     48     0      0     0      0
#> 0267ba9db97e00ae01fc2cccbb6263e0      0     0      0     0     78     0      0
#> be4f701c1dbfbda16c5150af29430141      0     0      0    86     61     0     39
#> 655a66e6566798ade91f73b2942c0b63      0     0      0     0      0     0      0
#> cd79f2696509d0ad906f5addba264079      0     0      0     2      0     0      0
#> 07ea3284bda2bf7daf1251ad2744e4c4      0     0      0     0      0     0      0
#> e3870da9de4a2b18a7ae9e9cae208e3e      0     0      0     0      0     0      0
#> e00bd59c1ea0ca9b1cc5f976801d4e1b      0     0      0     0      0     0      0
#> b39b381d0f533c9923c3b898ffab5614      0     0      0     0      0     0      0
#> e2959f0552fa6df5707d9ec8b57816d6      0     0      0    23      0     0      0
#> 53157acd1642f973842d01e1cd287421      0     0      0     0      0     0      0
#> 5991790d1c0aa8b73b2d7f906449a1f9      0     0      0     0      0     0      0
#> adee3552c211d9e33aaa28e82532bdea      0     0      0     0      0     0      0
#> 0089040d041888e3ae4aa4e0010a0e36      0     0      0     0      0     0      0
#> 8e876e13d051992c4c8616cf2db7d73a      0     0      0     0      0     0      0
#> 179abbc5c1f1549ffd0d7ae6963e3783      0     0      0     0      0     0      0
#> 340d095f9a9fb445582b60e1282d985d      0     0      0     0      0     0      0
#> 3d7e5432d8eaa3fa56ab7138605ea43d      0     0      0     0      0     0      0
#> d4b233740db907d69035915e8657fa8d      0     0      0     0      0     0      0
#> b22119df3b1de76745abaac585b646f8      0     0      0     0      0     0      0
#> 64761a24fc66e7c93145217dd153b1ee      0     0      0     0      0     0      0
#> eb6d7f9e0601207fdcfc4a88947d513c      0     0      0     0      0     0      0
#> 0d20271e3e2dbc68e1de868599a9a1ae      0     0      5     0      0     0      0
#> 369ea22772cb774950c830d0af7e8235      0     0      0     0      0     0      0
#> cae63066cfe088a49fcac50f119c9e2d     23     0      0     0      0     0      0
#> d9c57ec20ed90e2cebf6cacc4f9a5511      0     0      0     0      0     0      0
#> f5828cc3078b4c219a44c67e96e5d455      2     0      0     0      0     0      0
#> 73847ac2778d4cba27fe9758ca129d6c      0     0      0     0      0     0      0
#> d78003f13712da48d0bb36fd45116b91      0     0     11     0      0     0      0
#> eb40853d986ad72c82658d939f291ff2      0    20      0     0      0     0      0
#> 5d675a3518222fc99c8cba34feaeac81     60    39      0     0     58     0      0
#> 3ba390ae11c4e76375af082eefe7e41f      0     0      0     0      0     0      0
#> 3aa66daa4e560cf62c1833697bcd3837      0     0      0     0      0     0      0
#> 6f4388510be5ea70356f2257d8c57ddf      0     0      0     0      0     0      0
#> b33390d034148edf07d8a742778ddcbd      0     0      0     0      0     0      0
#> 09225ec7044394bbd3b6811a6799e481      0     0      0     0      0     0      0
#> fe6fbf9b86e0f487911a856009fa4a1e     45     0      0     0      0     0      0
#> 0aee6329e0b5b81449980151bf40205a      0     0      0     0      0     0      0
#> 0ef66a78236f90916afc2305f3ff6c80      0     0      0     0      0     0      0
#> fb1eb41359005c2e3e308fb668f07220      0     0      0     0      0     0      0
#> cb59c57800227006aead644df42b8986      0     0      0     0     21     0      0
#> 138da9dc244c3292536cacb8df79467b      0     0      0     0      0     0      0
#> 72f92361ec0696cb8df985b93d277a25      0     0      0    23     25     0     21
#> 4d5ab121763879fd59e8f00d35473540      0     0      0     0      0     0      0
#> 20032b2b586fee22805eccebe3059170      0    51      0     0      0     0      0
#> 710324ccc8cafa6fc9b38baf8299d8aa      0     0      0     0      0     0      0
#> 523e6b60340f4a215877b1c336822e50      0     0      0     0      0     0      0
#> 4501f98c83484fd902be11383c8d032b      0     0      0     0      0     0      0
#> d4cbf23757aad57c91b058cbb05566e4      0     0      0     0      0     0      0
#> de7964d5c22112693ed5796bb9f04b0b      0     0      0     0      0     0      0
#> f4a283be42989b920b83ea73ec72b25b      0     0      0     0      0     0      0
#> 2f3111649d0850ab86799f8bbadc28f8      0     0      0     0      0     0      0
#> b8acdd45cf33ab75fe4c0851a1787dc5      0     0      0     0      0     0      0
#> 5bb7fec39744a6d4e4d61596dabc6886      0     0      0     0      0     0      0
#> e8eddfb37c4215c86813163a27f202a5      0     0      0     0      0     0      0
#> 826e188c911c0846fd927300127c80c0      0     0      0     0      0     0     11
#> ff1fe599c44074ee6bdf919a119df93c      0     0      0     0      0     0      0
#> aff00055d3876f624d67fb80bdb59934      0     0      0     0      0     0      0
#> 1beded194f3ae4b57715b45e50adef46      0     0      0     0      0     0      0
#> 908d91feba7c1f1d8cc8a7145173d942      0     0      0     0      0     0      0
#> e7095743695f6d4337d01ceb9e29e224      0     0      0     0      0     0      0
#> 4cd420962855afdf0bc69bb268b846c7     41    96      0     0      0    59     44
#> 430ba5bea0fc08f65b69497bfbeb5c4b      0     0      0     0      0     0      0
#> e12cb2f0d2ce4401413176007dcd4158     20    22      0     0      0    26      0
#> ec615213f6c7febb65c9f2e45e438111      0     0      0     5      0     0      0
#> 412dd52cd54d0244ba8fc42ad140acda      0     0      0     0      0     0      0
#> 89192d9f779341d313f52aad0dc6a06f      0     0      0     0      0     0      0
#> ed543f92c61a06fc4e703655ce910a2f      0    15      0     0      0     0      0
#> fea1a70b0ef0f72e3bf27d8b295331a8      0     0      0    44      0     0     17
#> 6afe805b36f5aa29e7f2096aac42ace1      0     0      0     0      0     0      0
#> d9cc5a46831f24ec9192ae82dbc38b0a      0     0     22     0      0     0      0
#> 0d3eb38eb502c4ad13bd5a7060295411      0     0      0     0      0     0      0
#> 69a0ff5538a5490e291c75ce70325066      0     0      0     0      0     0      0
#> f91d545c815effa3c5afa83a8ad3397d      0     0      0     0      0     0      0
#> e6b032dea8e13c34af958fc48af09b58      0     0      0     0      0     0      0
#> 84a1b89045ab6285328f366224f3ec87      0     0      0     0      0     0      0
#> 4a6f0076e4ac71474b3c364294ad3b33      0     0      0     0      0     0      0
#> 7bdd39e4e3ccc7000862b487510e67d6      0     0      0     0      0     0      0
#> 998aedf12cdf054f5e504a2e21d20498      0     0      0    13      0     0      0
#> 567c2640108006c8eca081ad4180e346      0     0      0     0      0     0      0
#> e50a6625b9176b6856764b213013150f      0     0      0    62     79    67      0
#> 8d5a5ab7c784b0325a4a86a16fbba80e      0     0      0     0      0     0      0
#> f58ae28a309fa3d5898d2eb46d9eb9e0      0     0      0     0      0     0      0
#> 7ac65a8fd21111b68e7c2d6205112a83      0     0      0     0      0     0      0
#> d7add64c8700f5a8cadf40c8c7f5ff67      0     0      0     0      0     0      0
#> 6bad9d981d3cf0f90262547180dc6371      0     0      0     0      0     0      0
#> 149c0e63c83113e1246c0e041de67d9c      0     0      0     0      2     0      0
#> 6f196fae54ea90b5c11627c503c98940      0     0      0     0      0     0      0
#> b67c995d842e5b4b92fa0ce052603ef1      0     0      0     0      0     0      0
#> f44ffae67c7297e2c696b75fbd4e0f4c      0     0      0     0      0     0      0
#> 33997912cd952feaec18fbd11b31acac      0    24      0     0      0    34      0
#> 0a035423c4bf71a17566549f9b69448d      0     0      0     0      0     0      0
#> f65a6f1f1e5a4bc5c733ba364c0096d3      0    38      0    25      0     0      0
#> f15389a9f245465c33db3de72aeb17ba      0     0      0     0      0     0      0
#> cd2b676d3c0785a60359c0d416701ee4      0     0      0     0      0     0      0
#> 70ba2a378ebcd9a8f57b024a561cc6b1      0     0      0     0      0     0      0
#> 0273707384afdd41e9d022b539094211      0     0      0     0      0     0     45
#> f9ac98e876ecdb504202388d7ccdda5f      0     0      0     0      0     0      0
#> bdd9b87b782b8f935ecac62a9824e6a1      0     0      0     0      0     0      0
#> 1f12b97495668f1fc378a349546d1806      0     0      0     0      0     0      0
#> 83f1b4d776f6e08a6ba6f86634a1e5a9      0     0      0     0      0     0      0
#> 49f614420449f6c7ecc8c9f4e137ad05      0     0      0     0      0     0      0
#> aba5aa00a1b262f27517eca64ee3333f      0     0      0     0      0     0      0
#> 8971980484eeabeef00145b93b699e30      0     0      0     0      0     0      0
#> c1347ae54b0fc1dd385a1ccf3b8061ce      0     0      0    42      0     0      0
#> e9fe19783ebcc7e43a3a07ca8225db62      0     0      0     0      0     0      0
#> 7c523a5f0a7c736efd82afef84137901      0     0      0     0      0     0      0
#> e13c07ec336c89db5183cc6e9e57ee6e      0     0      0     0      0     0      0
#> a017ce431ccadea609b501d46e64e55b      0     0      0     0      0     0      0
#> fadf38d03d4b944f5782a759591ec3ae      0     0      0     0      0     0      0
#> 41bfda1c7f32205246e73a1e3e3ec4f4      0     0      0     0      0     0      0
#> 9e4f3e697b70758a7a996f3dd4533650      0     0      0     0      0     0      0
#> 39c4f6caa7510ed108f3ebcf72f79edf      0     0      0     0      0     0      0
#> 4aa083fe7c1b0639dbc48ec649e68a31      0     0      0     0      0     0      0
#> fd2ac2e14b72b0e22a052e3ab4580471      0     0      0     0      0     0      0
#> 4ee10294799bcfeab0c8ee567d5d9608      0     0      0     0      0     0      0
#> 46869f12c31de7a05781e930c25a1ccd      0     0      0     0      0     0      0
#> 5e2548d6a970645d0ca943a72d207489      0     0      0     0      0     0      0
#> 0538e1beed69f098c3f602130ee11a1c      0     0      0     0      0     0      0
#> 991c563278c4e8c91ffebce3560fe51b      0     0      0     0      0     0      0
#> ce05f58164e61517198f7ab0e3304cb5      7     0      0    12      0    10      0
#> fdbd0fed21bae91c05b2da268a025a89      0     0      0     0      0     0      0
#> 7412d11595d77f26792071cd28a26760      0     0      0     0      0     0      0
#> 4b087e13957eb38128f63e1c99df23c3      0     0      0     0      0     0      0
#> 310a6f0e6841db2ce75005cac942a972      0     0      0     0      0     0      0
#> 3df0fada803ac582802b47d4f4efe0c1      0     0      0     0      0     0      0
#> 3cd4e34c02ffdd1446df65bb54278c57      0     0      0     0      0     0      0
#> b14a1dea93554257791f16c2e94906b0      0     0      0     0      0     0      0
#> 902a4435fc3af101e8a9b5f51550fa4d      0     0      0     0      0     0      0
#> 1995185b31349529f58964c2e5a3becd      0     0      0     0      0     0      0
#> 69b247c565afdeacee0240482b1ce536      0     0      0     0      0     0      0
#> 95ff1b9528acca018e9c4697a8625504      0     0      0     0      0     0      0
#> 132d26d3ce6ac2ab25ccb4dd171879d7      0     0      0     0      0     0      0
#> bf1a6e4f94ee8ba60cf79089139641b6      0     0      0     0      0     0      0
#> 704189aa09f64b992b3a02ee001dc706      0     0      0     0      0     0      0
#> dea5ea27839ee67a31d18ba0bf0f7f65      0     0      0     0      0     0      0
#> 3840a4fff0250232f8570cae7c3514ba      0     0      0     0      0     0      0
#> 31ecabc032997aaf41dd0d8def3ec93f      0     0      0     0      0     0      0
#> 024a0e0348ed9abe9ae1c48b1b7ba09e      0     0      0     0      0     0      0
#> fc1d9940419420113ce3fbfacc8d703a     76   112      0    30      0    34      0
#> 0770fbf28cd307a08bf557af7ffecfd3      0     0      0     0      0     0      0
#> 883cbe8e4d1d98f468355e8ea974d222      0     0      0     0      0     0      0
#> cb6d6561383fb68f43f92a88a607749d      0     0      0     0      0     0      0
#> e07e0f1a9fbc696a9fe4f9b7739ae45b      0     0     82     0      0     0      0
#> 73d8d8f999acfa13cba2277f7526897e      0     0      0     0      0     0      0
#> 5218274c354d41965abbf4b0751d6f3a      0     0      0     0      0     0      0
#> 5c31570f706752fa4f7851aca5b5c291      0     0      0     0      0     0      0
#> 6a5e71c8b86a52ab1dd5e61164103940     12     8      0    42      0    20     18
#> 25641b0b803c7d8493ebc02cb104bf3d      0     0      0     0      0     0      0
#> 8844d70dd1fdb96f2a7faa65bccf78c6      0     0      0     0      0     0      0
#> 371b855afa8bed030e86ec78241b01a3      0     0      0     0      3     0      5
#> d0c180a47378d414a1d088ac14db236c      0     0      0     0      0     0      0
#> a072a588b3940fab3b227e4901c40f7f      0     0      0     0      0     0      0
#> 977a3523051e6b23da6221dbde795b34      0     0      0     0      0     0      0
#> 1c3adc87d2e93b08e576953508eb680d      0     0      0     2      0     0      0
#> 74eb0ddf379b00c5032b18fb9fbf1f32      0    10      5    19     37    11     23
#> c335b34a15db7d532f511e6c952abe62      0    10      0     0      0     6      0
#> 0c78285397fb853a0e39069213510c54      0     0      0     0      0     0      0
#> 779dec6de9c6a8b7df28fc8ec3c75bee      0     0      0     0      0     0      0
#> 5526624b473bb1f0541dc6fb95c98678      0     0      0     0      0     0      0
#> 649e7fa7b9ce5c5472d34675ed643639      4     0      0     7      0     0      0
#> a9c5addeea542953d1bae4d30d599b02     11    16      0     7      7     0      0
#> 68291fb3b2558d438da4810961c8d53a      0     0      0     0      0     0      0
#> 9c5d92bfefc7a190b93d5b9a6a77697b      9    12      8    27     26     0     10
#> 3ddab42d701754625cab0c236c5b6e6d      0     0      0     3      7     0      2
#> 287e840f9ecaef80d7a4762c8d6bf401      0     0      0     0      0     8      0
#> 885f411fa54c5576148678f665bb053f      0     0      0     0      0     0      0
#> b5709a57df9cbc8a741bb78a36f29c06     38    54     31   101    104    50     60
#> 26698cc60893017ef0b1ab173688d3de      0    18      0     0      0     0      0
#> db64f0028b6518839ab7129c4b7ef72d     56    40     66   169    135    63    105
#> f46fb43262da299ec7f0084298cafa50      0     0      0     0      0     0      0
#> c45e7087de669be270881fc8ace0fad8     49     0     36    60     63    56     39
#> 854a20d5019388e81b1e78419eb017fd      0     0      0    81     15     0     28
#> 9a2325ea015a8c5944047e3349cc1565     15     0      0    56     34     0     19
#> c5305ae3c6adc90814d70ab48f4f071a      0     0      0     0      0     0      0
#> e230123ec3ce8654942ee1a000a80010      0    24     17    61     61    35     31
#> 14e3c60cb885509667670151bc4a1df0      0     0      0     0      0     0      0
#> 7fd7b7bece92bbbfcceb2dcf1e835a8e      0     0      0     0      0     0      0
#> c098d1fb769dd975c2332dfb70da489e      0     0      0     0      0     0      0
#> 95dcdf1eb06b252e67dcdbdaffc505e8      0     0      0     0      0     0      0
#> c5399d0f9eda8f0815fd28a257952496      0     0      0     0      0     0      0
#> c7388cc314e992819b8f97ea35280982      0     0      0     0      0     0      0
#> ae8facc45d372f7108bd9499f74d1bfc      0     0      0     0     76     0      0
#> 16caf6f5bc653fb6ba668eb821076a21      0     0      0     0      0     0      0
#> e7e283b72ee15bcf0877c5eeb6137ebc      0     0      0     0      0     0      0
#> 1b3a779b384c73b33207dccd14b14fe6      0     0     62     0      0     0      0
#> c8d1afdf7b395ffc22193c5dc2fd5dc3      0    16      0    11      0     7      0
#> 0e3ffc2a68378086e2480831acd1a025      0     0      0     0      0     0      0
#> 0ddf7339b608a340252b33f63fcc0019      0    11      0    19     14    11     17
#> 5fd4b2f7dee3fb0bde6f856ae266dbcb      0     0      0     0      0     0      0
#> e2f4786bc66410c07bd8b96b5d0e2d92      0     0      0    17      0    25      0
#> c30aae8f5216c9536492d0beae74950f      0    30      0    92    171   123      0
#> 9087ecee0806a80861788a17dd7e6809      0     0      0     0      0     0     19
#> c29bc08418318360592b2b2717376a0a      0     0      0     0      0     0      0
#> f791c7db1a78d35e84803b5e03c01965      0     0      0    17      0     0      0
#> 44b44fc4055a58c01f6d1e9de5e0d1d9      0   132      0    40      0     0     58
#> 342df9401103d274851d3bd80704508d      0     0      0     0      0     0      0
#> 0590e6672bb2c55c0ff54999a2b6dcfe      0     0      0     0      0     0      0
#> d4e9147608bf883fdbbf9d23b17028d1     20     0      0     0      0     0      0
#> db24794eb8667aed224401386223f671      0    11     16    14     24    12     57
#> fe41e6fd9770030a21418b288930710f      0     0      0     0      0     4      0
#> 699f18f3b155bff97de941f2cfbbbadc      0     0      0     0     12     0     29
#> 13870cf319997750e3b96bbb0b2f1aea      0     0      0     0      0     0      0
#> 7ca4fb1b1740406241acb862445cf311      0    14     14     0     39     0      0
#> 9c0569fb95451d009fd3754303eaed12      0     0      0     0      0     0      0
#> 1ac26a3c12288d8fbadd6aa61366b1b7      0     0      0     0      0     0      0
#> 89f3ba99552605501b8acb816d288afe      0     0      0     0      0     0      0
#> b4051b5227746f04c1feb9d8979299ff      0     0      0     0      0     0      0
#> 087d358aaf77f81aaf50a623828865bb      0     0      0     0      0     0      0
#> e451c4b3ff1c3dabeee296529b3a7a07      0     0      0     0      0     0      0
#> 012871784f77c9dca7d446ddef2048d1      0     0      0     0      0     0      0
#> 95c6ce2788783bf2f83343deba23b8bd      0     0      0     0      0     0      0
#> 213227f244d97553a9a14c4e3fd04c98      0     0      0     0      0     0      0
#> c0e5e55c5cfbd576e7c085f1e4b69638     11     0      0     0      0     0      0
#> 51e98451d7e3f46127481ec5490f1c91      0     0      0     0      0     0      0
#> c75285a6c61d8fce20d823769c25c12a      0     0      0     0      0     0      0
#> 177966f7619eb24fdbce2266f5bed39a      0     0      0     0      0     6      0
#> b89347f0534679f890628364d2c7b4fe      0     0      0     0      0     0     22
#> 8dff1cb12c0beb5bc7b2546480a6a6c7      0     0      0     0      0     0      0
#> 84f2777cfcb0fbed6281917fcc765022      0     0      0     0      0    15      0
#> a5e4d106a2ab947984fdef187996a0b9      0     0      0     0      0     0      0
#> b6cbf6077b64acbb718971a1ffc97a27      0     0      0     0      0     0      0
#> 27e29495fda424728fd6a80b83a16772      0     0      0    16      0    13     12
#> 63c44ab3725dd0b803d43a29db9fd7fd      0    20      0     0      0     0      0
#> 9c931ab5bf18ef9b10996d59d292e6f5     15     0      0    48      0    61     25
#> 12fd8efb82c3d4909543f6db9853134e      6     0      0     0      0     0      0
#> a8a60c20f7ab7d0f2dadd9f8858c19be      0     0      0     0      0     0      0
#> d57e09eb1645f1abba6074d3bd7d4c52      0     0      0     0      0     0      0
#> 7e0f2b982ca2d35b6ac2cf251d78d3c4      0     0      0     0      0     0      0
#> 3e2e8e938a68a0318bcf7145ad5ed4a9      0     0      0    12      6     0      0
#> f2117b49bf57d9451bbef6de5b762483      0     0      0     0      0     0      0
#> be24dfe003725a0b0f28f4136106faf8      0     0      0     0      0     0      0
#> e63301d7d0af55d5c43db051dc51522f      0     0      0     0      0    25      0
#> a6e2e5501507562f2174f65e606b2d43     10    18      0     0     15     0      0
#> a06aa93cb50b9a39941a11c7242dde96      0     3      0     0      0     0      0
#> 6a0a8abb36aa3bfdee14edbe97e9d781      0     0      0     8      0     0      0
#> 4d3def5ffad0c392381591440bfc3668      0     0      0     0      0     0      0
#> ce437b81e7cf5ce986f7240c26464b0a      0     0      3     0      0     0      0
#> e57788f9d2fae76ed9ca566ad4d475a8      0     0      0     0      0     0      0
#> 8c525d0f63137bd068b75c57f5d04c09      0     0      0     0      0     0      0
#> 6d670aa40dafa152c2d59cf2be9bef74      0     0      0     0      0     0      0
#> e08f8e9e9ffb10a702f2ef617f70eecf      0     0      0     0      0    12      0
#> 7def1c807ca0767ccdfe5a6001dc51f1      0     0      0     0      0     0      0
#> 3b9a17275147adf1dc8c42eb00a06a36      0     0      0     0      0     0      0
#> 24be0b7e4c85ac53b09fddbca916ab96      0     0      2     0      0     0      0
#> f521a47d3bb4acfec6b261f880b4eea8      0     0      0     0      0     0      0
#> 91c24c9b6ffbedcc257d164176060fe0      0     0      0     3      0     0      0
#> 9b45f210647c5fe71b005eb9e25299a5      0     0      0     0      0     0      0
#> 974077d25fa6b32a7af22670b142a8f1      0     0      0     0      0     0      0
#> 203aa2535e54a11e4c48af4ab946686b      0     0      0     0      0     0      0
#> 9064afe569448558eb08e5bf5e458bce      0     0      0     0      0     0      0
#> fad432eeeff7f1d923e6d60aafe18dcd      0     0      0     2      0     0      0
#> a29db0e5f4a04c5c6937745ccccb167b      9     0      0    26      0    22      6
#> f179c10f8e86873c90ceb14b974a114e      0     0      0     0      0     0      0
#> a23e8588a2240d7af6f65519158c902c      0     0      0     0      0     0      0
#> 52034e0da6d5e7962369437956984e18      0     0      0     0      0     0      0
#> 522d20c27cf4237b60f36ffdca056c19      0     0      0     0      0     0      0
#> 957615aefb421e6ad4728e610b2bdfa6      0     0      0     0      0     3      0
#> 5ccc74e70b8a82a6eebea3385738ece3      0     0      0     3      0     0      0
#> e5924fd11a48419c6ec0e2a7cb4e43ee      0     0      0     0      0     0      0
#> eef7be4b2531cb7e909b58e91af5640d      0     0      0     0      0     0      0
#> 4416cbc2026260fa4e2ed0a8feff1443      0     0      0     0      0     0      0
#> 03f6091f5aacb69f06fc925acd1c9506      0     0      0     5      0     5      0
#> f655e5ddf5074d5e71bfa583926dec30      0     0      0     0      0     0      9
#> 9f733d39e588c7310499dd444cdd7812      3     0      0     0      0     0      0
#> 2b413de5e05a1385bab229f50879628f      0     0      0     0      0     0      0
#> 178b596f2fe5cb6cefc5513917988d84      0     0      0     0      0     0      0
#> aa94b014c726cf40983045bd606241c1      0     0      0     2      0     0      0
#> 55310b3aac4868c4e33b0a9492d7609a      0     0      0     0      0     0      0
#> c5d8d828df0e1baf5c2db05a6b8c4476      0     0      0     0      0     0      0
#> 6ceb26f60606233a7728d3457af84380      0     0      0     0      0     0      0
#> 1ff83868d1c9ecf52d9b776e0b3fd896      0     0      0     0      0     0      0
#> 3edfccfb9f81b340854e84ee5102a4c1      0     0      0     0      0     0      0
#> 963dc552c551872800f2fcbc5663e5fd      0     3      0     0      0     0      0
#> 3fda8c20e5da9fbc3fb7c5c864559019      0     0      0     0     10     5      0
#> 6fd88fa7e883bdead3c5ba89da853e8d      0     0      0     5      0     0      0
#> 35826fe826b7a9f202e8bef7b39d9a47      0     0      0    10      0     0      0
#> cc514f614ba146f489bf62b5a56074c5     42    13     33     0     20    12     60
#> 26260de3858081ce55e6f220b3ba25df      0     0      9     0     14     0      0
#> d6727a85782ffd2a069dccc23adaf521      0    35      0    70     70    42     59
#> eb1a74c2ca8831ef945775de0827ddda      0     0      0     0      3     0      0
#> 4ed3218615bbc1ed4da82d9806e4786b      0     0      0     0      0     0      0
#> 8a3367e36ed4e416d5a77a7946cf6fc5      0     0      0     0      7     0      7
#> 56f45eba0496ed29fb7171041a2a950d      2     3      0     0      0     7      0
#> 9d2782f19e60ec0cf72e1dc46b1d1766      0     0      0     0      0     0      0
#> 1559f8b64dc2bdc8cdc26ff8a5ffa4cb      0     0      0     9      0     0      0
#> 6d01f384a1cce140edea17dccd0f116d      0     0      0     0      0     0      0
#> 5950949b71082741302b712e97309b17      0     0      0     0      0     0      0
#> f6f238f95cc97995c9bf619b0c620b33      0     0      0     0      0     0      0
#> c7a112ffadf8009efd318625da44b18e      0     0      0     0      0     0      0
#> e79f4f30b0dacfcc50f3c1b12804ca3f      0    15      0     0      0     0      0
#> 297c825dc268f9a7b86b402edb6961b3      7     0      0     0      0     7     10
#> a1ed56ad4e86e0b4adc11baf2c4ab94b      0     0      0     0      0     0      0
#> 8afccb8f7f50c5603b471fb0c808c9e6      0     0      0     0      0     0      0
#> ce41c973bd4965a400989ec128890d55      0     0      8    10      0     0      0
#> 5d99ba689ad34e4393a905e8a725b72b      0     0      0     6      0     0      0
#> bfaa2d35e6dd13b9b32de811c6f92592      0     0      0     0      0     0      0
#> 0c544d5f3a49a010971a96d8c39e4905      0     0      0     0      0     0      0
#> bab3f08e1d75f42987012cae451639ff     18    44     14   158    122    78     41
#> fa97a81b8c33569e2aa1fc9a3137b435      0     0      0     0     12     0      0
#> a4d08886c3c6c62efd5678339a07828b      0     0      0     0     21     0      0
#> 11ed0972601847ec168473a6bc4eeba7      0     0      0     0     14     0      0
#> a6e4d88aa4b40402b1f94760d86b729c      0     0      0     0      0     0      0
#> 86deaff47d33946ede1b7c25c6bbe9f9      0     0      0     0      0     0      0
#> a1d3f1f8c151eaeb5cf12179861e0c87      0     0      0     0      0     0      0
#> 3980891c22ccbb387b65f9e437b2f94f      0     0      0     0      0     0      0
#> fc8ec321e0c2fc2a6de9abe98256b3f5      0     0      0     0      0    40     51
#> 84c11e6ffb29ebc2192f4b049b7dbddb      0     0      0     0      0     0      0
#> 82104585f0b2617c13eb69f0bbdf7d5a      0     0      0     0      2     0      0
#> 61cd4918428d3b36910b9cb2fc45946a      0     0      0     0      0     0      0
#> 0ab9a976935482bce5f8af4cc6685147      0     2      0     0      0     0      0
#> 87be3450216b45c53851fd173582a2d0      0     0      0     0     13     0      0
#> 93a51bb503fadf90121f518db980e158      0     0      0     0      7     0      0
#> b70cd843161996e075a853900e699d0f      0     0      0     0     56     0      0
#> 154b8be9eea158f93d27d314d8e0e2c8      0     0      0     8     31     3      0
#> 70691c2737cf7dc7c194de819f0e6cda      0     0      0    13      0     0     17
#> 2ae1d80ead92f5991f58da14d000acd6      0     0      0     0      0     0      0
#> 44fa299492152518fec34b9fa4554648      0     0      0     0      0     0      0
#> e15c8f9df87f363b0dbbf49de0bed249      0     0      0    15      0     0      0
#> 73689cad6eb4e86d9d0bb48c0da24463      0     0      0     0      0     0      0
#> 234e15e82d64f2e1d61802fa8ea1998a      0     0      0    69      0    38     23
#> 88c76ddb1b24acffe65d8128e2f48244      0     0      0     0      0     0      0
#> 5544139a16716cc6ef48b2c16c9613e9      0     0      0     0      0     0      0
#> 5b7fe68483e701059142cb531397b547      0     0      0     0      0     0      0
#> 0a0875a0acf4a28d1bc1c4d1d11aaa64      0     0      0     0      0     0      0
#> 6b47ec9d1c70f07d16834e9a4e54aa70      0     0      0     0      0     0      0
#> f8752b4824d239cfb19f3216cc34941b      0     0      0     0      0     0      0
#> ec8a0768cd2fb426d7e6bfc2c4f50e5b      0     0      0     0      0     0      0
#> 8a6b88f9515bf2dfac09e6a69045a4fe      0     0      0     0      0     0      0
#> f5151c67af19c8471743ef3e2914f6b1      0     0      0     0      0     0      0
#> c3842619fcbde52781c512c2561fca0b      0     0      0     0      0     0      0
#> b3b00faa9dee9fa6b3b33744f51da00e      0     0      0     0      0     0      0
#> c4f04e28c6267ac6c5b05d403e8f91c3      0     0      0     0      0     0      0
#> 7fc06f107a6ddf6486ef2d43fcf23626      0     0      0     4      0     0      0
#> 8d1dffcd377b338d2a0388637507b592      0     0      0    16      0     0      3
#> 15410d762889d51818986bfcdb5d8965      3     4      0     0      0     0      0
#> 48f68e44e258de63575b8bf6cf565d52      0     0      0     0      0     0      0
#> 4712befe3685661e0fc0968859c369d5      0     0      0     0      0     0      0
#> be9a3727a0422aea147100370d046fa8      5    21     11    40     45    14     23
#> 3c5942781761830f7b5406f551574424      0     0      0     2      0     0      0
#> 291d9f4d3d3c5cc3b0a1804770ed1e77      0     0      0     0      0     0      0
#> 784764c519a64dc9a147025a23222e93      0     0      0     0      0     0      0
#> b9e955ea5254dd68b54775aad2a67829      0     0      0     0      0     0      0
#> 8197f1594e8a17f3fc989fea64f54b8f      0     0      0     0      0     0      0
#> bf4a4f2383a16d2ecf27bf8e162b9ba4      0     0      0     0      0     0      0
#> 4e7d22d793d30bdedbd32a5b73b564b8      0     9     10    22     27     6     15
#> e88959b739c8cceb94b94b7887dd1df3      0     0      0     0      0     0      0
#> 5403445694339e2ae5ecac06d9f3ff47      0     0      0     0      0     0      0
#> ba3943fb075bbf63fb699638292826a1      0     0      0     0      0     0      0
#> 67b2d836db093cead5a8286142669613      0     0      0     0      0     0      0
#> 6c3a0ea164d3eb098dd2baac5519b38f      0     0      0    19      0     0      0
#> c643a776c36bcc00dfb5d790e1ea3dfe      0     0      0     0      0     0      0
#> d53cc81a9d57a64d194c876b76bd66e0      0     0      0     0      0     5      0
#> 7dc00d40971bb454b39fdd9b2dd9aa25      0     0      0     0      0     0      0
#> b4d0cec8543e940be1b5096370fdfb80      0     0      0     0      0     0      0
#> 6c54ceb94bc5e8fbd65224cffc60990a      0     0      0     0     14     8      0
#> b12ee9ca5a2320957213a114788cf2b7      0     0      0     0      0     0      0
#> 5e040c2397f959e7f0dfc6ec2855b812      0     0      0     0      0     0     11
#> 3105a594bba2604ae5ed19f5e0495d2e      0     0      0     0      0     0      0
#> d75c040bdaa1cca96cd0655c9af7ab28      0     0      0    19      0     0      0
#> 9b909edaebabc90b9c1a703626e6bd0c      0     0      0     0      0     0      0
#> a6593da81035cf121dbe46346a4ea960      0     0      0     0      0     0      0
#> 60e4c88400a750a8122b97476313a67e      0     0      0     0      0     0      0
#> 7854a78e3b881cee1c63ea0a50ba29ce      0     0      0    18      0     0      0
#> 1796734b56d45bdc997ccc09ab8731d1      7     0      0     6      0     0      0
#> 8243cd669ae7327f4326be78429ab05a      7    17      6     0     30     6      0
#> ebb93dc32c63402289114859d476cc5b      0     0      0     0      0     0      0
#> 8a01af5a531e851e2c50fc01131cda33      0     2      6    15      9     3      6
#> 03d2cc3f7b2bc352752d39e3738d0af9      0     0      0     0      0     0      0
#> ec378af6e32960b11b5f7391da4ae4fb      0     0      0     0      0     0      3
#> b407f2de868c064f65eb5ad8f82a8d0d      0     0      0     0      0     0      0
#> 73a95b6ba34fcd260989d848fe15d377      0     0      0     0      0     0      0
#> 3abb68202d15781590a2d1f9c5b37e15      0    10      0    17      0     0      0
#> 77923ef8e9ed01c9b6a8ac85eb028163      0     0      0     0      0     0      0
#> 9981b8fd20f400e052527f84f5db1b66      0     0      0     0      0    15      0
#> 597bdec49adf7675912f9fba4e262e70      0     0      0     0      0     0      0
#> 0fd3126118b13fb187bc55476a67c95a      0     0      0     0      0     0      0
#> 2a0c511ecc2513c42fe624e8e99d07da      0     0      0     0      0     0      0
#> d2912fea2dcf2e5b97e56c6a7d24295d      0     0      0     0      0     0      0
#> 2976af6aaf4261b6c274fcc3229a4d65      0     2      0     0      0     0      0
#> f097c45ce419da576aa111c98a4c7272      0     0      0     0      0     0      0
#> 3abc8260d34008f263e232ea25426c3f      0     0      0     0      0     0      0
#> d7e5b6432dc178aaeed9c81e0ec57afb      0     0      0     0      0     0      0
#> dfe71ce78ddfdbc7b453d17d628365f2      0     0     11    15      5     0     18
#> 289b18a86084273ce3b4038bf37e2840      0     0      0     0      0     0      0
#> 844a67032c52dc49019bf5dff248bf9c      0     0      0     3      0     0      0
#> 38cb3e3ce5ac5aed7dec50247ecfb267      0     0      0     0      0     0      0
#> 74e88e71e9929b1ddfcafccfab3b7bba      0     0      0     0      0     0      0
#> 66bc0737fc7fa1463ac7189ab7073998      0     0      0     0      0     0      0
#> b5a7aa588c9cae6ae8ca2f7d87ac53e8      0     0      0     7      0     0      0
#> db0540920457b21b170947aac9137771      0     0      0     0      0     0      0
#> ebbaf260510ab16ac6b8c2b9866d05b9      0     0      0    12      0     0      0
#> c1f461b67419d4a80a61db0867da9ed6      0     0      0     0      7     0      4
#> 0c6eba6452635eecfb865a1049299176      0     0      0     4      0     0      0
#> 7ecabb10f29a0c44eb53e9f9852fe22a      0     0      0     0      0     0      0
#> a33cb3037e7825920073b4dddaaa41e6      0     0      0     0     21     0      0
#> d625812c15c074c2e2b18099f97412c9      0     0      0     0      0     0      0
#> 17ef77f12e2a55717e910a5f6a625448      0     0      0     0      0     0      5
#> 7baf638fdd3f71166aafd22ecf410348      0     0      0     0      0     0      0
#> 9436cda211f60f722865cd23ef75d4ed      0     0      0     0      0     0      0
#> c8ab0b20e740eb6b35e68d58f69736e4     25     0      0    79      0    43    112
#> 5e43f2d8be71399fc8d6736754cb42ac      0     0      0     0      0     0      0
#> 19ad90ee1081a003406b02179ae63709      0     0      0     0     11     0      0
#> 0a25cdbd0660e68970c563c43a18f54a      0     0      0    18      0    14      0
#> 0aa6b210e902d8ab9d6f267407258ebb      0     0      0     0      0     0     35
#> a864c36d0d898cc5abd94d6d6f0bd5be      0     0      0     0      0     0      0
#> dea0dcb305ff0f27e3c36c7133cf933d      0     0      0     0      0     0      0
#> 463977baa69c7724fdc9407ad162714c      0     0      0     0      0     0      0
#> fcb72326c29342bcfb052953466eafc6      0    24     17    53      0     0      0
#> 5a554353f9ed728f0ab0631f29c8be2e      0     5      0     0      0     0      0
#> a9dd27b48e443415a743bf9902f515c4      0     0      0     0      0     0      0
#> 77a40bcf979c0099c52f19a3f1fe5217      0     4      0     0      0     0      0
#> b4007ac525a9b409d70426e24421cff0      0     0      0     0     19     0      0
#> 406ba18a175c4c79efb25f48a6beafe4      0     0      0     0      0     0      0
#> fe1c166755bc722cefd2e2a3b2b70672      0     0      0     0      0     0      0
#> 8afd274849ac07c0ca66f8f5f143df3c     22     0      0     0      0     0      0
#> 197301ed9d8728b5c061b99a17973161     31     0      0     0      0     0      0
#> 14779f9dc3505598ad3d799e4ee676c4      0     0      0     0      0     0      0
#> 4085c246640eee043d523ce2d799023a      0    24      0    41      0     0      0
#> bba09a403e67c1381de01c86511f1057      0     0      0     0      0     0      0
#> 392c59fee4e6b91d7654f2d3e3c65be8      0     0     21     0      0     0      0
#> d2b05ae73ad5ed296d7add9e2bb2cce4      0     0      0     0      6     0      0
#> 8306a9473dabfd5918ea869dbfcd1291      0     0     16    19     13     0     17
#> 06064c372bb624503c2905d13f7691f5     65    98      0   175      0   107     32
#> 3265692341e42dbd5cf4ed1938b906e3      0     0      0    23      0     0      0
#> da96a0797b9881a123795838b9f01140      0     0      0     0      0     0      0
#> 0bdddf10e577657becec88e42001fd16      0     0      0     0      0     0     25
#> 168999a950ebe17e6d2ea1f90e557d91      0     0      0     0      0     0      0
#> a266ffe9d38f626d106fa615da4618ef      0     0      0     0      0     0      0
#> 7706405d5bce7715f01d8b54188def42      0     0      0     0      0     0      0
#> 4b0cfc8899b466bd9992f7b22a5a06be      0     0      0     0      0     0      0
#> b81b2b4a7ca9eb2c23b170648dfaf22d      9     0      0     0      0     0      0
#> 75cf40872201f150cb5b191676f64b78     38    32     40    42     36     0      0
#> c8e5794c58de9326084a62d3282d257e      0     0      0     0      0     0      0
#> b0706a571f169c0bac7b4864052fcaa7      0     0      0    11      0     0      0
#> 0388f969e3527b2938acc676d6757bdb      0     0      0     0      0     0      0
#> 2dfe751390aedfbd9f62dbb77d1f10b4      0     0      0     7      0     0      0
#> 2548dc9697e1f9afe04a9012f40e2f0a      0     0      0     0      0     0      0
#> 2f019615d5141a4ab00b1aab98882a33     28    59      0     6      0     7      0
#> 629e4f0cf1826904ff74907b2507f58b      0     0      0     0      0     0      0
#> 774861f348cec8451102730c21a790b0      0     0      0     0      0     0      0
#> 11c0d6de28ab026a9ae17abe4efdd6bf      8     0      0     6      0     0      0
#> c11381dda0841c0d08cbfcef3a760aa3      0     0      0     0      0     0      0
#> 3c9619d4f9db4ac6ad7151a0c9d86ad3      0     0      0     0      0     0      0
#> 7e2c72204a8c45b1b0ac39e0e6047098      0     0      0     7      0     0      0
#> 8254e0f9f07709b62dd8483f3b25fad0      0     0      0    23      0    27      0
#> e78135f2886d427e678d12c446bf070e      0     0      0     0      0     0      0
#> abec938a8d53e7471f3f476d50694be9      0     0      0     0      0     0      0
#> 4963567d2979bb761291aa98e0b8595b      0     0      0     0     16     0      8
#> d4de9d1453102564f35262faa3952266      0     0      0     0      0     0      0
#> efa4a78c4b993e9d115841c8a02c795e      0     0      0     0      0     0      0
#> fcacd7424bf3543322ea3e9842d14e86      0     0     12    11     17     0      0
#> f941e15ea7f7560191d57abc4d2bedd8      0     0      0     0      0     0      0
#> ac294f5242de3c9b26e75502c48f6383      0     0      0     0      0     0      0
#> 7b6a8e3f8bfa1bcd1b5fc55db048d005      0     0      0     0      0     0      9
#> c9e776ab599fe75db286675b698c9067      0    33      0    17      0     0     29
#> ebe5cb945f6d6d849c07424591f80753      0     0      0     0      0     0      0
#> f116387f77802409103a559e4205df06     18    15      0     0      0     0      0
#> 93f9939b261cddc9a57071b42771c1d4      0     0      0    32      0     0      0
#> 971a5f48b3418449ca7cf3e190a5e831     34     0      0     0      0     0      0
#> 730de992ca06a95e7193a93146a97cf2      0     0      0    77      0    57      0
#> 8fa8b4ac69784c26901e036d13ab2332      0     0     15     0      0     0      0
#> c363c9efaa85023b488b732b97e47270      0     0      0     0      0     0      0
#> 2e2181cc893ca1e1e677e888bf8af441      0     0      0     0      0     0      0
#> 86fa662da6869634e7f90f85e9ee4933      0     0      0     0      0     0      0
#> 9a02177f2c235a6eb62984d5aab151ff      0     0      0     0      0     0      0
#> d6854be498e22b65da3c00572ee178f0      0     0      0     0      0     0      0
#> 015f364f150b230450da0a123e7bb8ab      0    30      0     0      0     0     13
#> af2c1b22758cd6653e8e841497d096e7      0     0      0     0      0     0      0
#> e7f2fab8c0e7c365768a43ef71406605      0    31      0    24      0     0      0
#> fb57cf192bce1876471f2f6181bf4b55      0     0      0     0      0     0      0
#> 6c740ce69c9f45bea684dc9e36054006      0     0      0     0      0     0      0
#> 3b27c7bf279535d8e39664d8dc1d0daf      0     0      0     0      0     0      0
#> ceea3c29ad012c2918fb6aa21a44e8cb      0     0      0     0      0     0      0
#> fa99ecdb27c1cc44fd866ee1da250c4c      0     0      0     0      0     0      0
#> 387cc83573db0e444e81c89a8a4f8f9c     24    21     24    95     54    35      0
#> c7a8d1af5f7f472671f700775de3e0a2      0     0      0     0     21     0     10
#> 98a6a84550e07ea252cfd8a0f51042d4      0     0      0     0      0    33      0
#> 08b809a68126c38db786cafe9710b013      0     0      0     0      0     0      0
#> 8db4cc0a5d2f03fce8efd2b8374ff359      0     0      0     0      0     0      0
#> 3e7c7d38423c93b43a1d9865cf0bc4a6      0    15      8    19     28     0      0
#> 6f02050760c030ed8cd1ed836697fa6b      0     0      0     0     16     0     22
#> 8fc1865cfcae32f72c1b1570299e2ef3      0    17      0    21     14     0      0
#> 38cc138fac95b006929361093d840a00      0     0      8     0      0    19      0
#> f5f8eb63e22362509e2c7c78271305f2      0    10      0     0      0     0      0
#> 0bb238aaee3ccb7dc8dfa86dbc1ba9f7      0     0      0     0      0     0      0
#> ca4546caf07edf5200cf79cea04aebcd      0     0      0     0      0     0      0
#> 893d12123eb7f623489a4f025a15ae7f      0     0      0     0      0     0      0
#> 192c5a817f48368339e519e0a0a907da      0     0     12     0      0     0      0
#> 3467185ca92b7bec2523214605adabf0      0     0      0     0      0     0      0
#> a2aae80d8fb922fe1074d243b59ce237      0    19     12    37     41    15      0
#> f11233f16d4177d38ab1c0bfc649e056      0     0      0     0      0     0      0
#> b219ca7793ab9f5e441f2a68bc76152c      0     0      0     0      0     0      0
#> e6e1848b51e705e55b3805703f088311      0    49      0     0      0     0      0
#> bd074b9bf365d11468bd76268d63b15b     13    17      8    35     31    15      0
#> f381b492d6651c52fd09ce9351e7eb2d      0     0      0     0      0     0      0
#> d9b31e1c30facc24fa778b1d65f5c457      0     0     61   135    280    81     46
#> 14a67c8c8cc31d203eac8109443fab5f      0     0      0     0      0     0      0
#> 026ad0ea53fb202ebb770acf65705d6e      9     0      0    15      0     0      0
#> 19ba9fea04969f975d6dc8d3da20d5e6     18     0      0     0      0     0      0
#> 0fc4426ed7dccc48ce27cfff8083433d      0   299    427   179    118   126     98
#> 45ad90e179842705fa984962f131fd61      0     0      0     0      0     0      0
#> d284afab226d22fd0596221d17db23d3      0    35      0    17      0     0     17
#> 38194655ba17dbfd7fcc1ed42a57c169     24     0      7    26     14    17     16
#> e30be5268bbebe597424bdf171cc1607      0     0      0     0      0     0      0
#> cb7ccdb4956134b452a3cb78d1751c57      0     0      0     0      0     0      0
#> b585436677ac88dd4ac6239655f2a268      0     0      0     0      0     0      0
#> f5d770eff5c9b9b0c859db57d175eb34      0     0      0     0      0     0      0
#> 1b45ffbc88c788dcee2239e68901cf05      0     0      0     0      0     0      0
#> 3fd951cf7cc3a0740a1f20f8a9c6d2c2      0     6      0     7      0     0      0
#> 2027d3e90ccc4003f1954322e99db6a9      0     0      0     0      0     0      0
#> 6157cb92771830955413f808e83fc420      0     0      0    21      0    17      0
#> 6349506cac32056a214b04239bdc3a33      0    53     11    67     43    49     62
#> 2254a8d37cf537b1e35097d70b6ac91b     46    44      0   209      0   136    254
#> 38f3bc8ef8836574de0fcbeafd348b65      0     0     40     0      0     0      0
#> fec23c70ae0c70ecceea1487cd978ee5      0     0      0     0      0     0      0
#> 503d123a6c78c6fb280a6472bd61b51e      0     0      0     6      0     0      0
#> 9f33de8570c639e727465f1ee4d750f6      0     0      0     0      0     0      0
#> d68e1899fa65515135d9bd8d356b7e3d      0     0      0     0      0     0      4
#> e15b88ccac7a4bef64215e915f6e0d89     18     6      0    16      0    15      0
#> 767042b4b19d8bfaaa426dad7dbce378      0     0      0     0      0     0      0
#> 4a3c4abee502d6185abd750522566bf3      0     0      0     0      0     0      7
#> c254aa5c00b7f27d756f7bbe96b3d4ad      8     0      0     0      0     0      0
#> fe8d8af791eb6e5d4335bb7604c74e9d      0     0      0     0      0     9      0
#> d81800d7e31b69877ce1eab9858b8f3a      8     0      0     0      0     0      0
#> c55c7cb80c5700ff9aded76192543b5a      0     0      0     0      0     0      0
#> 91dee7b3b90bbf2e53b2472b9bbde0e9      0     0      0     0      0     0      0
#> 6b28df61cb5dee74b76b1a13aaa0027d      0     0      0     0      0     0      0
#> 9a830fa26205ad91a02b5bcaebb0221a      0     0      5    22      8     0      4
#> a54b8142d2a09d058b81c377e8e286b8      0     0      0     0      0     0      0
#> 69ef3a03a4d733dbffd99fdaa1d2e62a      0     0      0     0      0     0      0
#> 083cbb131dd305a88b7d338bc371c35f      0    13      0     8      0     0      0
#> ef1a92caffb4d7cae286fdc77b0a93a6      0     0      0    21      0    17      0
#> 2cddd5dca166486f50a0828b28001e1e      0     0      0     0      0     0      0
#> bb899d1dcfd89cb200a76cb8951904f8      0     0      0    13      0     0      0
#> 7817f6851e9df03a3d992cf915310dba      0     0      0     0      0     0      0
#> 54eaf6571de743a86ccb531bee8aad7e      0     0      0     0      0     0      0
#>                                  6.16RI 2.7RI 2.11RO 2.7RO 2.15RO 2.15RI 3.14RO
#> 374a5eb6496d14e8e9d6f3d8c7e34d9d      0     0      0     0      0      0      0
#> 4c47d5cf81df1ac6e2a0560029a81578      0     0      0     0      0      0      0
#> 8fb6741a6685fad0374f85ad3e16156c      0    22      0     0      0      0      0
#> 310f3de009c95de5c938b8e811dfc96f      0     0      0     0      0      0      0
#> 7d3b39e8c0f0fdcb81b0c898a98c38f9      0     0      0     0      0      0      0
#> fc211549300b0954dbb4a4bf57e8a605      0     0      0     0      0      0      0
#> 53fad3d8c52f4b5a023ff92aee9c0f1a      0     0      0     0      0      0      0
#> 5e37342d650cde6cfc5bbe2726b12e3e      0     0      0     0      0      0      0
#> 56025f74cc8e54d4a054442d786ef096      0     0      0     0      0      0      0
#> d9e51d98f1c86a1c426552d252c944e4      0     0      0     0      0      0      0
#> bc905dcc609da83d3620f1986e63a9f7      0     7      0     0      0      0      0
#> 90009df763e86e6eca7baa54f01551e6      0     0      0     0      0      0      0
#> 6bbd12b392461134369d835715cd2765      0     0      0     0      0      0      0
#> 83c88a72c2dab171c38aadd98aa2d389      0    20      0     0      0     15      0
#> 2b93708de9f1c24ff6140814a599a3a8      0     0      0     0      0      0      0
#> 671a53a2acbf0c8c7e240e446daf0000      0     0      0     0      0      0      0
#> bfe35644221421300c0c17e00e465c05     46   111      0    23      0     76      0
#> d12bdf20f24b92aa8ec94c6c3f60497f      0    20      0     0      0      0      0
#> 3bdb72ce060bd89f107c5c2bc95399ed      0     0      0     0      0      0      0
#> e6a82be4b75aac8b928c1dfbe9ee9287      0    22      0     0      0     36      0
#> e3a63ba0a0ec40fbf270a57c180e0d3c      0     0      0     0      0      0      0
#> 2877cf981b128a18f788ba5437c1afa6      0     0      0     0      0      0      0
#> c015bc50696fe3e3d94ace01d9dcd691      0     5      0     0      0      0      0
#> eb2eeddf829080f8d9f25c6be5d507cd      0     0      0     0      0      0      0
#> b1fcda6d8df4c25c8e6d7124c076b1fe      0    24      0     0      0      0      0
#> de8bb5c39fb121fd048dfe0d639e6e8a      0     0      0     0      0      0      0
#> a234ac03223c275785dedf6df3b2aeff     13     0      0     0      0      0      0
#> 9ea1641eeed1f51fc228c653076d1d08      0     0      0     0      0      0      0
#> 694b926113eb8291fea166446890c8fe      0     0      0     0      0      0     15
#> b2af163540d17007833dc819704afb28      0     0      0     0      0      0      0
#> 6fc0ba1f8ff8f9259c9d492e59771755      0     0      0     0      0      0      0
#> 7b8c1dba059eac7ec7f3cefa203f6b2d      0    26      0     0      0      0      0
#> 0b06028193ef1ce19f5339a8ba7d4448      0     0      0     0      0      0      0
#> 56b4ec05ef3bbded515ec3dd2d573197      0     0      0     0      0      0      0
#> fc9c84f8767611251f0eb93af4ca10db      0     0      0     0      0      0      0
#> d3616bc0d8925a452274c484c73406c3     80     0      0     0      0      0      0
#> 7986b8b93f3afce920bd19a5d2229c09      0     0      0     0      0      0      0
#> bc78ee7297b6912452be0865d113f7ad      0     0      0     0      0      0      0
#> c472668b46b0d0f575108aa4170dd81f      0     0      0     0      0      0      0
#> 2658a78c5ff06760389ce2e3a548b4a2     48    45      0     0      0     46      0
#> 600c9e81fd1d1ddd5a37a5aa0f9f796d      0     0      0     0      0      0      0
#> f8d0bb97f093d8e8630ea62ac8a1bf13      0     0      0     0      0      0      0
#> 6c685b83d63d6d56a9dcaba34215ae6b     15    13      0     0      0      0      0
#> e0e7e6208280a5846a8b66f00667f8ca     42    84      0    18      0     56      0
#> 6d1aa1c1f93c423a856a5dd3b62e5c70      0     0      0     0      0      0      0
#> e553e766ea0e08f1f9fa2debe7946af6      0     0      0     0      0      0      0
#> 0f1ef43f6b23a567aba80599766d70bc      0     0      0     0      0      0      0
#> 3c4ccaa45e8d148c9a17b5d0f3720cd2     52    89      0    84      0      0    104
#> d9df6265d046a558d87d397b8a177d54      0     0      0     0      0      0      0
#> 4bbad82a39655a26a2de61b88c5af7fa      0     0      0     0      0     65      0
#> f5a6865451a3b138f4915c5ab445a0a7      0    21      0     0      0     13      0
#> 4a86fc550fa1ac37aa480d37edea875f      0     0      0     0      0      0      0
#> 02fc21e72f66adea757b355edea4a369      0     0      0     0      0      0      0
#> 25cf10f4f15dab00925cef52b9f08442      0     0      0    11      0      0      0
#> fa8c8b98d197c252368c160b43cefb98     34     0      0     0      0      0      0
#> c0b07d29aed1a875641aede53bc0a20c      0     0      0     0      0      0      0
#> 341433c1c8cc65ce47dcf32179d08fd8      0     0      0     0      0     13      0
#> 8bf175f329d16e4aa732cf2b32279df3      0     0      0     0      0      0      0
#> 062f7552b747dccab9586dff7275b5d9      0     0      0     0      0      0      0
#> 059e124ad6ba8f98c3d21db28a075ea8      0   119     78    78     46    114    121
#> 3d7a81faea3f783158352c9f739c8d17      0     0      0     0      0      0      0
#> 90be2fb2c82a98002d33ed7531e877e3      0     5      0     0      0      0      0
#> b0e263cdd1c189fc1757d8189db14be2    118   189      0     0      0    126      0
#> b962ce332fea338b05293d1a5a29cbca      0     0      0     0      0      0      0
#> 7a310a8f8e4bcbfb717263bd5376994d      0     0      0     0      0      0      0
#> 57136cbdae8fe066c2bb4a661003719b      0     0      0     0      0      0      0
#> 69fe46cab1cc9407470a2397bfc84103      0     0      0     0      0      0      0
#> b88efba97658a6936bb053985d026ace      0    12      0     0      0      0      0
#> 4d691fbbe2dbfb2815649fb72c757977      0     0      0     0      0      0      0
#> 49192cc8a222d97b965b839ff5d0f333      0     0      0     0      0      0      0
#> 859873fcd088c294a92a13340744eb6f      0     0      0     0      0      0      0
#> 0a8c44ca4d744e313e6dc5b239fca562      0     0      0     0      0     10      0
#> 09f94345824c82e8fe0d824cc4478fb5      0     0      0     0      0      0      0
#> 9dfbdb741fb0dc25d2206b8680bbf579      0     0      0     0      0      0      0
#> dc2bea463c720874cccd189ccdcf0e3e      0     0      0     0      0      0      0
#> 6273eb52991d528cf5281215f79b40a3      0     0      0    23      0      0      0
#> dc4bcbe74986e44ec3a040b866dfda9e      0     0      0     0      0      0      0
#> 96ec4781ee61f4083792f252f1b6cd3a      0     0      0     0      0      0      0
#> 9657d52acca51e701e7a3b7cedcfcc12      0     0      0     0      0      0      0
#> f65f183293d4710041057ab98f2c94ed      0     0      0     0      0      0      0
#> 8f6b3783b31a2dd8640de3df92dbb9c4      0     0      0     0      0      0      0
#> 0d90eb4842c9ad9126ed33c538c768ac      0     0      0    27      0      0      0
#> 8d3ee46aaaf7728594014be61c649347      0     0      0    13      0      0      0
#> e011a95fda1f1f2bafc4e5b83b52d989      0     0      0     0      0      0      0
#> ef8ea65ba4c2e2da3ba9556a87d2385d      0     0      0     0      0      0      0
#> dbbb92eb3a8a3feaf96092d46f2b6a45      0   109      0     0      0      0      0
#> 06cb03cddd8a12dbfaf27e8984b395c4     93   229     99    68     61    147     44
#> cdf134bb501c54d997945dabd44e35ff      0     0      0     0      0      0      0
#> f3c4561cc01a3b45ee3f134d709c6b94      0     0      0     0      0      0    127
#> b121442fd1eefe78ce4f4602aac1842e      0     0      0     0      0      0      0
#> ba8707aa88a2a62939d9338189b2a733      0     0      0     0      0      0      0
#> ecf7289628d49f9244469337ee1622cb     43    64      0     0      0     53     13
#> 9f06b400b93f859b11874b0bd7f09295      0     0      0     0      0      0      0
#> f775be5cf7c0dd2177006a44913503d5      0     0      0     0      0      0      0
#> e87ceb10f2becb2a6c4cf0692b1923ee      0     0      0     0      0      0     10
#> 036574bb8e05f3b9a11ab7a4090d6950      0     0      0     0      0      0      0
#> 0267ba9db97e00ae01fc2cccbb6263e0      0     0      0     0      0      0      0
#> be4f701c1dbfbda16c5150af29430141     45     0      0     0      0     43      0
#> 655a66e6566798ade91f73b2942c0b63      0     0      0     0      0      0      0
#> cd79f2696509d0ad906f5addba264079      0     0      0     0      0      0      0
#> 07ea3284bda2bf7daf1251ad2744e4c4      0    77      0     0      0      0      0
#> e3870da9de4a2b18a7ae9e9cae208e3e     29     0      0     0      0      0      0
#> e00bd59c1ea0ca9b1cc5f976801d4e1b      0     0      0     0      0      0      0
#> b39b381d0f533c9923c3b898ffab5614      0     0      0     0      0      2      0
#> e2959f0552fa6df5707d9ec8b57816d6      0    32      0     0      0     22      0
#> 53157acd1642f973842d01e1cd287421      0     0      0     0      0      0      0
#> 5991790d1c0aa8b73b2d7f906449a1f9      0     0      0     0      0      0      0
#> adee3552c211d9e33aaa28e82532bdea      0     0      0     0      0      0      0
#> 0089040d041888e3ae4aa4e0010a0e36      0     0      0     0      0      0      0
#> 8e876e13d051992c4c8616cf2db7d73a      0     9      0     0      0      0      0
#> 179abbc5c1f1549ffd0d7ae6963e3783      0     5      0     0      0      0      0
#> 340d095f9a9fb445582b60e1282d985d      0    25      0     0      0      0      0
#> 3d7e5432d8eaa3fa56ab7138605ea43d      0     0      0     0      0      0      0
#> d4b233740db907d69035915e8657fa8d      0     0      0     0      0      0      0
#> b22119df3b1de76745abaac585b646f8      0     0      0     0      0      0      0
#> 64761a24fc66e7c93145217dd153b1ee      0     0      0     0      0      0      0
#> eb6d7f9e0601207fdcfc4a88947d513c      0     0      0     0      0      0      0
#> 0d20271e3e2dbc68e1de868599a9a1ae      0     0      0     0      0     18      0
#> 369ea22772cb774950c830d0af7e8235      0     0      0     0      0      0      0
#> cae63066cfe088a49fcac50f119c9e2d      0     0      0     0      0      0      0
#> d9c57ec20ed90e2cebf6cacc4f9a5511      0     0      0     0      0      0      0
#> f5828cc3078b4c219a44c67e96e5d455      0     0      0     0      0      0      0
#> 73847ac2778d4cba27fe9758ca129d6c      0     0      0     0      0      0      0
#> d78003f13712da48d0bb36fd45116b91      0     0      0     0      0      0     24
#> eb40853d986ad72c82658d939f291ff2      0     0      0     0      0      0      0
#> 5d675a3518222fc99c8cba34feaeac81     63    92     34     0      0      0     52
#> 3ba390ae11c4e76375af082eefe7e41f      0     0      0     0      0     13      0
#> 3aa66daa4e560cf62c1833697bcd3837      0     0      0     0      0      0      0
#> 6f4388510be5ea70356f2257d8c57ddf      0     0      0     0      0      0      0
#> b33390d034148edf07d8a742778ddcbd      0    30      0     0      0      0      0
#> 09225ec7044394bbd3b6811a6799e481      0     0      3     0      0      0      0
#> fe6fbf9b86e0f487911a856009fa4a1e     75     0      0     0      0      0      0
#> 0aee6329e0b5b81449980151bf40205a      0     0      0     0      0      0      0
#> 0ef66a78236f90916afc2305f3ff6c80      0     0      0     0      0      0      0
#> fb1eb41359005c2e3e308fb668f07220      0     0      0     0      0      0      0
#> cb59c57800227006aead644df42b8986      0     0      0     0      0      0      0
#> 138da9dc244c3292536cacb8df79467b      0     0      0     0      0      0      0
#> 72f92361ec0696cb8df985b93d277a25     21    19      5     9     15     40      0
#> 4d5ab121763879fd59e8f00d35473540      0     0      0     0      0      0      0
#> 20032b2b586fee22805eccebe3059170      0     0      0     0      0      0      0
#> 710324ccc8cafa6fc9b38baf8299d8aa      0    27      0     0      0      0      0
#> 523e6b60340f4a215877b1c336822e50      0     0      0     0      0      0      0
#> 4501f98c83484fd902be11383c8d032b      0     0      0     0      0     19      0
#> d4cbf23757aad57c91b058cbb05566e4      0     0      0     0      0      0      0
#> de7964d5c22112693ed5796bb9f04b0b      0     0      0     0      0      0      0
#> f4a283be42989b920b83ea73ec72b25b      0     0      0     0      0      0      0
#> 2f3111649d0850ab86799f8bbadc28f8      0     0      0     0      0      0      0
#> b8acdd45cf33ab75fe4c0851a1787dc5      0     0      0     0      0      0      0
#> 5bb7fec39744a6d4e4d61596dabc6886      0    31      0     0      0      0      0
#> e8eddfb37c4215c86813163a27f202a5      0     0      0     0      0      0      0
#> 826e188c911c0846fd927300127c80c0      0     0      0     0      0      0      8
#> ff1fe599c44074ee6bdf919a119df93c      0     0      0     0      0      0      0
#> aff00055d3876f624d67fb80bdb59934      0     0      0     0      0      6      0
#> 1beded194f3ae4b57715b45e50adef46      0     0      0     0      0      0      0
#> 908d91feba7c1f1d8cc8a7145173d942      0     0      0     0      0      0      0
#> e7095743695f6d4337d01ceb9e29e224      0     0      0     0      0      5      0
#> 4cd420962855afdf0bc69bb268b846c7     47   125      0     0      0     80      0
#> 430ba5bea0fc08f65b69497bfbeb5c4b      0     0      0     0      0      0      0
#> e12cb2f0d2ce4401413176007dcd4158     62    85      0     0      0      0      0
#> ec615213f6c7febb65c9f2e45e438111      0     0      0     0      0      0      0
#> 412dd52cd54d0244ba8fc42ad140acda      0     0      0     0      0      0      0
#> 89192d9f779341d313f52aad0dc6a06f      0     0      0     0      0      0      0
#> ed543f92c61a06fc4e703655ce910a2f      0     0      0     0      0      0      0
#> fea1a70b0ef0f72e3bf27d8b295331a8      0    27      0     0      0      0      0
#> 6afe805b36f5aa29e7f2096aac42ace1      0    41      0     0      0      0      0
#> d9cc5a46831f24ec9192ae82dbc38b0a      0     0      0     0      0      0      0
#> 0d3eb38eb502c4ad13bd5a7060295411      0     0      0     0      0     20      0
#> 69a0ff5538a5490e291c75ce70325066      0     0      0     0      0      0      0
#> f91d545c815effa3c5afa83a8ad3397d      0     0      0     0      0      0      0
#> e6b032dea8e13c34af958fc48af09b58      0     5      0     0      0      0      0
#> 84a1b89045ab6285328f366224f3ec87      0     0      0     0      0      0      0
#> 4a6f0076e4ac71474b3c364294ad3b33      0     0      0     0      0      0      0
#> 7bdd39e4e3ccc7000862b487510e67d6      0     0      0     0      0      0      0
#> 998aedf12cdf054f5e504a2e21d20498      0     0      0     0      0      0      0
#> 567c2640108006c8eca081ad4180e346      0     0      0     0      0      0      0
#> e50a6625b9176b6856764b213013150f      0     0      0     0      0      0      0
#> 8d5a5ab7c784b0325a4a86a16fbba80e      0     0      0     0      0      0      0
#> f58ae28a309fa3d5898d2eb46d9eb9e0      0    23      0     0      0      0      0
#> 7ac65a8fd21111b68e7c2d6205112a83      0     0      0     0      0     11      0
#> d7add64c8700f5a8cadf40c8c7f5ff67      0     0      0     0      0      0      0
#> 6bad9d981d3cf0f90262547180dc6371      0     0      0     0      0      0      7
#> 149c0e63c83113e1246c0e041de67d9c      0     0      0     0      0      0      0
#> 6f196fae54ea90b5c11627c503c98940      0     0      0     0      0      0      0
#> b67c995d842e5b4b92fa0ce052603ef1      0     0      0     0      0      0      0
#> f44ffae67c7297e2c696b75fbd4e0f4c      0     0      0     0      0      0      0
#> 33997912cd952feaec18fbd11b31acac      0     0      0     0      0      0      0
#> 0a035423c4bf71a17566549f9b69448d      0     0      0     0      0      0      0
#> f65a6f1f1e5a4bc5c733ba364c0096d3      0     0      0     0      0      0      0
#> f15389a9f245465c33db3de72aeb17ba      0     0      0     0      0      0      0
#> cd2b676d3c0785a60359c0d416701ee4      0     0      0     0      4      0      0
#> 70ba2a378ebcd9a8f57b024a561cc6b1      0     0      0     0      0      0      0
#> 0273707384afdd41e9d022b539094211      0     0      0     0      0      0      0
#> f9ac98e876ecdb504202388d7ccdda5f      0     0      0     0      0     13      0
#> bdd9b87b782b8f935ecac62a9824e6a1     16     0      0     0      0      0      0
#> 1f12b97495668f1fc378a349546d1806      0     0      0     0      0      0      0
#> 83f1b4d776f6e08a6ba6f86634a1e5a9      0     0     15     0      0      0      0
#> 49f614420449f6c7ecc8c9f4e137ad05      0     0      0     0      0      0      0
#> aba5aa00a1b262f27517eca64ee3333f      0     0      0     0      0      0      0
#> 8971980484eeabeef00145b93b699e30      0     0      0     0      0      0      0
#> c1347ae54b0fc1dd385a1ccf3b8061ce      0     0      0     0      0      0     57
#> e9fe19783ebcc7e43a3a07ca8225db62      0     0      0     0      0      0      0
#> 7c523a5f0a7c736efd82afef84137901      0     0      0     0      0     25      0
#> e13c07ec336c89db5183cc6e9e57ee6e      0     0      0     0      0      0      0
#> a017ce431ccadea609b501d46e64e55b      0     0      0     0      0      0      0
#> fadf38d03d4b944f5782a759591ec3ae      0     0      0     0      0      0      0
#> 41bfda1c7f32205246e73a1e3e3ec4f4      0     0      0     0      0      0      0
#> 9e4f3e697b70758a7a996f3dd4533650      0     0      0     0      0      0      0
#> 39c4f6caa7510ed108f3ebcf72f79edf      0     0      0     0      0      0      0
#> 4aa083fe7c1b0639dbc48ec649e68a31      0    14      0     0      0      0      0
#> fd2ac2e14b72b0e22a052e3ab4580471      0     0      0     0      0     24      0
#> 4ee10294799bcfeab0c8ee567d5d9608      0    12      0     0      0      0      0
#> 46869f12c31de7a05781e930c25a1ccd      0     0      0     0      0      0      0
#> 5e2548d6a970645d0ca943a72d207489      0     0      0     0      0      0      0
#> 0538e1beed69f098c3f602130ee11a1c      0     0      0     0      0      0      0
#> 991c563278c4e8c91ffebce3560fe51b      0     0      0     0      0      0      0
#> ce05f58164e61517198f7ab0e3304cb5     34    17      0     0      0     31      9
#> fdbd0fed21bae91c05b2da268a025a89      0     0      0     0      0      0      0
#> 7412d11595d77f26792071cd28a26760      0     0      0     0      0      0      0
#> 4b087e13957eb38128f63e1c99df23c3      0     0      0     0      0      0      0
#> 310a6f0e6841db2ce75005cac942a972      0     0      0     0      0      0    153
#> 3df0fada803ac582802b47d4f4efe0c1      0     0      0     0      0      0      0
#> 3cd4e34c02ffdd1446df65bb54278c57      0     0      0     0      0      0      0
#> b14a1dea93554257791f16c2e94906b0      0     0      0     0      0      0      0
#> 902a4435fc3af101e8a9b5f51550fa4d      0     2      0     0      0      0      0
#> 1995185b31349529f58964c2e5a3becd      0     0      0     0      0      0      0
#> 69b247c565afdeacee0240482b1ce536      0     0      0     0      0      0      0
#> 95ff1b9528acca018e9c4697a8625504      0     0      0     0      0      0      0
#> 132d26d3ce6ac2ab25ccb4dd171879d7      0     0      0     0      0      0      0
#> bf1a6e4f94ee8ba60cf79089139641b6      0     0      0     0      0      0      0
#> 704189aa09f64b992b3a02ee001dc706      0     0      0     0      0      0      0
#> dea5ea27839ee67a31d18ba0bf0f7f65      0    27      0     0      0     28      0
#> 3840a4fff0250232f8570cae7c3514ba      0     0      0     0      0      0      0
#> 31ecabc032997aaf41dd0d8def3ec93f      0     0      0     0      0      0      0
#> 024a0e0348ed9abe9ae1c48b1b7ba09e      0     0      0     0      0      0      0
#> fc1d9940419420113ce3fbfacc8d703a     43   100      0     0      0      0      0
#> 0770fbf28cd307a08bf557af7ffecfd3      0     0      0     0      0      0      0
#> 883cbe8e4d1d98f468355e8ea974d222      0     0      0     0      0      0      0
#> cb6d6561383fb68f43f92a88a607749d      0     0      0     0      0      0      0
#> e07e0f1a9fbc696a9fe4f9b7739ae45b      0     0      0     0      0      0      0
#> 73d8d8f999acfa13cba2277f7526897e      0    67      0     0      0      0      0
#> 5218274c354d41965abbf4b0751d6f3a      0     0      0     0      0      0      0
#> 5c31570f706752fa4f7851aca5b5c291      0    14      0     0      0      0      0
#> 6a5e71c8b86a52ab1dd5e61164103940     26    54      0     8      0     42     11
#> 25641b0b803c7d8493ebc02cb104bf3d      0     0      0     0      0      0      0
#> 8844d70dd1fdb96f2a7faa65bccf78c6      0     0      0     0      0     16      0
#> 371b855afa8bed030e86ec78241b01a3     12     5      0     0      0     10      0
#> d0c180a47378d414a1d088ac14db236c      0    22      0     0      0      0      0
#> a072a588b3940fab3b227e4901c40f7f      0     0      0     0      0      0      0
#> 977a3523051e6b23da6221dbde795b34      0     0      0     0      0      0      0
#> 1c3adc87d2e93b08e576953508eb680d      0     6      0     0      0      0      0
#> 74eb0ddf379b00c5032b18fb9fbf1f32      0    34      0    10     10     61      0
#> c335b34a15db7d532f511e6c952abe62     10    25      0     0      0      0      0
#> 0c78285397fb853a0e39069213510c54      0     0      0     0      0      0      0
#> 779dec6de9c6a8b7df28fc8ec3c75bee      0     0      0     0      0      0      0
#> 5526624b473bb1f0541dc6fb95c98678      0     0      0     0      0      7      0
#> 649e7fa7b9ce5c5472d34675ed643639      0     0      0     0      0     17      0
#> a9c5addeea542953d1bae4d30d599b02     13    31      0     0      2      0      4
#> 68291fb3b2558d438da4810961c8d53a      0     0      0     0      0      0      0
#> 9c5d92bfefc7a190b93d5b9a6a77697b     22    28      0    10      5     22      0
#> 3ddab42d701754625cab0c236c5b6e6d      0     0      0     0      0      5      0
#> 287e840f9ecaef80d7a4762c8d6bf401     30    28      0     0      0      0      0
#> 885f411fa54c5576148678f665bb053f      0     0      0     0      0      0      0
#> b5709a57df9cbc8a741bb78a36f29c06    108   169     24    50     35    171     71
#> 26698cc60893017ef0b1ab173688d3de     50    35      0     0      0     24      0
#> db64f0028b6518839ab7129c4b7ef72d    215   267     37    73     55    322     75
#> f46fb43262da299ec7f0084298cafa50      0     0      0     0      0      0      0
#> c45e7087de669be270881fc8ace0fad8     58    47      0    44     28     93     40
#> 854a20d5019388e81b1e78419eb017fd     69    31      0     0      0     46      0
#> 9a2325ea015a8c5944047e3349cc1565     24    64      0    13      0     43     34
#> c5305ae3c6adc90814d70ab48f4f071a      0     0      0     0      0      0      0
#> e230123ec3ce8654942ee1a000a80010     70    62     23    22      0     52     27
#> 14e3c60cb885509667670151bc4a1df0      0     0      0     0      0     20      0
#> 7fd7b7bece92bbbfcceb2dcf1e835a8e      0     0      0     0      0      0      0
#> c098d1fb769dd975c2332dfb70da489e      0     0      0     0      0      0      0
#> 95dcdf1eb06b252e67dcdbdaffc505e8      0     0      0     0      0      0      0
#> c5399d0f9eda8f0815fd28a257952496      0    46      0     0      0      0      0
#> c7388cc314e992819b8f97ea35280982      0    35      0     0      0     22      0
#> ae8facc45d372f7108bd9499f74d1bfc      0    36      0     0      0      0      0
#> 16caf6f5bc653fb6ba668eb821076a21      0     0      0     0      0      0      0
#> e7e283b72ee15bcf0877c5eeb6137ebc      0     0     38     0      0      0     36
#> 1b3a779b384c73b33207dccd14b14fe6      0     0      0     0      0      0      0
#> c8d1afdf7b395ffc22193c5dc2fd5dc3      0     0      0     0      0      0      0
#> 0e3ffc2a68378086e2480831acd1a025     15    26      0     0      0      0      0
#> 0ddf7339b608a340252b33f63fcc0019     15    23      0     0      0     23      0
#> 5fd4b2f7dee3fb0bde6f856ae266dbcb      0     0      0    12      0      0      0
#> e2f4786bc66410c07bd8b96b5d0e2d92     24    39      0     0      0      0      0
#> c30aae8f5216c9536492d0beae74950f    103   186      0    19     36    160      0
#> 9087ecee0806a80861788a17dd7e6809      0     0      0     0      0     24      0
#> c29bc08418318360592b2b2717376a0a      0     0      0     0      0      0      0
#> f791c7db1a78d35e84803b5e03c01965      0     0      0     0      0     23      0
#> 44b44fc4055a58c01f6d1e9de5e0d1d9      0     0     80    82      0    120      0
#> 342df9401103d274851d3bd80704508d      0     0      0     0      0      0      0
#> 0590e6672bb2c55c0ff54999a2b6dcfe      0     0      0     0      0      0      0
#> d4e9147608bf883fdbbf9d23b17028d1      0     0      0     0      0      0      0
#> db24794eb8667aed224401386223f671     25    42     35    24     17     80      0
#> fe41e6fd9770030a21418b288930710f      0     0      0     0      0      0      0
#> 699f18f3b155bff97de941f2cfbbbadc      0    49      0     0      0     40      0
#> 13870cf319997750e3b96bbb0b2f1aea      0     0      0    57      0     16      0
#> 7ca4fb1b1740406241acb862445cf311      0     4      0     0      0     29      0
#> 9c0569fb95451d009fd3754303eaed12      0     0      0     0      0      0      0
#> 1ac26a3c12288d8fbadd6aa61366b1b7      0     0      0     0      0      0      0
#> 89f3ba99552605501b8acb816d288afe      0    56      0     0     13     37      0
#> b4051b5227746f04c1feb9d8979299ff      0     0      0     0      0      0      0
#> 087d358aaf77f81aaf50a623828865bb      0     0      0     0      0      0      0
#> e451c4b3ff1c3dabeee296529b3a7a07      8    30      0     0      0      5      0
#> 012871784f77c9dca7d446ddef2048d1      0     0      0     0      0      0      0
#> 95c6ce2788783bf2f83343deba23b8bd      0     3      0     0      0      0      0
#> 213227f244d97553a9a14c4e3fd04c98      0     0      0     0      0      0      0
#> c0e5e55c5cfbd576e7c085f1e4b69638      0     0      0     0      0      0      0
#> 51e98451d7e3f46127481ec5490f1c91      0     0      0     0      0      0      0
#> c75285a6c61d8fce20d823769c25c12a      0     0      0    16      0     14      0
#> 177966f7619eb24fdbce2266f5bed39a      0     0      0     0      0      0      0
#> b89347f0534679f890628364d2c7b4fe      0     0     27     0      0     20      0
#> 8dff1cb12c0beb5bc7b2546480a6a6c7      0     0      0     0      0     21      0
#> 84f2777cfcb0fbed6281917fcc765022      0    53      0     0      0      0      0
#> a5e4d106a2ab947984fdef187996a0b9      0     0      0     0      0      0      0
#> b6cbf6077b64acbb718971a1ffc97a27      0     0      0     0      0      0      0
#> 27e29495fda424728fd6a80b83a16772     34    40      0     0      0     38     11
#> 63c44ab3725dd0b803d43a29db9fd7fd      0    39      0     0      0      0      9
#> 9c931ab5bf18ef9b10996d59d292e6f5      0    72      8     0      0     24      0
#> 12fd8efb82c3d4909543f6db9853134e     11     0      0     0      0      2      0
#> a8a60c20f7ab7d0f2dadd9f8858c19be      0     0      0     0      0      0      4
#> d57e09eb1645f1abba6074d3bd7d4c52      0     0      0     0      0      0      0
#> 7e0f2b982ca2d35b6ac2cf251d78d3c4      0     2      0     0      0      0      0
#> 3e2e8e938a68a0318bcf7145ad5ed4a9      0     0      0     0      0      0      0
#> f2117b49bf57d9451bbef6de5b762483      0    16      0     0      0      0      0
#> be24dfe003725a0b0f28f4136106faf8     11    24      0     0      0      0      0
#> e63301d7d0af55d5c43db051dc51522f     11    29      0     0      0     16      0
#> a6e2e5501507562f2174f65e606b2d43     14    14      0    13      6     17      0
#> a06aa93cb50b9a39941a11c7242dde96      0     5      0     0      0      3      4
#> 6a0a8abb36aa3bfdee14edbe97e9d781      0     3      0     0      0      0      0
#> 4d3def5ffad0c392381591440bfc3668      0     0      0     0      0      0      0
#> ce437b81e7cf5ce986f7240c26464b0a      0     0      0     0      0      0      0
#> e57788f9d2fae76ed9ca566ad4d475a8      0     0      0     0      0      0      0
#> 8c525d0f63137bd068b75c57f5d04c09      0     0      0     0      0      0      0
#> 6d670aa40dafa152c2d59cf2be9bef74      0     0      0     0      0      0      0
#> e08f8e9e9ffb10a702f2ef617f70eecf      0     6      0     0      0      0      0
#> 7def1c807ca0767ccdfe5a6001dc51f1      0     0      0     0      0      0      0
#> 3b9a17275147adf1dc8c42eb00a06a36      0     0      0     0      0      0      0
#> 24be0b7e4c85ac53b09fddbca916ab96      0     0      0     0      0      0      7
#> f521a47d3bb4acfec6b261f880b4eea8      0     0      0     0      0      0      0
#> 91c24c9b6ffbedcc257d164176060fe0      0     0      0     0      0      0      0
#> 9b45f210647c5fe71b005eb9e25299a5      0    11      0     0      0      0      0
#> 974077d25fa6b32a7af22670b142a8f1      0     0      0     0      0      0      0
#> 203aa2535e54a11e4c48af4ab946686b      0     4      0     0      0      0      0
#> 9064afe569448558eb08e5bf5e458bce      0     0      0     0      0      0      0
#> fad432eeeff7f1d923e6d60aafe18dcd      0     0      0     0      0      0      0
#> a29db0e5f4a04c5c6937745ccccb167b     32    43      0     0      5      7      0
#> f179c10f8e86873c90ceb14b974a114e      0     0      0     0      0      0      0
#> a23e8588a2240d7af6f65519158c902c      0     0      0     0      0      0      0
#> 52034e0da6d5e7962369437956984e18      0     2      0     0      0      0      0
#> 522d20c27cf4237b60f36ffdca056c19      0     0      0     0      0     19      0
#> 957615aefb421e6ad4728e610b2bdfa6      0     0      0     0      0      0      0
#> 5ccc74e70b8a82a6eebea3385738ece3      8    10      0     0      0      6      0
#> e5924fd11a48419c6ec0e2a7cb4e43ee      0     9      0     0      0      0      0
#> eef7be4b2531cb7e909b58e91af5640d      0     0      0     0      0      0      0
#> 4416cbc2026260fa4e2ed0a8feff1443      0     0      0     0      0      0      0
#> 03f6091f5aacb69f06fc925acd1c9506      7     9      0     0      0      6      0
#> f655e5ddf5074d5e71bfa583926dec30      0    14      0     9      0     10      0
#> 9f733d39e588c7310499dd444cdd7812      0     0      0     0      0      0      0
#> 2b413de5e05a1385bab229f50879628f      4     7      0     0      0      0      0
#> 178b596f2fe5cb6cefc5513917988d84      0     0      0     0      0      0      0
#> aa94b014c726cf40983045bd606241c1      0     0      0     0      0      0      0
#> 55310b3aac4868c4e33b0a9492d7609a      0    12      0     0      0      5      0
#> c5d8d828df0e1baf5c2db05a6b8c4476      0     0      0     0      0      0      0
#> 6ceb26f60606233a7728d3457af84380      0     0      0     0      0      0      0
#> 1ff83868d1c9ecf52d9b776e0b3fd896      0     0      0     0      0      0      0
#> 3edfccfb9f81b340854e84ee5102a4c1      0    22      0     0      0     14      0
#> 963dc552c551872800f2fcbc5663e5fd      0     0      0     0      0      0      0
#> 3fda8c20e5da9fbc3fb7c5c864559019      0    11      0     0      0      0      0
#> 6fd88fa7e883bdead3c5ba89da853e8d      0     0      0     0      0     17      0
#> 35826fe826b7a9f202e8bef7b39d9a47      0    21      0     6      0      0      0
#> cc514f614ba146f489bf62b5a56074c5     37    40     12    21     13     91     11
#> 26260de3858081ce55e6f220b3ba25df      0     0      0     0      0     18      0
#> d6727a85782ffd2a069dccc23adaf521    138   193     30    42     19    100     32
#> eb1a74c2ca8831ef945775de0827ddda      0     0      0     0      0      0      0
#> 4ed3218615bbc1ed4da82d9806e4786b      0     0      0     0      0      0      0
#> 8a3367e36ed4e416d5a77a7946cf6fc5      0    16      0     0      2      0      0
#> 56f45eba0496ed29fb7171041a2a950d      0    17      0     0      0     16      0
#> 9d2782f19e60ec0cf72e1dc46b1d1766      0     0      0     0      0      0      0
#> 1559f8b64dc2bdc8cdc26ff8a5ffa4cb      0     0      0     0      0      0      0
#> 6d01f384a1cce140edea17dccd0f116d      0     0      0     0      0      4      0
#> 5950949b71082741302b712e97309b17      0     0      0     0      0      0      0
#> f6f238f95cc97995c9bf619b0c620b33      0     0      0     0      0      0      0
#> c7a112ffadf8009efd318625da44b18e      0    14      0     0      0      0      0
#> e79f4f30b0dacfcc50f3c1b12804ca3f      0     0      0     0      0      0      0
#> 297c825dc268f9a7b86b402edb6961b3      0     0      4     0      0      0      0
#> a1ed56ad4e86e0b4adc11baf2c4ab94b      0     2      0     0      0      0      0
#> 8afccb8f7f50c5603b471fb0c808c9e6      0     0      0     0      0      0      0
#> ce41c973bd4965a400989ec128890d55      0     6      0     0      0      0      0
#> 5d99ba689ad34e4393a905e8a725b72b      3     0      0     0      0      0      0
#> bfaa2d35e6dd13b9b32de811c6f92592      0     0      0     0      0      0      0
#> 0c544d5f3a49a010971a96d8c39e4905      0     0      0     0      0      0      0
#> bab3f08e1d75f42987012cae451639ff     79   215     30    42      0     79    107
#> fa97a81b8c33569e2aa1fc9a3137b435      0     7      0     0      0      0      0
#> a4d08886c3c6c62efd5678339a07828b      0     2      0     0      0      0      6
#> 11ed0972601847ec168473a6bc4eeba7      0     0      0     0      0      0      0
#> a6e4d88aa4b40402b1f94760d86b729c      0     0      0     0      0      0      0
#> 86deaff47d33946ede1b7c25c6bbe9f9      0     0      0     0      0      0      0
#> a1d3f1f8c151eaeb5cf12179861e0c87      0     0      0     0      0      0      0
#> 3980891c22ccbb387b65f9e437b2f94f      0     0      0     0      0      0      0
#> fc8ec321e0c2fc2a6de9abe98256b3f5     34   102      0     0      0     11      0
#> 84c11e6ffb29ebc2192f4b049b7dbddb      0     0      0     0      0      0      0
#> 82104585f0b2617c13eb69f0bbdf7d5a      0     5      0     0      0      0      0
#> 61cd4918428d3b36910b9cb2fc45946a      0     0      0     0      0      0      0
#> 0ab9a976935482bce5f8af4cc6685147      0     0      0     0      0      0      0
#> 87be3450216b45c53851fd173582a2d0      0     0      0     0      0      0      6
#> 93a51bb503fadf90121f518db980e158      0     0      0     0      0      0      2
#> b70cd843161996e075a853900e699d0f      0     0      0     6      0      0     10
#> 154b8be9eea158f93d27d314d8e0e2c8      0     0      0     0      0      0      7
#> 70691c2737cf7dc7c194de819f0e6cda      9    22      3     0      5      0      0
#> 2ae1d80ead92f5991f58da14d000acd6      0     0      0     0      0      0      0
#> 44fa299492152518fec34b9fa4554648      0     0      0     0      0      0      0
#> e15c8f9df87f363b0dbbf49de0bed249      0     6      0     0      0      0      0
#> 73689cad6eb4e86d9d0bb48c0da24463      0     0      0     0      0      0      0
#> 234e15e82d64f2e1d61802fa8ea1998a      0    28      0    21      0     24      0
#> 88c76ddb1b24acffe65d8128e2f48244      0     0      0     0      0      0      0
#> 5544139a16716cc6ef48b2c16c9613e9      0     0      0     0      0      0      0
#> 5b7fe68483e701059142cb531397b547      0     0      0     0      0      0      0
#> 0a0875a0acf4a28d1bc1c4d1d11aaa64      0    14      0     0      0      0      0
#> 6b47ec9d1c70f07d16834e9a4e54aa70      0     0      0     0      0      0      0
#> f8752b4824d239cfb19f3216cc34941b      0     0      0     0      0      0      0
#> ec8a0768cd2fb426d7e6bfc2c4f50e5b      0     0      0     0      0      0      0
#> 8a6b88f9515bf2dfac09e6a69045a4fe      0     0      0     0      0      0      0
#> f5151c67af19c8471743ef3e2914f6b1      0     7      0     0      0      0      0
#> c3842619fcbde52781c512c2561fca0b      0     0      0     0      0      0      0
#> b3b00faa9dee9fa6b3b33744f51da00e      0     0      0     0      0      0      0
#> c4f04e28c6267ac6c5b05d403e8f91c3      0     0      0     0      0      0      0
#> 7fc06f107a6ddf6486ef2d43fcf23626      0     0      0     0      0      0      0
#> 8d1dffcd377b338d2a0388637507b592      0     0      0     0      0      0      0
#> 15410d762889d51818986bfcdb5d8965      0     0      0     0      0      0      0
#> 48f68e44e258de63575b8bf6cf565d52      0    20      0     0      0      0      0
#> 4712befe3685661e0fc0968859c369d5      0     0      0     0      0     14      0
#> be9a3727a0422aea147100370d046fa8     21   122     10    13      0     45      6
#> 3c5942781761830f7b5406f551574424      8     0      0     0      0      0      0
#> 291d9f4d3d3c5cc3b0a1804770ed1e77      0     5      0     0      0      0      0
#> 784764c519a64dc9a147025a23222e93      0     0      0     0      0      0      4
#> b9e955ea5254dd68b54775aad2a67829      0     0      0     0      0      0      0
#> 8197f1594e8a17f3fc989fea64f54b8f      0     0      0     0      0      0      0
#> bf4a4f2383a16d2ecf27bf8e162b9ba4      0     0      0     0      0      0      0
#> 4e7d22d793d30bdedbd32a5b73b564b8      0    17      4     0      0     27     11
#> e88959b739c8cceb94b94b7887dd1df3      0     0      0     0      0      0      0
#> 5403445694339e2ae5ecac06d9f3ff47      0    19      0     0      6     12      0
#> ba3943fb075bbf63fb699638292826a1      0     0      0     0      0      0      0
#> 67b2d836db093cead5a8286142669613      0     0      0     0      0      0      0
#> 6c3a0ea164d3eb098dd2baac5519b38f      0     0      0     0      0      0      0
#> c643a776c36bcc00dfb5d790e1ea3dfe      0     0      0     0      0      0      0
#> d53cc81a9d57a64d194c876b76bd66e0      0     0      0     0      0      0      0
#> 7dc00d40971bb454b39fdd9b2dd9aa25      0     0      0     4      0      0      0
#> b4d0cec8543e940be1b5096370fdfb80      0     0      0     0      0      0      0
#> 6c54ceb94bc5e8fbd65224cffc60990a      0    26      0     5      0      0     16
#> b12ee9ca5a2320957213a114788cf2b7      0     0      0     0      0     10      0
#> 5e040c2397f959e7f0dfc6ec2855b812      0     0      0     5      0      0      0
#> 3105a594bba2604ae5ed19f5e0495d2e      0     0      0     0      0      0      0
#> d75c040bdaa1cca96cd0655c9af7ab28      0     0      0     0      0      0      0
#> 9b909edaebabc90b9c1a703626e6bd0c      0     0      0     0      0      0      0
#> a6593da81035cf121dbe46346a4ea960      0    28      0     0      0      9      0
#> 60e4c88400a750a8122b97476313a67e      0     0      0     0      0      0      0
#> 7854a78e3b881cee1c63ea0a50ba29ce      0     0      0     0      2      8      0
#> 1796734b56d45bdc997ccc09ab8731d1     25    16      0     0      0      0      0
#> 8243cd669ae7327f4326be78429ab05a     32     0      0     0      0      0      0
#> ebb93dc32c63402289114859d476cc5b      0     0      0     0      6      7      0
#> 8a01af5a531e851e2c50fc01131cda33     10    33      2     0      0     27      4
#> 03d2cc3f7b2bc352752d39e3738d0af9      0     0      0     0      0      0      0
#> ec378af6e32960b11b5f7391da4ae4fb      0     0      0     0      0      0      0
#> b407f2de868c064f65eb5ad8f82a8d0d      0     2      0     0      0      0      0
#> 73a95b6ba34fcd260989d848fe15d377      0     0      0     0      0      0      0
#> 3abb68202d15781590a2d1f9c5b37e15     11     0      0     0      0      0      0
#> 77923ef8e9ed01c9b6a8ac85eb028163      0    35      0     0      0      0      0
#> 9981b8fd20f400e052527f84f5db1b66      0     0      0     0      0      0     13
#> 597bdec49adf7675912f9fba4e262e70      0     0      0     0      0      0      0
#> 0fd3126118b13fb187bc55476a67c95a      0     0      0     0      0      0      0
#> 2a0c511ecc2513c42fe624e8e99d07da      0     0      0     0      0      0      0
#> d2912fea2dcf2e5b97e56c6a7d24295d      0     0      0     0      0      0      0
#> 2976af6aaf4261b6c274fcc3229a4d65      0     0      0     0      0      0      0
#> f097c45ce419da576aa111c98a4c7272      0     5      0     0      0      0      0
#> 3abc8260d34008f263e232ea25426c3f      0     0      0     0      0      0      0
#> d7e5b6432dc178aaeed9c81e0ec57afb      0     3      0     0      0      0      0
#> dfe71ce78ddfdbc7b453d17d628365f2     36    48      0     0      7     29      9
#> 289b18a86084273ce3b4038bf37e2840      0     0      0     0      0      2      0
#> 844a67032c52dc49019bf5dff248bf9c      0     0      0     0      0      0      0
#> 38cb3e3ce5ac5aed7dec50247ecfb267      0    12      0     0      0     12      0
#> 74e88e71e9929b1ddfcafccfab3b7bba      0     0      0     0      0      0      0
#> 66bc0737fc7fa1463ac7189ab7073998      0     0      0     0      0      0      9
#> b5a7aa588c9cae6ae8ca2f7d87ac53e8      0     0      0     0      0      9      0
#> db0540920457b21b170947aac9137771      5     0      0     0      0      0      0
#> ebbaf260510ab16ac6b8c2b9866d05b9      0     7      0     0      0     14     14
#> c1f461b67419d4a80a61db0867da9ed6     25    17      0     0      0      9      5
#> 0c6eba6452635eecfb865a1049299176     61     8      0     0      0      0      0
#> 7ecabb10f29a0c44eb53e9f9852fe22a      0     0      0     0      0      0      0
#> a33cb3037e7825920073b4dddaaa41e6      0     0      0     0      0      0      0
#> d625812c15c074c2e2b18099f97412c9      0     0      0     0      0     12      0
#> 17ef77f12e2a55717e910a5f6a625448      9    16      0     0      0      0      0
#> 7baf638fdd3f71166aafd22ecf410348      0     0      2     0      0      0      0
#> 9436cda211f60f722865cd23ef75d4ed      0     0      0     0      0      0      0
#> c8ab0b20e740eb6b35e68d58f69736e4     24    55     19    43      0     77      0
#> 5e43f2d8be71399fc8d6736754cb42ac      0     0      0     0      0      0      0
#> 19ad90ee1081a003406b02179ae63709      0     0      0     0      0      0      0
#> 0a25cdbd0660e68970c563c43a18f54a      0    44     32     0      0     21      0
#> 0aa6b210e902d8ab9d6f267407258ebb      0     0      0     0      0     58      0
#> a864c36d0d898cc5abd94d6d6f0bd5be      0    16      0     0     12      0      0
#> dea0dcb305ff0f27e3c36c7133cf933d      0     0      0     0      0      0      0
#> 463977baa69c7724fdc9407ad162714c      0     0      0     0      0      0      0
#> fcb72326c29342bcfb052953466eafc6     22    37      0     0      0     12     17
#> 5a554353f9ed728f0ab0631f29c8be2e      0     0      0     0      0      0      0
#> a9dd27b48e443415a743bf9902f515c4      0     0      0     0      0      0      0
#> 77a40bcf979c0099c52f19a3f1fe5217      0     0      0     0      0      0      0
#> b4007ac525a9b409d70426e24421cff0     27    17      0     0      0      9      0
#> 406ba18a175c4c79efb25f48a6beafe4      0     0      0     0      0      0      0
#> fe1c166755bc722cefd2e2a3b2b70672      0     0      0     0      0      0      0
#> 8afd274849ac07c0ca66f8f5f143df3c      0    24     13     0      0      0     23
#> 197301ed9d8728b5c061b99a17973161     51    74      0     0      0     21      0
#> 14779f9dc3505598ad3d799e4ee676c4      0     0      0     0      0      0      0
#> 4085c246640eee043d523ce2d799023a      0    27      0    26      0      0      0
#> bba09a403e67c1381de01c86511f1057      0     0      0     0      0      0     16
#> 392c59fee4e6b91d7654f2d3e3c65be8      0     0      0     0      0     24      0
#> d2b05ae73ad5ed296d7add9e2bb2cce4      0     0      0     0      0     15      0
#> 8306a9473dabfd5918ea869dbfcd1291      0    19      0     0      0     13      5
#> 06064c372bb624503c2905d13f7691f5     47   185     56    43      0     40    220
#> 3265692341e42dbd5cf4ed1938b906e3      0     0      0     0      0      0      0
#> da96a0797b9881a123795838b9f01140      0     0      0     0      0     14      0
#> 0bdddf10e577657becec88e42001fd16      0     0      0     0      0      0      0
#> 168999a950ebe17e6d2ea1f90e557d91      0     0      0     0      0      0      0
#> a266ffe9d38f626d106fa615da4618ef      0     0      0     0      0      0     15
#> 7706405d5bce7715f01d8b54188def42      0     0      0     0      0      0      0
#> 4b0cfc8899b466bd9992f7b22a5a06be      0     0      0     0      0      0      0
#> b81b2b4a7ca9eb2c23b170648dfaf22d     19    12      0     0      0     15      0
#> 75cf40872201f150cb5b191676f64b78     66    53      0    37      0     38     17
#> c8e5794c58de9326084a62d3282d257e      0     0      0     0      0      0      0
#> b0706a571f169c0bac7b4864052fcaa7      0     0      0     0      0      0      0
#> 0388f969e3527b2938acc676d6757bdb      0    30      0     0      0      0      0
#> 2dfe751390aedfbd9f62dbb77d1f10b4      0     0      0     0      0      0      0
#> 2548dc9697e1f9afe04a9012f40e2f0a      0     0      0     0      0      0      0
#> 2f019615d5141a4ab00b1aab98882a33     57    83      0     0      5     18      0
#> 629e4f0cf1826904ff74907b2507f58b      9     5      0     0      0      0      0
#> 774861f348cec8451102730c21a790b0      0    28      0     0      0      0      0
#> 11c0d6de28ab026a9ae17abe4efdd6bf     13     0      0     0      0      0      0
#> c11381dda0841c0d08cbfcef3a760aa3      0     0      0     0      0      0      0
#> 3c9619d4f9db4ac6ad7151a0c9d86ad3      0    16      0     0      0      0      0
#> 7e2c72204a8c45b1b0ac39e0e6047098     13    31      0     0      0     17    128
#> 8254e0f9f07709b62dd8483f3b25fad0     28    66      0     0      0     39      0
#> e78135f2886d427e678d12c446bf070e      0     0      0     0      0      0      0
#> abec938a8d53e7471f3f476d50694be9      0     0      0     0      0      0      0
#> 4963567d2979bb761291aa98e0b8595b     11    15      0     0      0     11      8
#> d4de9d1453102564f35262faa3952266      0     0      0     0      0      0      0
#> efa4a78c4b993e9d115841c8a02c795e      0     0      0     0      0      0      0
#> fcacd7424bf3543322ea3e9842d14e86      0     0      0     0      0      0      0
#> f941e15ea7f7560191d57abc4d2bedd8      0     0      0     0      0      0      0
#> ac294f5242de3c9b26e75502c48f6383      0     6      0     0      0      0      0
#> 7b6a8e3f8bfa1bcd1b5fc55db048d005      0     0      0     0      0     11      0
#> c9e776ab599fe75db286675b698c9067    115    55      0    32      0     89     34
#> ebe5cb945f6d6d849c07424591f80753      0     0      0     0      0      0      0
#> f116387f77802409103a559e4205df06      0    41      0     0     18     59      0
#> 93f9939b261cddc9a57071b42771c1d4      0     0      0     0      0      0      0
#> 971a5f48b3418449ca7cf3e190a5e831     22   103     17    26      0      9      0
#> 730de992ca06a95e7193a93146a97cf2      0    12      0    18      0      0      0
#> 8fa8b4ac69784c26901e036d13ab2332      0     0      0     0      0      0      0
#> c363c9efaa85023b488b732b97e47270      0     0      0     0      0      0      0
#> 2e2181cc893ca1e1e677e888bf8af441      0     0      0     0      0      0      0
#> 86fa662da6869634e7f90f85e9ee4933      0     0      0     0      0      0      0
#> 9a02177f2c235a6eb62984d5aab151ff      0     0      0     0      0      0      0
#> d6854be498e22b65da3c00572ee178f0      0     0      0     0      0     15      0
#> 015f364f150b230450da0a123e7bb8ab      0     0      0     0      0      0      0
#> af2c1b22758cd6653e8e841497d096e7      0     0      0     0      0      0      0
#> e7f2fab8c0e7c365768a43ef71406605      0     0      0     0      0      0     45
#> fb57cf192bce1876471f2f6181bf4b55      0    10      0     9      0      0      0
#> 6c740ce69c9f45bea684dc9e36054006      0     0      7     0      0      0      0
#> 3b27c7bf279535d8e39664d8dc1d0daf      0     0      0     0      0      0      0
#> ceea3c29ad012c2918fb6aa21a44e8cb      0     0      0     0      0      0      0
#> fa99ecdb27c1cc44fd866ee1da250c4c      0    33      0     0      0      0      0
#> 387cc83573db0e444e81c89a8a4f8f9c     26    78      0     0     17     35     40
#> c7a8d1af5f7f472671f700775de3e0a2     18     0      0     0      0     20      0
#> 98a6a84550e07ea252cfd8a0f51042d4      0    62      0     0      0     35      0
#> 08b809a68126c38db786cafe9710b013      0    24      0     0      0      0      0
#> 8db4cc0a5d2f03fce8efd2b8374ff359      0     0      0     0      0      0      0
#> 3e7c7d38423c93b43a1d9865cf0bc4a6     55    83      0     0     15     40      0
#> 6f02050760c030ed8cd1ed836697fa6b     10    32      0    11      0     34      0
#> 8fc1865cfcae32f72c1b1570299e2ef3     12    44      0     0      0     13     10
#> 38cc138fac95b006929361093d840a00      0    11      0     0      0     12      0
#> f5f8eb63e22362509e2c7c78271305f2      0    10      0     0      0      0      0
#> 0bb238aaee3ccb7dc8dfa86dbc1ba9f7      0     0      0     0      0      0      0
#> ca4546caf07edf5200cf79cea04aebcd      0     4      0     0      0      0      0
#> 893d12123eb7f623489a4f025a15ae7f      0     0      0     0      0      0      0
#> 192c5a817f48368339e519e0a0a907da      0     0      0     0      0      0      0
#> 3467185ca92b7bec2523214605adabf0      0     0      0     0      0      0      0
#> a2aae80d8fb922fe1074d243b59ce237      0    20      0     0      0     27     32
#> f11233f16d4177d38ab1c0bfc649e056      0     0      0     0      0      0      0
#> b219ca7793ab9f5e441f2a68bc76152c      0     0      0     0      0      0      0
#> e6e1848b51e705e55b3805703f088311      0     0      0     0      0      0      0
#> bd074b9bf365d11468bd76268d63b15b     31    36      0     0      0     30     20
#> f381b492d6651c52fd09ce9351e7eb2d      0     0      0     0      0      0      3
#> d9b31e1c30facc24fa778b1d65f5c457    275     0     16   109      0     17      0
#> 14a67c8c8cc31d203eac8109443fab5f      0     0      0     0      0      0      0
#> 026ad0ea53fb202ebb770acf65705d6e      0    10     17    14      0      0      0
#> 19ba9fea04969f975d6dc8d3da20d5e6     26     0      0     0      0      0      0
#> 0fc4426ed7dccc48ce27cfff8083433d     96   134      0    35      9     52   1012
#> 45ad90e179842705fa984962f131fd61      0     0      0     0      0      0      0
#> d284afab226d22fd0596221d17db23d3     39    32      0    20      0      0     37
#> 38194655ba17dbfd7fcc1ed42a57c169     32    64      8    14     11     74     12
#> e30be5268bbebe597424bdf171cc1607      0     0      0     0      0      0      0
#> cb7ccdb4956134b452a3cb78d1751c57      0     0      0     0      0      0      0
#> b585436677ac88dd4ac6239655f2a268      0     0      0     0      0      0      3
#> f5d770eff5c9b9b0c859db57d175eb34      0    10      3     0      0      0      0
#> 1b45ffbc88c788dcee2239e68901cf05      0     0      0     0      0      0      0
#> 3fd951cf7cc3a0740a1f20f8a9c6d2c2      0     0     11     0      0      0      0
#> 2027d3e90ccc4003f1954322e99db6a9      0    35      0     0      0      0      0
#> 6157cb92771830955413f808e83fc420      0    28      0     0      0     20      0
#> 6349506cac32056a214b04239bdc3a33     61   129      0    26      0     81     62
#> 2254a8d37cf537b1e35097d70b6ac91b     28   264     48    41     26    172     52
#> 38f3bc8ef8836574de0fcbeafd348b65      0     0      0     0      0      0     27
#> fec23c70ae0c70ecceea1487cd978ee5      0    29      0     0      0     37      0
#> 503d123a6c78c6fb280a6472bd61b51e      0     0      0     0      0      0      0
#> 9f33de8570c639e727465f1ee4d750f6      0     0      0     0      0      0      0
#> d68e1899fa65515135d9bd8d356b7e3d      0     0      0     0      0      0      0
#> e15b88ccac7a4bef64215e915f6e0d89      0    14      3     0      0     11      0
#> 767042b4b19d8bfaaa426dad7dbce378      0     0      0     0      0      0      0
#> 4a3c4abee502d6185abd750522566bf3      0     0      0     0      0      0      0
#> c254aa5c00b7f27d756f7bbe96b3d4ad      0     0      0     0      0     20      0
#> fe8d8af791eb6e5d4335bb7604c74e9d      4     0      0     0      2      0      0
#> d81800d7e31b69877ce1eab9858b8f3a      0    12      0     0      0     22      0
#> c55c7cb80c5700ff9aded76192543b5a      0     0      0     0      0      0      0
#> 91dee7b3b90bbf2e53b2472b9bbde0e9      0     0      0     0      0      0      0
#> 6b28df61cb5dee74b76b1a13aaa0027d      0     0      0     0      0      0      0
#> 9a830fa26205ad91a02b5bcaebb0221a     14    16      0     0      0     10      0
#> a54b8142d2a09d058b81c377e8e286b8      0     0      0     0      0      0      0
#> 69ef3a03a4d733dbffd99fdaa1d2e62a      0     0      0     0      0      0      0
#> 083cbb131dd305a88b7d338bc371c35f      0    12      3     0      0      5      0
#> ef1a92caffb4d7cae286fdc77b0a93a6      0     0      0     0      0      0      0
#> 2cddd5dca166486f50a0828b28001e1e      0     0      0     0      0      0      0
#> bb899d1dcfd89cb200a76cb8951904f8     33    35      0     0      0      0      0
#> 7817f6851e9df03a3d992cf915310dba      0     0      0     0      0      0      0
#> 54eaf6571de743a86ccb531bee8aad7e      0     0      0     0      0      0      0
#>                                  3.14RI 3.16RO 4.16RO 7.1RO 6.13RI
#> 374a5eb6496d14e8e9d6f3d8c7e34d9d      0     19      0     0      0
#> 4c47d5cf81df1ac6e2a0560029a81578      0      0      0     0      0
#> 8fb6741a6685fad0374f85ad3e16156c      0      0      0     0      0
#> 310f3de009c95de5c938b8e811dfc96f      0      0      0     0      0
#> 7d3b39e8c0f0fdcb81b0c898a98c38f9      0      0      0     0      0
#> fc211549300b0954dbb4a4bf57e8a605      0      0      0     0      0
#> 53fad3d8c52f4b5a023ff92aee9c0f1a      0     47      0     0      0
#> 5e37342d650cde6cfc5bbe2726b12e3e      0      0      0     0      0
#> 56025f74cc8e54d4a054442d786ef096      0      0      0     0      0
#> d9e51d98f1c86a1c426552d252c944e4      0      0      0     0      0
#> bc905dcc609da83d3620f1986e63a9f7      0      0      0     0      0
#> 90009df763e86e6eca7baa54f01551e6      0      0      0     0      0
#> 6bbd12b392461134369d835715cd2765      0      0      0     0      0
#> 83c88a72c2dab171c38aadd98aa2d389      0      0      0     0     10
#> 2b93708de9f1c24ff6140814a599a3a8      0      0      0     0      0
#> 671a53a2acbf0c8c7e240e446daf0000      0     36      0     0      0
#> bfe35644221421300c0c17e00e465c05     20     60      0    31     48
#> d12bdf20f24b92aa8ec94c6c3f60497f      0      0      0     0      0
#> 3bdb72ce060bd89f107c5c2bc95399ed      0      0      0     0      0
#> e6a82be4b75aac8b928c1dfbe9ee9287      0      0      0     0      8
#> e3a63ba0a0ec40fbf270a57c180e0d3c      0      0      0     0      0
#> 2877cf981b128a18f788ba5437c1afa6      0      0      0     0      0
#> c015bc50696fe3e3d94ace01d9dcd691      0      0      0     0      0
#> eb2eeddf829080f8d9f25c6be5d507cd      0      0      0     0      6
#> b1fcda6d8df4c25c8e6d7124c076b1fe      0      0      0     0      0
#> de8bb5c39fb121fd048dfe0d639e6e8a      0      0      0     0      0
#> a234ac03223c275785dedf6df3b2aeff      0      0      0     0      0
#> 9ea1641eeed1f51fc228c653076d1d08      0      0      0     0      0
#> 694b926113eb8291fea166446890c8fe      0      0      0     0      0
#> b2af163540d17007833dc819704afb28      0      0      0     0      0
#> 6fc0ba1f8ff8f9259c9d492e59771755      0      0      0     0      0
#> 7b8c1dba059eac7ec7f3cefa203f6b2d      0      0      0     0      0
#> 0b06028193ef1ce19f5339a8ba7d4448      0      0      0     0      0
#> 56b4ec05ef3bbded515ec3dd2d573197      0      0      0     0      0
#> fc9c84f8767611251f0eb93af4ca10db      0      0      0     0      0
#> d3616bc0d8925a452274c484c73406c3      0      0      0     0    103
#> 7986b8b93f3afce920bd19a5d2229c09      0      0      0     0      0
#> bc78ee7297b6912452be0865d113f7ad      0     63      0     0     19
#> c472668b46b0d0f575108aa4170dd81f      0      0      0     0      0
#> 2658a78c5ff06760389ce2e3a548b4a2      0      0      0     0      0
#> 600c9e81fd1d1ddd5a37a5aa0f9f796d      0      0      0     0      0
#> f8d0bb97f093d8e8630ea62ac8a1bf13      0     22      0     0      0
#> 6c685b83d63d6d56a9dcaba34215ae6b      0      0      0     0      0
#> e0e7e6208280a5846a8b66f00667f8ca      0    110      0     0      0
#> 6d1aa1c1f93c423a856a5dd3b62e5c70      0      0      0     0      0
#> e553e766ea0e08f1f9fa2debe7946af6      2      0      0     0      0
#> 0f1ef43f6b23a567aba80599766d70bc      0      0      0     0      0
#> 3c4ccaa45e8d148c9a17b5d0f3720cd2    103      0      0     0     91
#> d9df6265d046a558d87d397b8a177d54      0      0      0     0      0
#> 4bbad82a39655a26a2de61b88c5af7fa      0      0      0     0      0
#> f5a6865451a3b138f4915c5ab445a0a7      0      0      0     0     14
#> 4a86fc550fa1ac37aa480d37edea875f      0      0      0     0      0
#> 02fc21e72f66adea757b355edea4a369      0      0      0     0      0
#> 25cf10f4f15dab00925cef52b9f08442      0      0      0     0      0
#> fa8c8b98d197c252368c160b43cefb98      0      0      0     0      0
#> c0b07d29aed1a875641aede53bc0a20c      0      0      0     0      0
#> 341433c1c8cc65ce47dcf32179d08fd8      0      0      0     0      0
#> 8bf175f329d16e4aa732cf2b32279df3      0      0      0     0      0
#> 062f7552b747dccab9586dff7275b5d9      0      0      0     0      0
#> 059e124ad6ba8f98c3d21db28a075ea8     81    252      0    53     81
#> 3d7a81faea3f783158352c9f739c8d17      0    138      0     0      0
#> 90be2fb2c82a98002d33ed7531e877e3      0      0      0     0      0
#> b0e263cdd1c189fc1757d8189db14be2      0     92     44     0    130
#> b962ce332fea338b05293d1a5a29cbca      0     30      0     0      0
#> 7a310a8f8e4bcbfb717263bd5376994d      0      0      0     0      0
#> 57136cbdae8fe066c2bb4a661003719b      0      0      0     0      0
#> 69fe46cab1cc9407470a2397bfc84103      0      0      0     0      0
#> b88efba97658a6936bb053985d026ace      0      0      0     0      0
#> 4d691fbbe2dbfb2815649fb72c757977      0      0      0     0     16
#> 49192cc8a222d97b965b839ff5d0f333      0      0      0     0      0
#> 859873fcd088c294a92a13340744eb6f      0      0      0     0      0
#> 0a8c44ca4d744e313e6dc5b239fca562      0      0      0     0      0
#> 09f94345824c82e8fe0d824cc4478fb5      0      0      0     0      0
#> 9dfbdb741fb0dc25d2206b8680bbf579      0      0      0     0     38
#> dc2bea463c720874cccd189ccdcf0e3e      0      0      0     0      0
#> 6273eb52991d528cf5281215f79b40a3      0      0      0     0      0
#> dc4bcbe74986e44ec3a040b866dfda9e      0      0      0     0      0
#> 96ec4781ee61f4083792f252f1b6cd3a      0      0      0     0      0
#> 9657d52acca51e701e7a3b7cedcfcc12      0      0      0     0      0
#> f65f183293d4710041057ab98f2c94ed      0      0      0     0      0
#> 8f6b3783b31a2dd8640de3df92dbb9c4      0      0      0     0      0
#> 0d90eb4842c9ad9126ed33c538c768ac      0      0      0     0      0
#> 8d3ee46aaaf7728594014be61c649347      0      0      0     0      0
#> e011a95fda1f1f2bafc4e5b83b52d989      0      0      0     0      0
#> ef8ea65ba4c2e2da3ba9556a87d2385d      0      0      0     0      0
#> dbbb92eb3a8a3feaf96092d46f2b6a45      0      0      0     0      0
#> 06cb03cddd8a12dbfaf27e8984b395c4     96    101      0     0    130
#> cdf134bb501c54d997945dabd44e35ff      0      0      0     0      0
#> f3c4561cc01a3b45ee3f134d709c6b94      0      0      0     0      0
#> b121442fd1eefe78ce4f4602aac1842e      0      0      0     0      0
#> ba8707aa88a2a62939d9338189b2a733      0      0      0     0      0
#> ecf7289628d49f9244469337ee1622cb      0     48     19    27     25
#> 9f06b400b93f859b11874b0bd7f09295      0      0      0     0      0
#> f775be5cf7c0dd2177006a44913503d5      0      0      0     0      0
#> e87ceb10f2becb2a6c4cf0692b1923ee     45     41      0     0      0
#> 036574bb8e05f3b9a11ab7a4090d6950      0      0      0     0      0
#> 0267ba9db97e00ae01fc2cccbb6263e0      0      0      0     0      0
#> be4f701c1dbfbda16c5150af29430141      0      0      0     0     64
#> 655a66e6566798ade91f73b2942c0b63      0      0      0     0      0
#> cd79f2696509d0ad906f5addba264079      0      0      0     0      0
#> 07ea3284bda2bf7daf1251ad2744e4c4      0      0      0     0      0
#> e3870da9de4a2b18a7ae9e9cae208e3e      0      0      0     0      0
#> e00bd59c1ea0ca9b1cc5f976801d4e1b      0      0      0     0      0
#> b39b381d0f533c9923c3b898ffab5614      0      0      0     0      0
#> e2959f0552fa6df5707d9ec8b57816d6      0     39      0     0      0
#> 53157acd1642f973842d01e1cd287421      0      0      0     0      0
#> 5991790d1c0aa8b73b2d7f906449a1f9      0      0      0     0      0
#> adee3552c211d9e33aaa28e82532bdea      0      0      0     0      0
#> 0089040d041888e3ae4aa4e0010a0e36      0      0      0     0      0
#> 8e876e13d051992c4c8616cf2db7d73a      0      0      0     0      0
#> 179abbc5c1f1549ffd0d7ae6963e3783      0      0      0     0      0
#> 340d095f9a9fb445582b60e1282d985d      0      0      0     0      0
#> 3d7e5432d8eaa3fa56ab7138605ea43d      0      0      0     0      0
#> d4b233740db907d69035915e8657fa8d      0      0      0     0      0
#> b22119df3b1de76745abaac585b646f8      0      0      0     0      0
#> 64761a24fc66e7c93145217dd153b1ee      0      0      0     0      0
#> eb6d7f9e0601207fdcfc4a88947d513c      0      0      0     0      0
#> 0d20271e3e2dbc68e1de868599a9a1ae      0      0      0     0      0
#> 369ea22772cb774950c830d0af7e8235      0      0      0     0      0
#> cae63066cfe088a49fcac50f119c9e2d      0      0      0     0      0
#> d9c57ec20ed90e2cebf6cacc4f9a5511      0     52      0     0      0
#> f5828cc3078b4c219a44c67e96e5d455      0      0      0     0      0
#> 73847ac2778d4cba27fe9758ca129d6c      0      0      0     0      0
#> d78003f13712da48d0bb36fd45116b91      0      0      0     0      0
#> eb40853d986ad72c82658d939f291ff2      0      0      0     0      0
#> 5d675a3518222fc99c8cba34feaeac81     50     74      0     0      0
#> 3ba390ae11c4e76375af082eefe7e41f      0      0      0     0      0
#> 3aa66daa4e560cf62c1833697bcd3837      0      0      0     0      0
#> 6f4388510be5ea70356f2257d8c57ddf      0      0      0     0      0
#> b33390d034148edf07d8a742778ddcbd      0      0      0     0      0
#> 09225ec7044394bbd3b6811a6799e481      0      0      0     0      0
#> fe6fbf9b86e0f487911a856009fa4a1e      0      0      0     0      0
#> 0aee6329e0b5b81449980151bf40205a      0      0      0     0      0
#> 0ef66a78236f90916afc2305f3ff6c80      0      0      0     0      0
#> fb1eb41359005c2e3e308fb668f07220      0      0      0     0      0
#> cb59c57800227006aead644df42b8986      0      0      0     0      0
#> 138da9dc244c3292536cacb8df79467b      0      0      0     0      0
#> 72f92361ec0696cb8df985b93d277a25      0      0      0     0     25
#> 4d5ab121763879fd59e8f00d35473540      0      0      0     0      0
#> 20032b2b586fee22805eccebe3059170      0     70      0     0      0
#> 710324ccc8cafa6fc9b38baf8299d8aa      0      0      0     0     27
#> 523e6b60340f4a215877b1c336822e50      0      0      0    28      0
#> 4501f98c83484fd902be11383c8d032b      0      0      0     0      0
#> d4cbf23757aad57c91b058cbb05566e4      0      0      0     0      0
#> de7964d5c22112693ed5796bb9f04b0b      0      0      0     0      0
#> f4a283be42989b920b83ea73ec72b25b      0      0      0     0      0
#> 2f3111649d0850ab86799f8bbadc28f8      0      0      0     0      0
#> b8acdd45cf33ab75fe4c0851a1787dc5      0      0      0     0     43
#> 5bb7fec39744a6d4e4d61596dabc6886      0      0      0     0      0
#> e8eddfb37c4215c86813163a27f202a5      0      0      0     0      0
#> 826e188c911c0846fd927300127c80c0      0     18      0     0     18
#> ff1fe599c44074ee6bdf919a119df93c      0      0      0     0      0
#> aff00055d3876f624d67fb80bdb59934      0      0      0     0      0
#> 1beded194f3ae4b57715b45e50adef46      0      0      0     0      0
#> 908d91feba7c1f1d8cc8a7145173d942      9      0      0     0      0
#> e7095743695f6d4337d01ceb9e29e224      0      0      0     0      0
#> 4cd420962855afdf0bc69bb268b846c7      0    137      0    76     79
#> 430ba5bea0fc08f65b69497bfbeb5c4b      0      0      0     0      0
#> e12cb2f0d2ce4401413176007dcd4158      0     33      0     0      0
#> ec615213f6c7febb65c9f2e45e438111      0      0      0     0      3
#> 412dd52cd54d0244ba8fc42ad140acda      0      0      0     0      0
#> 89192d9f779341d313f52aad0dc6a06f      0      0      0     0      0
#> ed543f92c61a06fc4e703655ce910a2f      0      0      0     0      0
#> fea1a70b0ef0f72e3bf27d8b295331a8     18      0      0     0     23
#> 6afe805b36f5aa29e7f2096aac42ace1      0      0      0     0      0
#> d9cc5a46831f24ec9192ae82dbc38b0a      0      0      0     0      0
#> 0d3eb38eb502c4ad13bd5a7060295411      0      0      0     0      0
#> 69a0ff5538a5490e291c75ce70325066      0      0      0     0      0
#> f91d545c815effa3c5afa83a8ad3397d      0      0      0     0      0
#> e6b032dea8e13c34af958fc48af09b58      0      0      0     0      0
#> 84a1b89045ab6285328f366224f3ec87      0      0      0     0      0
#> 4a6f0076e4ac71474b3c364294ad3b33      0      0      0     0     26
#> 7bdd39e4e3ccc7000862b487510e67d6      0      0      0     0      0
#> 998aedf12cdf054f5e504a2e21d20498      0      0      0     0      0
#> 567c2640108006c8eca081ad4180e346      0      0      0     0      0
#> e50a6625b9176b6856764b213013150f      0      0      0     0      0
#> 8d5a5ab7c784b0325a4a86a16fbba80e      0      0      0     0      0
#> f58ae28a309fa3d5898d2eb46d9eb9e0      0      0      0     0      0
#> 7ac65a8fd21111b68e7c2d6205112a83      0      0      0     0      0
#> d7add64c8700f5a8cadf40c8c7f5ff67      0      0      0     0      0
#> 6bad9d981d3cf0f90262547180dc6371      0      0      0     0      0
#> 149c0e63c83113e1246c0e041de67d9c      0      0      0     0      0
#> 6f196fae54ea90b5c11627c503c98940      0      0      0     0     21
#> b67c995d842e5b4b92fa0ce052603ef1      0      0      0     0      0
#> f44ffae67c7297e2c696b75fbd4e0f4c      0      0      0     0      0
#> 33997912cd952feaec18fbd11b31acac      0     22      0     0      0
#> 0a035423c4bf71a17566549f9b69448d      0      0      0     0      0
#> f65a6f1f1e5a4bc5c733ba364c0096d3      0      0      0     0      0
#> f15389a9f245465c33db3de72aeb17ba      0      0      0     0      0
#> cd2b676d3c0785a60359c0d416701ee4      0      0      0     0      0
#> 70ba2a378ebcd9a8f57b024a561cc6b1      0      0      0     0      0
#> 0273707384afdd41e9d022b539094211      0      0      0     0     50
#> f9ac98e876ecdb504202388d7ccdda5f      0      0      0     0      0
#> bdd9b87b782b8f935ecac62a9824e6a1      0      0      0     0      0
#> 1f12b97495668f1fc378a349546d1806      0      0      0     0     34
#> 83f1b4d776f6e08a6ba6f86634a1e5a9      0      0      0     0      0
#> 49f614420449f6c7ecc8c9f4e137ad05      0      0      0     0      0
#> aba5aa00a1b262f27517eca64ee3333f      0      0      0     0      0
#> 8971980484eeabeef00145b93b699e30      0      0      0     0      0
#> c1347ae54b0fc1dd385a1ccf3b8061ce      0      0      0     0      0
#> e9fe19783ebcc7e43a3a07ca8225db62      0      0      0     0      0
#> 7c523a5f0a7c736efd82afef84137901      0     22      0     0      0
#> e13c07ec336c89db5183cc6e9e57ee6e      0      0      0     0      0
#> a017ce431ccadea609b501d46e64e55b      0      0      0     0      0
#> fadf38d03d4b944f5782a759591ec3ae      0      0     28     0      0
#> 41bfda1c7f32205246e73a1e3e3ec4f4      0      0      0     0      0
#> 9e4f3e697b70758a7a996f3dd4533650      0      0      0     0      0
#> 39c4f6caa7510ed108f3ebcf72f79edf      0      0      0     0      0
#> 4aa083fe7c1b0639dbc48ec649e68a31      0      0      0     0      0
#> fd2ac2e14b72b0e22a052e3ab4580471      0      9      0     0      0
#> 4ee10294799bcfeab0c8ee567d5d9608      0      0      0     0      0
#> 46869f12c31de7a05781e930c25a1ccd      0      0      0     0      0
#> 5e2548d6a970645d0ca943a72d207489      0      0      0     0     23
#> 0538e1beed69f098c3f602130ee11a1c      0      0      0     0      0
#> 991c563278c4e8c91ffebce3560fe51b      0      0      0     0      0
#> ce05f58164e61517198f7ab0e3304cb5     11     10      0     0     17
#> fdbd0fed21bae91c05b2da268a025a89      0      0      0     0      0
#> 7412d11595d77f26792071cd28a26760      0      0      0     0      0
#> 4b087e13957eb38128f63e1c99df23c3      0      0      0     0      0
#> 310a6f0e6841db2ce75005cac942a972      0      0      0     0      0
#> 3df0fada803ac582802b47d4f4efe0c1      0      0      0     0      0
#> 3cd4e34c02ffdd1446df65bb54278c57      0      0      0     0      0
#> b14a1dea93554257791f16c2e94906b0      0      0      0     0      0
#> 902a4435fc3af101e8a9b5f51550fa4d      0      0      0     0      0
#> 1995185b31349529f58964c2e5a3becd      0      0      0     0      0
#> 69b247c565afdeacee0240482b1ce536      0      8      0     0      0
#> 95ff1b9528acca018e9c4697a8625504     38      0      0     0      0
#> 132d26d3ce6ac2ab25ccb4dd171879d7      0      0      0     0      0
#> bf1a6e4f94ee8ba60cf79089139641b6      0      0      0     0      0
#> 704189aa09f64b992b3a02ee001dc706      0      0      0     0      0
#> dea5ea27839ee67a31d18ba0bf0f7f65      0      0      0     0      0
#> 3840a4fff0250232f8570cae7c3514ba      0      0      0     0      0
#> 31ecabc032997aaf41dd0d8def3ec93f      0      0      0     0      0
#> 024a0e0348ed9abe9ae1c48b1b7ba09e      0      0      0     0      0
#> fc1d9940419420113ce3fbfacc8d703a      0    145      0     0      0
#> 0770fbf28cd307a08bf557af7ffecfd3      0      0      0     0      0
#> 883cbe8e4d1d98f468355e8ea974d222      0      0      0     0      0
#> cb6d6561383fb68f43f92a88a607749d      0      0      0     0      0
#> e07e0f1a9fbc696a9fe4f9b7739ae45b      0      0      0     0      0
#> 73d8d8f999acfa13cba2277f7526897e      0      0      0     0      0
#> 5218274c354d41965abbf4b0751d6f3a      0      0      0     0      0
#> 5c31570f706752fa4f7851aca5b5c291      0      0      0     0      0
#> 6a5e71c8b86a52ab1dd5e61164103940      0     35      0     0     22
#> 25641b0b803c7d8493ebc02cb104bf3d      0      0      0     0      0
#> 8844d70dd1fdb96f2a7faa65bccf78c6      0      0      0     0      0
#> 371b855afa8bed030e86ec78241b01a3      7      5      0     0     12
#> d0c180a47378d414a1d088ac14db236c      0      0      0     0      0
#> a072a588b3940fab3b227e4901c40f7f      0      0      0     0      0
#> 977a3523051e6b23da6221dbde795b34      0      0      0     0      5
#> 1c3adc87d2e93b08e576953508eb680d      0      0      0     0      0
#> 74eb0ddf379b00c5032b18fb9fbf1f32      6     20      6     0     44
#> c335b34a15db7d532f511e6c952abe62      0      0      0     0      9
#> 0c78285397fb853a0e39069213510c54      0      0      0     0      0
#> 779dec6de9c6a8b7df28fc8ec3c75bee      0      0      0     0      0
#> 5526624b473bb1f0541dc6fb95c98678      0      0      0     0      0
#> 649e7fa7b9ce5c5472d34675ed643639      0     15      0     0      0
#> a9c5addeea542953d1bae4d30d599b02      6     28      0     0      0
#> 68291fb3b2558d438da4810961c8d53a      0      0      0     0      0
#> 9c5d92bfefc7a190b93d5b9a6a77697b      6     25      0     0     23
#> 3ddab42d701754625cab0c236c5b6e6d      0      0      0     0      4
#> 287e840f9ecaef80d7a4762c8d6bf401      0      0      0     4      0
#> 885f411fa54c5576148678f665bb053f      0      0      0     0      0
#> b5709a57df9cbc8a741bb78a36f29c06     50    142      0     0    121
#> 26698cc60893017ef0b1ab173688d3de      0     28      0     0     15
#> db64f0028b6518839ab7129c4b7ef72d     58    181      3     5    244
#> f46fb43262da299ec7f0084298cafa50      0     72      0     0      0
#> c45e7087de669be270881fc8ace0fad8     59     44      0     0     81
#> 854a20d5019388e81b1e78419eb017fd      0      0      0     0      0
#> 9a2325ea015a8c5944047e3349cc1565     15     29      0     0     48
#> c5305ae3c6adc90814d70ab48f4f071a      0      0      0     0      0
#> e230123ec3ce8654942ee1a000a80010     43     67      0     0     59
#> 14e3c60cb885509667670151bc4a1df0      0      0      0     0     31
#> 7fd7b7bece92bbbfcceb2dcf1e835a8e      0      0      0     0      0
#> c098d1fb769dd975c2332dfb70da489e      0      0      0     0      0
#> 95dcdf1eb06b252e67dcdbdaffc505e8      0      0      0     0      0
#> c5399d0f9eda8f0815fd28a257952496      0      0      0     0     72
#> c7388cc314e992819b8f97ea35280982      0     40      0     0      0
#> ae8facc45d372f7108bd9499f74d1bfc      0      0      0     0      0
#> 16caf6f5bc653fb6ba668eb821076a21      0     54      0     0      0
#> e7e283b72ee15bcf0877c5eeb6137ebc      0      0      0     0      0
#> 1b3a779b384c73b33207dccd14b14fe6      0      0      0     0      0
#> c8d1afdf7b395ffc22193c5dc2fd5dc3      0     12      0     0      0
#> 0e3ffc2a68378086e2480831acd1a025      0      0      0     0      0
#> 0ddf7339b608a340252b33f63fcc0019     16     18      0     0     32
#> 5fd4b2f7dee3fb0bde6f856ae266dbcb      0      0      0     0      0
#> e2f4786bc66410c07bd8b96b5d0e2d92      0      0      0     0      0
#> c30aae8f5216c9536492d0beae74950f     82    103     21     0     69
#> 9087ecee0806a80861788a17dd7e6809      0     16      0     0     17
#> c29bc08418318360592b2b2717376a0a      0      0      0     0      0
#> f791c7db1a78d35e84803b5e03c01965      0     19      0     0      0
#> 44b44fc4055a58c01f6d1e9de5e0d1d9     87     93     68    78      0
#> 342df9401103d274851d3bd80704508d      0      0      0     0      0
#> 0590e6672bb2c55c0ff54999a2b6dcfe      0      0      0     0      0
#> d4e9147608bf883fdbbf9d23b17028d1      0      0      0     0      0
#> db24794eb8667aed224401386223f671     14     21      0     0    136
#> fe41e6fd9770030a21418b288930710f      0      0      0     0      0
#> 699f18f3b155bff97de941f2cfbbbadc      0     26      0     0     48
#> 13870cf319997750e3b96bbb0b2f1aea      0      0      0     0      9
#> 7ca4fb1b1740406241acb862445cf311      0     56      0     0     56
#> 9c0569fb95451d009fd3754303eaed12      0      0      0     0      0
#> 1ac26a3c12288d8fbadd6aa61366b1b7      0      0      0     0      0
#> 89f3ba99552605501b8acb816d288afe      0      0      0     0      0
#> b4051b5227746f04c1feb9d8979299ff      0      0      0     0      0
#> 087d358aaf77f81aaf50a623828865bb      0      0      0     0      0
#> e451c4b3ff1c3dabeee296529b3a7a07      0     12      0     0      0
#> 012871784f77c9dca7d446ddef2048d1      0      0      0     0      0
#> 95c6ce2788783bf2f83343deba23b8bd      0      0      0     0      0
#> 213227f244d97553a9a14c4e3fd04c98      0      0      0     0      0
#> c0e5e55c5cfbd576e7c085f1e4b69638      0      0      0     0      0
#> 51e98451d7e3f46127481ec5490f1c91      0      0      0     0      0
#> c75285a6c61d8fce20d823769c25c12a      0      0      0     0      0
#> 177966f7619eb24fdbce2266f5bed39a      0      0      0     0      0
#> b89347f0534679f890628364d2c7b4fe      0      0      0     0     24
#> 8dff1cb12c0beb5bc7b2546480a6a6c7      0      0      4     0      0
#> 84f2777cfcb0fbed6281917fcc765022      0      0      0     0     21
#> a5e4d106a2ab947984fdef187996a0b9      0      0      0     0      0
#> b6cbf6077b64acbb718971a1ffc97a27      0      0      0     0      0
#> 27e29495fda424728fd6a80b83a16772      0      0      0    12     20
#> 63c44ab3725dd0b803d43a29db9fd7fd      0     42      0     0      0
#> 9c931ab5bf18ef9b10996d59d292e6f5     52      0      0     0     36
#> 12fd8efb82c3d4909543f6db9853134e      0     17      0     0      2
#> a8a60c20f7ab7d0f2dadd9f8858c19be      0      0      0     0      0
#> d57e09eb1645f1abba6074d3bd7d4c52      0      0      0     0      0
#> 7e0f2b982ca2d35b6ac2cf251d78d3c4      0      0      0     0      0
#> 3e2e8e938a68a0318bcf7145ad5ed4a9      6      0      0     0      0
#> f2117b49bf57d9451bbef6de5b762483      0      0      0     0     15
#> be24dfe003725a0b0f28f4136106faf8      0     50      0     0      0
#> e63301d7d0af55d5c43db051dc51522f      0      0      0     0      0
#> a6e2e5501507562f2174f65e606b2d43     11     18      0     0      7
#> a06aa93cb50b9a39941a11c7242dde96      0      5      0     0      0
#> 6a0a8abb36aa3bfdee14edbe97e9d781      0      0      0     0      0
#> 4d3def5ffad0c392381591440bfc3668      0     19      0     0      0
#> ce437b81e7cf5ce986f7240c26464b0a      0      0      0     0      5
#> e57788f9d2fae76ed9ca566ad4d475a8      0      0      0     0      0
#> 8c525d0f63137bd068b75c57f5d04c09      0      0      0     0      0
#> 6d670aa40dafa152c2d59cf2be9bef74      0      0      0     0      0
#> e08f8e9e9ffb10a702f2ef617f70eecf      0      0      0     0      0
#> 7def1c807ca0767ccdfe5a6001dc51f1      0      7      0     0      0
#> 3b9a17275147adf1dc8c42eb00a06a36      0      0      0     0      0
#> 24be0b7e4c85ac53b09fddbca916ab96      0      0      0     0      0
#> f521a47d3bb4acfec6b261f880b4eea8      0      0      0     0      0
#> 91c24c9b6ffbedcc257d164176060fe0      0      0      0     0      0
#> 9b45f210647c5fe71b005eb9e25299a5      0      0      0     0      0
#> 974077d25fa6b32a7af22670b142a8f1      0      0      0     0      0
#> 203aa2535e54a11e4c48af4ab946686b      0      0      0     0      0
#> 9064afe569448558eb08e5bf5e458bce      0      0      0     0      0
#> fad432eeeff7f1d923e6d60aafe18dcd      0      0      0     0      0
#> a29db0e5f4a04c5c6937745ccccb167b      6      0      0     0     10
#> f179c10f8e86873c90ceb14b974a114e      0      0      0     0      0
#> a23e8588a2240d7af6f65519158c902c      0      3      0     0      0
#> 52034e0da6d5e7962369437956984e18      0      0      0     0      0
#> 522d20c27cf4237b60f36ffdca056c19      0      0      0     0     17
#> 957615aefb421e6ad4728e610b2bdfa6      0      0      0     0      4
#> 5ccc74e70b8a82a6eebea3385738ece3      0      0      0     0      0
#> e5924fd11a48419c6ec0e2a7cb4e43ee      0      0      0     0      0
#> eef7be4b2531cb7e909b58e91af5640d      0      0      0     0      0
#> 4416cbc2026260fa4e2ed0a8feff1443      0      0      0     0      0
#> 03f6091f5aacb69f06fc925acd1c9506      0      0      0     0      0
#> f655e5ddf5074d5e71bfa583926dec30      4     11      0     0     17
#> 9f733d39e588c7310499dd444cdd7812      0      0      0     0      0
#> 2b413de5e05a1385bab229f50879628f      0      4      0     0      5
#> 178b596f2fe5cb6cefc5513917988d84      0      0      0     0      0
#> aa94b014c726cf40983045bd606241c1      0      0      0     0      0
#> 55310b3aac4868c4e33b0a9492d7609a      0      0      0     0      0
#> c5d8d828df0e1baf5c2db05a6b8c4476      0      0      0     0      0
#> 6ceb26f60606233a7728d3457af84380      0      0      0     0      0
#> 1ff83868d1c9ecf52d9b776e0b3fd896      0      0      0     0      0
#> 3edfccfb9f81b340854e84ee5102a4c1      0      4      0     0      0
#> 963dc552c551872800f2fcbc5663e5fd      0      4      0     0      0
#> 3fda8c20e5da9fbc3fb7c5c864559019      0      0      0     0      0
#> 6fd88fa7e883bdead3c5ba89da853e8d     10     13      0     0      5
#> 35826fe826b7a9f202e8bef7b39d9a47      0      0      0     0      0
#> cc514f614ba146f489bf62b5a56074c5     26     25      0     0    124
#> 26260de3858081ce55e6f220b3ba25df      0      0      0     0      0
#> d6727a85782ffd2a069dccc23adaf521      0     70      0     0    144
#> eb1a74c2ca8831ef945775de0827ddda      0      0      0     0      0
#> 4ed3218615bbc1ed4da82d9806e4786b      0      0      0     0      0
#> 8a3367e36ed4e416d5a77a7946cf6fc5      0      0      0     0      0
#> 56f45eba0496ed29fb7171041a2a950d      0     10      0     4      0
#> 9d2782f19e60ec0cf72e1dc46b1d1766      0      0      0     0      0
#> 1559f8b64dc2bdc8cdc26ff8a5ffa4cb      0      0      0     0      0
#> 6d01f384a1cce140edea17dccd0f116d      0      2      0     0      0
#> 5950949b71082741302b712e97309b17      0      0      0     0      0
#> f6f238f95cc97995c9bf619b0c620b33      0      0      0     0      0
#> c7a112ffadf8009efd318625da44b18e      0      0      0     0      0
#> e79f4f30b0dacfcc50f3c1b12804ca3f      0     13      0     0      0
#> 297c825dc268f9a7b86b402edb6961b3      0      0      0     0      0
#> a1ed56ad4e86e0b4adc11baf2c4ab94b      0      0      0     0      0
#> 8afccb8f7f50c5603b471fb0c808c9e6      0      0      0     0      0
#> ce41c973bd4965a400989ec128890d55      3      5      0     0      0
#> 5d99ba689ad34e4393a905e8a725b72b      0      6      0     0      0
#> bfaa2d35e6dd13b9b32de811c6f92592      0      0      0     0      0
#> 0c544d5f3a49a010971a96d8c39e4905      0      0      0     0      0
#> bab3f08e1d75f42987012cae451639ff     82    109      0     0     50
#> fa97a81b8c33569e2aa1fc9a3137b435      0      0      0     0      0
#> a4d08886c3c6c62efd5678339a07828b      0      0      0     0      0
#> 11ed0972601847ec168473a6bc4eeba7      0      0      0     0      0
#> a6e4d88aa4b40402b1f94760d86b729c      0      0      0     0      0
#> 86deaff47d33946ede1b7c25c6bbe9f9      0      0      0     0      0
#> a1d3f1f8c151eaeb5cf12179861e0c87      0      0      0     0      0
#> 3980891c22ccbb387b65f9e437b2f94f      0      0      0     0      0
#> fc8ec321e0c2fc2a6de9abe98256b3f5     27      0      0    13    127
#> 84c11e6ffb29ebc2192f4b049b7dbddb      0      0      0     0      0
#> 82104585f0b2617c13eb69f0bbdf7d5a      2     23      0     0      0
#> 61cd4918428d3b36910b9cb2fc45946a      0      0      0     0      0
#> 0ab9a976935482bce5f8af4cc6685147      0      0      0     0      0
#> 87be3450216b45c53851fd173582a2d0      0      0      0     0      0
#> 93a51bb503fadf90121f518db980e158      0      0      0     0      0
#> b70cd843161996e075a853900e699d0f      0      0      0     0      0
#> 154b8be9eea158f93d27d314d8e0e2c8      6      0      0     0      0
#> 70691c2737cf7dc7c194de819f0e6cda     21      0      2     0     19
#> 2ae1d80ead92f5991f58da14d000acd6      0      0      0     0      0
#> 44fa299492152518fec34b9fa4554648      0      0      0     0      0
#> e15c8f9df87f363b0dbbf49de0bed249      0      0      0     0      0
#> 73689cad6eb4e86d9d0bb48c0da24463      0     14      0     0      0
#> 234e15e82d64f2e1d61802fa8ea1998a     57      0      0     0     48
#> 88c76ddb1b24acffe65d8128e2f48244     15      0      0     0      0
#> 5544139a16716cc6ef48b2c16c9613e9      0      0      0     0      0
#> 5b7fe68483e701059142cb531397b547      0      0      0     0      0
#> 0a0875a0acf4a28d1bc1c4d1d11aaa64      0      0      0     0      0
#> 6b47ec9d1c70f07d16834e9a4e54aa70      0      0      0     0      0
#> f8752b4824d239cfb19f3216cc34941b      0      0      0     0      0
#> ec8a0768cd2fb426d7e6bfc2c4f50e5b      0      0      0     0      0
#> 8a6b88f9515bf2dfac09e6a69045a4fe      0      0      0     0      0
#> f5151c67af19c8471743ef3e2914f6b1      0      0      0     0      0
#> c3842619fcbde52781c512c2561fca0b      3     52      9     0      0
#> b3b00faa9dee9fa6b3b33744f51da00e      0      0      0     0      0
#> c4f04e28c6267ac6c5b05d403e8f91c3      0      0      0     0      0
#> 7fc06f107a6ddf6486ef2d43fcf23626      0      0      0     0      0
#> 8d1dffcd377b338d2a0388637507b592      0      0      0     0      2
#> 15410d762889d51818986bfcdb5d8965      0      0      0     0      0
#> 48f68e44e258de63575b8bf6cf565d52      0      0      0     0      0
#> 4712befe3685661e0fc0968859c369d5      0      0      0     0      0
#> be9a3727a0422aea147100370d046fa8     19     55      0     0     47
#> 3c5942781761830f7b5406f551574424      0      0      0     0      0
#> 291d9f4d3d3c5cc3b0a1804770ed1e77      0      0      0     0      0
#> 784764c519a64dc9a147025a23222e93      0      3      0     0      0
#> b9e955ea5254dd68b54775aad2a67829      0      0      0     0      0
#> 8197f1594e8a17f3fc989fea64f54b8f      0      0      0     0      0
#> bf4a4f2383a16d2ecf27bf8e162b9ba4      0      0      0     0      0
#> 4e7d22d793d30bdedbd32a5b73b564b8     10     10      0     0     19
#> e88959b739c8cceb94b94b7887dd1df3      0      0      0     0      0
#> 5403445694339e2ae5ecac06d9f3ff47      0      0      0     0     11
#> ba3943fb075bbf63fb699638292826a1      0      0      0     0      0
#> 67b2d836db093cead5a8286142669613      0      0      0     2      0
#> 6c3a0ea164d3eb098dd2baac5519b38f      9      0      0     0      0
#> c643a776c36bcc00dfb5d790e1ea3dfe      0      0      0     0      0
#> d53cc81a9d57a64d194c876b76bd66e0      0      0      0     0      0
#> 7dc00d40971bb454b39fdd9b2dd9aa25      0      0      0     0      0
#> b4d0cec8543e940be1b5096370fdfb80      0      0      0     0      0
#> 6c54ceb94bc5e8fbd65224cffc60990a      0      0      0     0      0
#> b12ee9ca5a2320957213a114788cf2b7      0      0      0     0      0
#> 5e040c2397f959e7f0dfc6ec2855b812     14      8      0     0     19
#> 3105a594bba2604ae5ed19f5e0495d2e      0     22      0     0      0
#> d75c040bdaa1cca96cd0655c9af7ab28      0      0      0     0      0
#> 9b909edaebabc90b9c1a703626e6bd0c      0      0      0     0      0
#> a6593da81035cf121dbe46346a4ea960      0      0      0     0      0
#> 60e4c88400a750a8122b97476313a67e      0      0      0    12      0
#> 7854a78e3b881cee1c63ea0a50ba29ce      0      0      0     0      0
#> 1796734b56d45bdc997ccc09ab8731d1      0     24      0     0      0
#> 8243cd669ae7327f4326be78429ab05a     16      4      0     0      0
#> ebb93dc32c63402289114859d476cc5b      0      0      0     0      0
#> 8a01af5a531e851e2c50fc01131cda33      0     10      0     0     15
#> 03d2cc3f7b2bc352752d39e3738d0af9      0      0      0     0      0
#> ec378af6e32960b11b5f7391da4ae4fb      0      0      0     0      0
#> b407f2de868c064f65eb5ad8f82a8d0d      0      0      0     3      0
#> 73a95b6ba34fcd260989d848fe15d377      0      0      0     0      0
#> 3abb68202d15781590a2d1f9c5b37e15      0      0      0     0      0
#> 77923ef8e9ed01c9b6a8ac85eb028163      0      0      0     0      0
#> 9981b8fd20f400e052527f84f5db1b66      0      9      0     0      5
#> 597bdec49adf7675912f9fba4e262e70      0      0      0     0      0
#> 0fd3126118b13fb187bc55476a67c95a      0      0      0     0      0
#> 2a0c511ecc2513c42fe624e8e99d07da      0      0      0     0      0
#> d2912fea2dcf2e5b97e56c6a7d24295d      0      0      0     0      0
#> 2976af6aaf4261b6c274fcc3229a4d65      0      3      0     0      0
#> f097c45ce419da576aa111c98a4c7272      0      0      0     0      0
#> 3abc8260d34008f263e232ea25426c3f      0      0      0     0      0
#> d7e5b6432dc178aaeed9c81e0ec57afb      0      0      0     0      0
#> dfe71ce78ddfdbc7b453d17d628365f2     13     27      0     0     28
#> 289b18a86084273ce3b4038bf37e2840      0      0      0     0      0
#> 844a67032c52dc49019bf5dff248bf9c      0      0      0     0      0
#> 38cb3e3ce5ac5aed7dec50247ecfb267      0      0      0     0      0
#> 74e88e71e9929b1ddfcafccfab3b7bba      0      0      0     0      0
#> 66bc0737fc7fa1463ac7189ab7073998      0      0      0     0      0
#> b5a7aa588c9cae6ae8ca2f7d87ac53e8      0      0      0     0     33
#> db0540920457b21b170947aac9137771      0      4      0     0      0
#> ebbaf260510ab16ac6b8c2b9866d05b9      0      0      0     4     15
#> c1f461b67419d4a80a61db0867da9ed6      0     18      0     0      8
#> 0c6eba6452635eecfb865a1049299176      0     11      0     0     17
#> 7ecabb10f29a0c44eb53e9f9852fe22a      0      0      0     0      0
#> a33cb3037e7825920073b4dddaaa41e6      0      0      0     0      0
#> d625812c15c074c2e2b18099f97412c9      0      0      0     0      0
#> 17ef77f12e2a55717e910a5f6a625448      0      0      0     0      0
#> 7baf638fdd3f71166aafd22ecf410348      0      0      0     0      0
#> 9436cda211f60f722865cd23ef75d4ed      0      0      0     0      0
#> c8ab0b20e740eb6b35e68d58f69736e4    109      0      0     0    213
#> 5e43f2d8be71399fc8d6736754cb42ac      0     22      0     0      0
#> 19ad90ee1081a003406b02179ae63709      8      0      0     0      0
#> 0a25cdbd0660e68970c563c43a18f54a      0      0      0     0     19
#> 0aa6b210e902d8ab9d6f267407258ebb     13      0      0     0     46
#> a864c36d0d898cc5abd94d6d6f0bd5be      0      0      0     0      0
#> dea0dcb305ff0f27e3c36c7133cf933d      0      0      0     0      0
#> 463977baa69c7724fdc9407ad162714c      0      0      0     0      0
#> fcb72326c29342bcfb052953466eafc6      0     18      0    18     24
#> 5a554353f9ed728f0ab0631f29c8be2e      0      0      0     0      0
#> a9dd27b48e443415a743bf9902f515c4      0      0      0     0      0
#> 77a40bcf979c0099c52f19a3f1fe5217      0      0      0     0      0
#> b4007ac525a9b409d70426e24421cff0      0     14      0     0      0
#> 406ba18a175c4c79efb25f48a6beafe4      0      0      0     0      0
#> fe1c166755bc722cefd2e2a3b2b70672      0      0      0     0      0
#> 8afd274849ac07c0ca66f8f5f143df3c     27      0      0     0     21
#> 197301ed9d8728b5c061b99a17973161      0     58      0     0     27
#> 14779f9dc3505598ad3d799e4ee676c4      0      0      0     0      0
#> 4085c246640eee043d523ce2d799023a     24     43      0    21      0
#> bba09a403e67c1381de01c86511f1057      0      0      0    18      0
#> 392c59fee4e6b91d7654f2d3e3c65be8     12      0      0     0     16
#> d2b05ae73ad5ed296d7add9e2bb2cce4      0      8      0     0     15
#> 8306a9473dabfd5918ea869dbfcd1291      0      0      0     0      0
#> 06064c372bb624503c2905d13f7691f5     90    142     48    27     92
#> 3265692341e42dbd5cf4ed1938b906e3      0      0      0     0      0
#> da96a0797b9881a123795838b9f01140      0     23      0     0      0
#> 0bdddf10e577657becec88e42001fd16      0      0      0     0     24
#> 168999a950ebe17e6d2ea1f90e557d91      0      0      0     0     16
#> a266ffe9d38f626d106fa615da4618ef      0      0      0     0      0
#> 7706405d5bce7715f01d8b54188def42      0     23      0     0      0
#> 4b0cfc8899b466bd9992f7b22a5a06be      0      0      0     0      0
#> b81b2b4a7ca9eb2c23b170648dfaf22d      0     21      0     0      0
#> 75cf40872201f150cb5b191676f64b78      0     76      0     0     52
#> c8e5794c58de9326084a62d3282d257e      0      0      0     0      0
#> b0706a571f169c0bac7b4864052fcaa7      0      0      0     0      0
#> 0388f969e3527b2938acc676d6757bdb      0      0      0     0      0
#> 2dfe751390aedfbd9f62dbb77d1f10b4      0      0      0     0      0
#> 2548dc9697e1f9afe04a9012f40e2f0a      0      0      0     0      0
#> 2f019615d5141a4ab00b1aab98882a33      0     17      0     0     12
#> 629e4f0cf1826904ff74907b2507f58b      0      6      0     0      0
#> 774861f348cec8451102730c21a790b0      0     12      0     0      0
#> 11c0d6de28ab026a9ae17abe4efdd6bf      8      0      0     0      0
#> c11381dda0841c0d08cbfcef3a760aa3      0      0      0     0      0
#> 3c9619d4f9db4ac6ad7151a0c9d86ad3      0      0      0     0      0
#> 7e2c72204a8c45b1b0ac39e0e6047098      0    133      0     0      0
#> 8254e0f9f07709b62dd8483f3b25fad0     16     65      0     0     28
#> e78135f2886d427e678d12c446bf070e      0      0      0     0      0
#> abec938a8d53e7471f3f476d50694be9      0      0      0     0      0
#> 4963567d2979bb761291aa98e0b8595b      0     26      0     0      0
#> d4de9d1453102564f35262faa3952266      0      0      0     0      0
#> efa4a78c4b993e9d115841c8a02c795e      0      0      0     0      0
#> fcacd7424bf3543322ea3e9842d14e86     11     14      0     0     28
#> f941e15ea7f7560191d57abc4d2bedd8     13      0      0     5      6
#> ac294f5242de3c9b26e75502c48f6383      0      0      0     0      0
#> 7b6a8e3f8bfa1bcd1b5fc55db048d005      0      0      0     0      0
#> c9e776ab599fe75db286675b698c9067     30     24      0    34      0
#> ebe5cb945f6d6d849c07424591f80753      0     12      0     0      0
#> f116387f77802409103a559e4205df06     11      0      0     0      0
#> 93f9939b261cddc9a57071b42771c1d4      0      0      0     0      0
#> 971a5f48b3418449ca7cf3e190a5e831      0      0      0     0      0
#> 730de992ca06a95e7193a93146a97cf2     41      0      0     0      0
#> 8fa8b4ac69784c26901e036d13ab2332      0      0      0     0      0
#> c363c9efaa85023b488b732b97e47270      0      0      0     0      0
#> 2e2181cc893ca1e1e677e888bf8af441      0      0      0     0      0
#> 86fa662da6869634e7f90f85e9ee4933      0      0      0     0      0
#> 9a02177f2c235a6eb62984d5aab151ff      0      0      0     0      0
#> d6854be498e22b65da3c00572ee178f0      0      0      0     0      0
#> 015f364f150b230450da0a123e7bb8ab      0      0      0     0      0
#> af2c1b22758cd6653e8e841497d096e7      0      0      0     0      0
#> e7f2fab8c0e7c365768a43ef71406605      0     59      0    10      0
#> fb57cf192bce1876471f2f6181bf4b55      0      0      0     0      0
#> 6c740ce69c9f45bea684dc9e36054006      0      0      0     0      0
#> 3b27c7bf279535d8e39664d8dc1d0daf      0      0      0    31      0
#> ceea3c29ad012c2918fb6aa21a44e8cb      0     45     10     0      0
#> fa99ecdb27c1cc44fd866ee1da250c4c      0      0      0     0      0
#> 387cc83573db0e444e81c89a8a4f8f9c     32     67      0     0     42
#> c7a8d1af5f7f472671f700775de3e0a2     18      0      0     0     11
#> 98a6a84550e07ea252cfd8a0f51042d4      0     70      0     0      0
#> 08b809a68126c38db786cafe9710b013      0     15      0     0      0
#> 8db4cc0a5d2f03fce8efd2b8374ff359      0      0      0     0      0
#> 3e7c7d38423c93b43a1d9865cf0bc4a6      0     43      8     0     22
#> 6f02050760c030ed8cd1ed836697fa6b      0     12      0     0     26
#> 8fc1865cfcae32f72c1b1570299e2ef3      0      0      0     0      0
#> 38cc138fac95b006929361093d840a00      0     20      0     0      0
#> f5f8eb63e22362509e2c7c78271305f2      0      0      0     0      0
#> 0bb238aaee3ccb7dc8dfa86dbc1ba9f7      0      0      0     0      0
#> ca4546caf07edf5200cf79cea04aebcd      0      0      0     0      0
#> 893d12123eb7f623489a4f025a15ae7f      0      0      0     0      0
#> 192c5a817f48368339e519e0a0a907da      0      0      0     0      0
#> 3467185ca92b7bec2523214605adabf0      0      0      0     0      0
#> a2aae80d8fb922fe1074d243b59ce237     20     42      0     0      0
#> f11233f16d4177d38ab1c0bfc649e056      0      0      0     0      0
#> b219ca7793ab9f5e441f2a68bc76152c      0      0      0     0      0
#> e6e1848b51e705e55b3805703f088311      0    124      0     0      0
#> bd074b9bf365d11468bd76268d63b15b     20     35      0     0     19
#> f381b492d6651c52fd09ce9351e7eb2d      0      0      0     0      0
#> d9b31e1c30facc24fa778b1d65f5c457     90      0      0     0     83
#> 14a67c8c8cc31d203eac8109443fab5f      0      0      0     0      0
#> 026ad0ea53fb202ebb770acf65705d6e      0     20      0    23      0
#> 19ba9fea04969f975d6dc8d3da20d5e6      0      0      0     0      0
#> 0fc4426ed7dccc48ce27cfff8083433d    117    295     56   357    161
#> 45ad90e179842705fa984962f131fd61      0      0      0     0      0
#> d284afab226d22fd0596221d17db23d3      0     21      0     0      0
#> 38194655ba17dbfd7fcc1ed42a57c169     11     34      0     0     42
#> e30be5268bbebe597424bdf171cc1607      0      0      0     0      0
#> cb7ccdb4956134b452a3cb78d1751c57      0      0      0     0      0
#> b585436677ac88dd4ac6239655f2a268      0      0      0     0      0
#> f5d770eff5c9b9b0c859db57d175eb34      0      0      0     0      0
#> 1b45ffbc88c788dcee2239e68901cf05      0      0      0     0      0
#> 3fd951cf7cc3a0740a1f20f8a9c6d2c2      0      5      0     0      0
#> 2027d3e90ccc4003f1954322e99db6a9      0      9      0     0      0
#> 6157cb92771830955413f808e83fc420      0     27      0     0      0
#> 6349506cac32056a214b04239bdc3a33     40     66      0     0    123
#> 2254a8d37cf537b1e35097d70b6ac91b     97     43     24    18    512
#> 38f3bc8ef8836574de0fcbeafd348b65      0      0      0     0      0
#> fec23c70ae0c70ecceea1487cd978ee5      0      0      0     0      0
#> 503d123a6c78c6fb280a6472bd61b51e      0      0      0     0     10
#> 9f33de8570c639e727465f1ee4d750f6      0      0      0     0      0
#> d68e1899fa65515135d9bd8d356b7e3d      0      0      0     0      0
#> e15b88ccac7a4bef64215e915f6e0d89      0      0      0     0      7
#> 767042b4b19d8bfaaa426dad7dbce378      0      0      0     0      0
#> 4a3c4abee502d6185abd750522566bf3      0      0      0     0      5
#> c254aa5c00b7f27d756f7bbe96b3d4ad      0      0      0     0      0
#> fe8d8af791eb6e5d4335bb7604c74e9d      0      0      0     0      0
#> d81800d7e31b69877ce1eab9858b8f3a      0      0      0     0     10
#> c55c7cb80c5700ff9aded76192543b5a      0      0      0     0      0
#> 91dee7b3b90bbf2e53b2472b9bbde0e9      0      0      0     0      0
#> 6b28df61cb5dee74b76b1a13aaa0027d      0      0      0     0      0
#> 9a830fa26205ad91a02b5bcaebb0221a      0      0      0     0      0
#> a54b8142d2a09d058b81c377e8e286b8      0      0      0     0      0
#> 69ef3a03a4d733dbffd99fdaa1d2e62a      0      0      0     0      0
#> 083cbb131dd305a88b7d338bc371c35f      0     23      0     0      0
#> ef1a92caffb4d7cae286fdc77b0a93a6      0      0      0     0      0
#> 2cddd5dca166486f50a0828b28001e1e      0      0      0     0      0
#> bb899d1dcfd89cb200a76cb8951904f8     13      0      0     0      0
#> 7817f6851e9df03a3d992cf915310dba      0      8      0     0      0
#> 54eaf6571de743a86ccb531bee8aad7e      0      0      0     0      0
#> 
#> $long_format
#> # A tibble: 29,072 × 3
#>    taxonomy                                               SAMPLEID Counts
#>    <chr>                                                  <chr>     <int>
#>  1 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 3.16RI        0
#>  2 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 2.14RO        0
#>  3 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 3.10RO        0
#>  4 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 3.10RI        0
#>  5 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 4.16RI       25
#>  6 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 1.8RO         0
#>  7 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 4.1RI        35
#>  8 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 4.1RO         0
#>  9 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 2.16RI        0
#> 10 d__Bacteria; p__Verrucomicrobiota; c__Verrucomicrobiae 1.2RI         0
#> # ℹ 29,062 more rows
#> 
```
