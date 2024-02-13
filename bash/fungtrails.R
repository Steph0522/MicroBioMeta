tabla_guilds<- read.delim("/Documents and Settings/HP/Downloads/table.taxa.guilds.txt")
tab<-tabla_guilds %>% filter(!trophicMode=="na")

taxa<- read.delim("/Documents and Settings/HP/Downloads/table.taxa.txt")

taxa_fun<- taxa %>% inner_join(data, by = "Genus") %>% dplyr::select(-n)

trails<- read.delim("/Documents and Settings/HP/Downloads/fungtrails.txt")
dt<- trails %>% filter(!primary_lifestyle=="#N/D" ) %>% filter(!primary_lifestyle==0 )

write_tsv(dt,"/Documents and Settings/HP/Downloads/fungtrails_matched.txt")

install.packages("devtools")
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("traitecoevo/fungaltraits")
library(fungaltraits)
a<-fungal_traits()
head(fungal_traits())
