library(readr)
library(dplyr)

genomes <- read_csv(snakemake@input[['genome_taxonomy']])
phix <- data.frame(ident = 'GCF_000819615.1',
                   superkingdom = "d__Virus",
                   phylum = 'p__Phixviricota',
                   class = 'c__Malgrandaviricetes',
                   order = 'o__Petitvirales',
                   family = 'f__Microviridae',
                   genus = 'g__Sinsheimervirus',
                   species = 's__PhiX')
contams <- read_csv(snakemake@input[['picklist']]) %>%
  mutate(ident = gsub(" .*", "", name))

gtdb <- read_csv(snakemake@input[['gtdb_taxonomy']]) %>%
  filter(ident %in% contams$ident)

taxonomy <- bind_rows(gtdb, genomes, phix)

write_csv(taxonomy, snakemake@output[['taxonomy']])
